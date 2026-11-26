"""Three-stage diagrams: structure, reaction progress, typed electron curves."""
import copy,math,textwrap,time
from pathlib import Path
from xml.etree import ElementTree as ET

from ..benchmark.umr_bridge import digest,parse_document
from ..umr.document import run_document
from ..umr.exporting import write_json
from .geometry import add,collisions,distance,head,mul,normal,route,samples,Spatial,sub,unit
from .palette import PURPLE, TEAL, TEAL_TEXT
from .models import RenderConfig
from .materials import display_plan,phase_inventory
from .ports import requests_for
from .electron_layout import place_electrons
from .structure import place_layout,structure
from .svg import element,label,svg_to_png,text_box,write_svg
from .language import event_name,conditions_text


def arrow_line(group,a,b,*,resonance=False):
    color=PURPLE if resonance else TEAL;u=unit(sub(b,a));v=normal(u)
    group.append(element('line',x1=a[0],y1=a[1],x2=b[0],y2=b[1],stroke=color,stroke_width=2.5))
    for p,direction in [(b,u)]+([(a,mul(u,-1))] if resonance else []):
        one=add(add(p,mul(direction,-11)),mul(v,5));two=add(add(p,mul(direction,-11)),mul(v,-5))
        group.append(element('path',d=f'M {p[0]} {p[1]} L {one[0]} {one[1]} L {two[0]} {two[1]} Z',fill=color))


def marker_layer(view,arrows,coords,parent):
    objects=[];atoms={a.atom_id:a for a in view.atoms};selected={}
    for r in arrows:
        if r['source']['type'] in {'lone_pair','radical','metal_owned_pair'}:
            selected[(r['source']['atom_ids'][0],r['source'].get('slot',0),r['source']['type'])]=r['source']['position']
    donors={r['source']['atom_ids'][0] for r in arrows if r['source']['type'] in {'lone_pair','metal_owned_pair'}}
    for i,a in atoms.items():
        if a.radical_electrons is None:continue
        kinds=[('radical',a.radical_electrons)]
        if i in donors:
            kind='metal_owned_pair' if a.element in {'Pd','Zn'} else 'lone_pair'
            count=sum(1 for x in selected if x[0]==i and x[2]==kind) if kind=='metal_owned_pair' else max(0,((a.nonbonding_electrons or 0)-a.radical_electrons)//2)
            kinds.append((kind,count))
        for kind,count in kinds:
            for slot in range(count):
                angle=slot*math.pi/2-math.pi/2
                pos=selected.get((i,slot,kind),add(coords[i],(22*math.cos(angle),22*math.sin(angle))))
                tangent=normal(unit(sub(pos,coords[i])));points=[pos] if kind=='radical' else [add(pos,mul(tangent,-2.8)),add(pos,mul(tangent,2.8))]
                for p in points:parent.append(element('circle',cx=p[0],cy=p[1],r=2,fill=PURPLE if kind=='radical' else TEAL_TEXT,data_atom_id=i))
                objects.append({'id':f'electron.{kind}.{i}.{slot}','kind':'electron','atom_ids':[i],'box':[pos[0]-6,pos[1]-6,pos[0]+6,pos[1]+6]})
    return objects


def draw_panel(source,target,actions,event_id,title,conditions,output,config,*,resonance=False,inactive=(),state_only=False,required_hydrogens=(),spin_transition_request=False,spin_transition_executed=False,declared_pair_donors=(),
               event=None,released_fragments=(),definitions=()):
    started=time.perf_counter();frozen=[digest(v.model_dump(mode='json')) for v in (source,target) if v]
    plan,touched,pairs=display_plan(source,target,actions,event,released_fragments,inactive,config,state_only)
    visible={side:{i for component in rows if component['displayed'] for i in component['atom_ids']}
             for side,rows in plan.items()}
    attempts=[];chosen=None
    for retry in range(config.max_relayouts+1):
        options=dict(touched=touched,expansion=1+retry*.35,inactive=inactive,required_hydrogens=required_hydrogens,
                     contact_pairs=pairs,definitions=definitions)
        left=structure(source,config,visible_ids=visible['input'],**options)
        right=structure(target,config,visible_ids=visible['output'],align=left['mol_coords'],**options) if target else None
        metal_notes=[(side,v,m) for side,v in [('INPUT',source),('OUTPUT',target)] if v for m in getattr(v,'metal_centers',[]) if m.atom_id in visible[side.lower()]]
        metal_footer=15*len(metal_notes)
        lw=left['width']+130;rw=(right['width'] if right else 520)+100;width=lw+rw+200;height=max(left['height'],right['height'] if right else 280)+300+metal_footer
        if state_only:width=lw+80
        if width*height>config.max_panel_pixels*2:raise ValueError('Event canvas exceeds pixel budget')
        root=element('svg',width=width,height=height,viewBox=f'0 0 {width} {height}')
        root.append(element('rect',width='100%',height='100%',fill='white'))
        groups={name:element('g',id=name,**({'class':name} if name else {})) for name in ['molecular-structures','reaction-flow','electron-arrows','review-annotations','atom-ids','electron-markers']}
        for g in groups.values():root.append(g)
        groups['reaction-flow'].set('data-activation',str(conditions.get('activation','unknown')))
        if not config.show_atom_ids:groups['atom-ids'].set('style','display:none')
        label(groups['review-annotations'],title,24,28,size=18,fill=PURPLE)
        if hasattr(source,'electronic_state'):
            if spin_transition_request or config.spin_display=='details':
                def spin_text(v):
                    return ' / '.join(dict.fromkeys(x for x in (v.electronic_state.character,v.electronic_state.label)
                                                     if x and x!='unknown'))
                before=spin_text(source);after=spin_text(target) if target else ''
                if before or after:
                    note='电子态  '+((before or '待定')+' → '+(after or '待定') if target and before!=after else before or after)
                    label(groups['review-annotations'],note,24,51,size=12,fill='#713688')
            if spin_transition_request:
                groups['reaction-flow'].append(element('line',x1=lw-20,x2=lw+150,y1=290,y2=290,stroke='#713688',stroke_width=2,stroke_dasharray='7 5'))
                label(groups['review-annotations'],'电子态声明已衔接 · 未计算自旋动力学' if spin_transition_executed else '态间转换请求 · 未执行，不作为共振',lw-70,265,size=14,fill='#713688')
        label(groups['atom-ids'],f'{event_id} | {source.state_id}/{source.view_id}'+(f' → {target.state_id}/{target.view_id}' if target else ''),24,52,size=12)
        condition_text=conditions_text(conditions)
        for k,line in enumerate(textwrap.wrap(condition_text,width=150)[:3]):label(groups['reaction-flow'],line,24,74+16*k,size=12,fill='#444')
        coords,obstacles,atoms=place_layout(left,groups['molecular-structures'],(58,148),source,'input',event_id,groups['atom-ids'],config.show_atom_ids)
        if right:
            _,obs,rs=place_layout(right,groups['molecular-structures'],(lw+200,148),target,'output',event_id,groups['atom-ids'],config.show_atom_ids);obstacles+=obs;atoms+=rs
        # Stage 2: explicit molecule separators, progress arrow and conditions.
        forming_ports=[]
        for action in actions:
            raw_action=action.model_dump(mode='json') if hasattr(action,'model_dump') else action
            pair=(raw_action['b'],raw_action['c']) if raw_action.get('op')=='SHIFT' else (raw_action['destination']['a'],raw_action['destination']['b']) if raw_action.get('op')=='FISHHOOK' and raw_action['destination']['kind']=='bond' else None
            if pair and all(i in coords for i in pair):forming_ports.append(mul(add(coords[pair[0]],coords[pair[1]]),.5))
        for side,layout,offset in [('input',left,(58,148))]+([('output',right,(lw+200,148))] if right else []):
            boxes=layout['components']
            for one,two in zip(boxes,boxes[1:]):
                a,b=one['box'],two['box']
                if abs((a[1]+a[3])/2-(b[1]+b[3])/2)<70 and b[0]>a[2]:
                    x=offset[0]+(a[2]+b[0])/2;y=offset[1]+(a[1]+a[3]+b[1]+b[3])/4+6
                    if side=='input' and any(distance((x,y),port)<36 for port in forming_ports):
                        y=offset[1]+max(a[3],b[3])+45
                    label(groups['reaction-flow'],'+',x,y,size=22,fill='#444',data_separator_side=side)
            for box in boxes:
                if box['inactive']:
                    b=box['box'];label(groups['review-annotations'],'NOT ADDED / inventory',offset[0]+b[0],offset[1]+b[1]-8,size=11,fill='#8a6500')
        if not state_only and not spin_transition_request:
            arrow_line(groups['reaction-flow'],(lw+20,height/2),(lw+165,height/2),resonance=resonance)
            label(groups['reaction-flow'],'共振' if resonance else '反应',lw+90,height/2-16,size=12,text_anchor='middle',fill=PURPLE if resonance else '#333')
        if not target and not state_only:label(groups['review-annotations'],'审查不通过：未得到可确认的输出结构',lw+220,height/2,size=13,fill='#a25b3f')
        label(groups['review-annotations'],'完整箭头：双电子 · 鱼钩箭头：单电子 · 紫色双向箭头：共振',24,height-28,size=12,fill='#75848a')
        for i,(side,v,m) in enumerate(metal_notes):
            show=lambda x:'unknown' if x is None else x
            values={'oxidation_state':m.oxidation_state,'d_electrons':m.d_electrons,'valence_electrons':m.valence_electrons}
            if hasattr(v,'element_rule_version'):
                from ..umr.v05.coordination import checks
                if v.element_rule_version=='umr-elements/0.7':
                    from ..umr.v07.elements import coordination as checks
                if v.element_rule_version=='umr-elements/0.8':
                    from ..umr.v08.elements import coordination as checks
                if hasattr(v,'material_boundaries'):
                    from ..umr.v09.coordination import checks
                    from ..umr.v07.ledger import project_symbols
                    from ..umr.chemistry import FormalError
                    try:values.update(next(c for c in checks(project_symbols(v))[1] if c['metal_atom_id']==m.atom_id))
                    except FormalError:values={key:None for key in values}
                else:values.update(next(c for c in checks(v)[1] if c['metal_atom_id']==m.atom_id))
            from .coordination_layout import GEOMETRIES
            graph=left if side=='INPUT' else right
            projection=next((r for r in graph['coordination'] if r['metal_atom_id']==m.atom_id),{})
            shape=GEOMETRIES.get(projection.get('geometry'),('',))[0]
            description=''
            if shape:
                description=' · '+shape
                if projection['status']=='schematic':description+='示意（构型未定）'
                elif not projection['sites']:description+='（投影待核对）'
                elif projection.get('site_assignment')=='schematic':description+='（位点示意）'
            label(groups['review-annotations'],f'{side} metal atom {m.atom_id}: formal q={next(a.formal_charge for a in v.atoms if a.atom_id==m.atom_id)}; OS={show(values["oxidation_state"])}; d={show(values["d_electrons"])}; VE={show(values["valence_electrons"])} (ionic)'+description,24,height-70-metal_footer+i*15,size=11,fill=TEAL_TEXT)
        # Text/legend/main-arrow obstacles use actual canvas positions.
        for j,node in enumerate(groups['reaction-flow']):
            if node.tag.endswith('text'):
                x=float(node.get('x'));y=float(node.get('y'));obstacles.append({'id':f'flow.text.{j}','kind':'text','atom_ids':[],
                    'box':text_box(node.text or '',x,y,float(node.get('font-size',12)))})
        input_obstacles=[o for o in obstacles if not o['id'].startswith('output.')]
        reqs,unsupported=requests_for(source,actions,coords,event_id,inactive,input_obstacles,declared_pair_donors)
        electron_objects,electron_display=place_electrons(left['display_view'],reqs,coords,groups['electron-markers'],input_obstacles,config)
        routed,diagnostics=route(reqs,[o for o in obstacles+electron_objects if not o['id'].startswith('output.')],(12,127,lw-8,height-70),config)
        if target:place_electrons(right['display_view'],[],{a['atom_id']:a['canvas_xy'] for a in atoms if a['side']=='output'},groups['electron-markers'],obstacles,config,side='output')
        final_spatial=Spatial(obstacles+electron_objects)
        final_hits=[]
        for req in routed:
            hits,_=collisions(req['control_points'],req,final_spatial,config.collision_padding)
            final_hits += [{'arrow_id':req['arrow_id'],**h} for h in hits]
        diagnostics['collisions']=final_hits
        if final_hits:diagnostics['status']='layout_warning'
        hydrogen_overlaps=[{'side':side,'display_id':h['display_id'],'overlaps':h['layout_overlap_count']}
            for side,layout in [('input',left),('output',right)] if layout
            for h in layout['hydrogen_display']['display_hydrogens'] if h['layout_overlap_count']]
        diagnostics['hydrogen_layout_overlaps']=hydrogen_overlaps
        if hydrogen_overlaps:diagnostics['status']='layout_warning'
        if unsupported:diagnostics['status']='render_unsupported'
        for req in routed:
            c=req['control_points'];color=TEAL_TEXT if req['electrons']==2 else PURPLE
            group=element('g',id=req['arrow_id'],data_action_id=req['action_id'],data_event_id=event_id,data_electrons=req['electrons'])
            path=element('path',d=f'M {c[0][0]:.2f} {c[0][1]:.2f} C {c[1][0]:.2f} {c[1][1]:.2f}, {c[2][0]:.2f} {c[2][1]:.2f}, {c[3][0]:.2f} {c[3][1]:.2f}',stroke=color,stroke_width=2.3,fill='none',stroke_linecap='round',stroke_linejoin='round')
            group.append(path)
            tips=head(c,req['electrons']);end=c[-1]
            if req['electrons']==2:
                p,q=tips[0][1],tips[1][1];group.append(element('path',d=f'M {end[0]} {end[1]} L {p[0]} {p[1]} L {q[0]} {q[1]} Z',fill=color))
            else:
                wing=tips[0][1];u=unit(sub(c[-1],c[-2]));control=add(add(end,mul(u,-4)),mul(normal(u),6))
                group.append(element('path',d=f'M {end[0]} {end[1]} Q {control[0]} {control[1]}, {wing[0]} {wing[1]}',stroke=color,stroke_width=2.5,fill='none',stroke_linecap='round'))
            mid=samples(c)[14];label(group,str(req['action_index']+1),mid[0]+5,mid[1]-7,size=11,fill=color,**{'class':'arrow-ids'})
            groups['electron-arrows'].append(group)
        if diagnostics['status']!='rendered':label(groups['review-annotations'],'绘图待调整：保留全部动作与问题定位',24,height-8,size=11,fill='#a04d00')
        if config.diagnostics:
            for o in obstacles+electron_objects:
                b=o['box'];groups['review-annotations'].append(element('rect',x=b[0],y=b[1],width=b[2]-b[0],height=b[3]-b[1],fill='none',stroke='#cc9999',stroke_width=.5))
        diagnostics.update({'relayout_attempt':retry,'unsupported_actions':unsupported,'actions_total':len(actions),'arrows_drawn':len(routed),
            'depiction_version':'umr-depiction/1.2','presentation':{
                s:{'abbreviations':l['labels'],'scaffolds':l['scaffolds'],'partner_layout':l['partner_layout'],
                   'components':l['components'],'coordination':l['coordination'],'chelate_links':l['chelate_links']}
                for s,l in [('input',left),('output',right)] if l},
            'hydrogen_display':{'input':left['hydrogen_display'],'output':right['hydrogen_display'] if right else None},
            'electron_display':electron_display,'undisplayed_lone_pairs_retained_in_ledger':True,
            'semantic_correspondence_complete':len(actions)==len(routed),'atom_coordinate_count':len(atoms),
            'bounds_input_only':[12,127,lw-8,height-70],
            'pipeline':['1_RDKit_structure_and_layout','2_main_reaction_separators_conditions','3_typed_electron_overlay'],
            'chemical_signature_unchanged':left['chemical_signature_unchanged'] and (right['chemical_signature_unchanged'] if right else True)})
        result=(root,routed,atoms,diagnostics,left)
        attempts.append({'attempt':retry,'status':diagnostics['status'],'collisions':len(final_hits),'arrow_crossings':len(diagnostics.get('arrow_crossings',[]))})
        quality=(len(actions)-len(routed))*10000+len(final_hits)+len(diagnostics.get('arrow_crossings',[]))*3+len(unsupported)*100+sum(h['overlaps'] for h in hydrogen_overlaps)
        if chosen is None or quality<chosen[0]:chosen=(quality,result)
        if diagnostics['status'] in {'rendered','render_unsupported'}:break
    root,arrows,atoms,diagnostics,_=chosen[1];diagnostics['relayouts']=attempts;diagnostics['elapsed_seconds']=time.perf_counter()-started
    if frozen!=[digest(v.model_dump(mode='json')) for v in (source,target) if v]:raise ValueError('Panel generation changed frozen UMR')
    output.parent.mkdir(parents=True,exist_ok=True)
    from .figure_keys import append_keys
    diagnostics['abbreviation_key']=append_keys(root,[f for side in diagnostics['presentation'].values() for f in side['abbreviations']])
    from .local_rules import append_rules
    local_ids=append_rules(root,[v.view_id for v in (source,target) if v], [event] if event else [])
    if local_ids:diagnostics['local_rule_ids']=local_ids
    write_svg(root,output.with_suffix('.svg'));raster=svg_to_png(output.with_suffix('.svg'),output.with_suffix('.png'))
    inventory_file=output.with_suffix('.inventory.json')
    write_json(inventory_file,{'schema_version':'umr-drawing-inventory/1','event_id':event_id,
        'sides':{s:{'components':plan[s],'view':v.model_dump(mode='json'),'view_sha256':digest(v.model_dump(mode='json')),
                    'abbreviations':diagnostics['presentation'][s]['abbreviations']}
                 for s,v in [('input',source),('output',target)] if v}})
    metadata={'schema_version':'umr-arrow-layout/0.1','event_id':event_id,'source_view_id':source.view_id,'target_view_id':target.view_id if target else None,
        'is_resonance':resonance,'inventory':inventory_file.name,'source_sha256':frozen[0],'target_sha256':frozen[1] if target else None,
        'raw_actions':[a.model_dump(mode='json') if hasattr(a,'model_dump') else a for a in actions], 'arrows':arrows,'atoms':atoms,'rasterization':raster,'layout':diagnostics}
    write_json(output.with_suffix('.coordinates.json'),metadata)
    write_json(output.with_suffix('.layout_diagnostics.json'),diagnostics)
    return {'event_id':event_id,'kind':'spin_transition_request' if spin_transition_request else 'state_only' if state_only else 'resonance' if resonance else 'reaction','status':diagnostics['status'],
        'svg':output.with_suffix('.svg').name,'png':output.with_suffix('.png').name,'metadata':output.with_suffix('.coordinates.json').name,
        'diagnostics':output.with_suffix('.layout_diagnostics.json').name,'inventory':inventory_file.name,'actions':len(actions),'arrows':len(arrows),'view_id':source.view_id}


