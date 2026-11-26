"""A folded, complete chain of actual input views and the final executed view."""
import copy,json,math
from xml.etree import ElementTree as ET
from .svg import element,label,write_svg,svg_to_png
from .language import event_name


def draw_chain(directory,panels,document):
    from ..umr.v09.local_rules import use_rules
    with use_rules(document):
        return _draw_chain(directory,panels,document)

def _draw_chain(directory,panels,document):
    events=[p for p in panels if p.get('kind')=='reaction' and p.get('metadata') and p.get('executed')]
    if not events:return None
    frames=[];transitions=[];handoffs=[]
    event_order={e.event_id:k for k,e in enumerate(document.events)}
    titles=document.metadata.get('event_labels',{})
    if not isinstance(titles,dict):
        titles={}
    arrow_titles={e.event_id:titles.get(e.event_id) or e.conditions.get('step_title') for e in document.events}
    views={v.view_id:v for v in document.views}
    previous=None
    for index,p in enumerate(events):
        meta=json.loads((directory/p['metadata']).read_text(encoding='utf-8'))
        event=next(e for e in document.events if e.event_id==p['event_id'])
        if previous:
            pp,pm,pe=previous
            before=views.get(pm['target_view_id']);after=views.get(meta['source_view_id'])
            contiguous=event_order[event.event_id]==event_order[pe.event_id]+1
            same=(contiguous and before is not None and after is not None and before.view_id==after.view_id
                  and pm.get('target_sha256')==meta.get('source_sha256'))
            witness=next((t for t in document.resonance_transforms if before and after and
                before.state_id==after.state_id==t.state_id and t.source_view_id==before.view_id and
                t.target_view_id==after.view_id),None)
            if witness and contiguous and not same:
                # Use only an executed/validated view switch. Existence of an
                # unverified resonance assertion alone cannot connect two chains.
                try:
                    import importlib
                    version=document.schema_version.split('/')[-1].replace('.','')
                    if version=='10':version='09'
                    mod='mechanism_workbench.umr.'+('' if version=='01' else 'v'+version+'.')+'resonance'
                    importlib.import_module(mod).validate_resonance(before,after,witness)
                except (ValueError,KeyError,ImportError,AttributeError):
                    witness=None
            kind='shared_state' if same else 'resonance' if witness and contiguous else 'gap'
            if same:
                transitions.append('reaction')
            else:
                if before is not None:
                    frames.append((pp,pm,'output','前一步产物'))
                    transitions.append('reaction')
                transitions.append(kind)
            def participating(panel_meta,side):
                if panel_meta.get('inventory'):
                    pool=json.loads((directory/panel_meta['inventory']).read_text(encoding='utf-8'))
                    return {i for c in pool['sides'].get(side,{}).get('components',[]) if c['displayed'] for i in c['atom_ids']}
                return {a['atom_id'] for a in panel_meta['atoms'] if a['side']==side}
            outgoing=participating(pm,'output');incoming=participating(meta,'input')
            handoffs.append({'from_event':pe.event_id,'to_event':event.event_id,'product_view':pm['target_view_id'],
                'reactant_view':meta['source_view_id'],'kind':kind,'state_id':after.state_id if after else None,
                'atom_ids':sorted(a.atom_id for a in after.atoms) if same else [],
                'continued_atom_ids':sorted(outgoing & incoming) if same else [],
                'from_pool_atom_ids':sorted(incoming-outgoing) if same else [],
                'to_pool_atom_ids':sorted(outgoing-incoming) if same else [],
                'resonance_transform':witness.transform_id if kind=='resonance' else None})
        frames.append((p,meta,'input',f'{index+1:02d}  '+event_name(event)))
        previous=p,meta,event
    if any(a['side']=='output' for a in frames[-1][1]['atoms']):
        p,meta,_,_=frames[-1];frames.append((p,meta,'output','最终状态'));transitions.append('reaction')
    tw,th,gap,pad=580,310,100,26;columns=4 if len(frames)>4 else 2
    if any(arrow_titles.values()):
        gap=max(gap, min(290,24+max(len(str(t)) for t in arrow_titles.values() if t)*10))
    rows=math.ceil(len(frames)/columns);width=columns*tw+(columns-1)*gap+2*pad;height=rows*(th+88)+40
    root=element('svg',width=width,height=height,viewBox=f'0 0 {width} {height}')
    root.append(element('rect',width='100%',height='100%',fill='white'))
    defs=element('defs');root.append(defs)
    positions=[]
    folds=[]
    for index,(panel,meta,side,title) in enumerate(frames):
        row=index//columns;col=index%columns if row%2==0 else columns-1-index%columns;x=pad+col*(tw+gap);y=28+row*(th+88)
        positions.append((x,y))
        label(root,title,x+18,y+28,size=18,fill='#660874')
        pts=[a['canvas_xy'] for a in meta['atoms'] if a['side']==side and a.get('visible',True)]
        parents={a['atom_id']:a for a in meta['atoms'] if a['side']==side}
        presentation=meta['layout'].get('presentation',{}).get(side,{})
        folds.extend(presentation.get('abbreviations',[]))
        for link in presentation.get('chelate_links',[]):
            anchor=parents[link['backbone_atom_ids'][0]]
            dx,dy=anchor['transform']['translate']
            x0,y0,x1,y1=link['box']
            pts.extend([(x0+dx,y0+dy),(x1+dx,y1+dy)])
        for h in meta['layout'].get('hydrogen_display',{}).get(side,{}).get('display_hydrogens',[]):
            offset=parents[h['parent_atom_id']]['transform']['translate']
            pts.append([h['local_xy'][k]+offset[k] for k in (0,1)])
        if side=='input':
            pts+=[xy for a in meta['arrows'] for xy in a['control_points']]
        x0,y0=min(p[0] for p in pts)-63,min(p[1] for p in pts)-64
        x1,y1=max(p[0] for p in pts)+63,max(p[1] for p in pts)+64
        original=ET.parse(directory/panel['svg']).getroot()
        scale=min((tw-24)/(x1-x0),(th-12)/(y1-y0))
        tx=x+12+((tw-24)-(x1-x0)*scale)/2;ty=y+45+((th-12)-(y1-y0)*scale)/2
        clip=element('clipPath',id=f'chain-clip-{index}',clipPathUnits='userSpaceOnUse');clip.append(element('rect',x=x+12,y=y+45,width=tw-24,height=th-12));defs.append(clip)
        view_id=meta['source_view_id'] if side=='input' else meta['target_view_id']
        handoff=next((h for h in handoffs if side=='input' and h['to_event']==panel['event_id']),None)
        outer=element('g',clip_path=f'url(#chain-clip-{index})',data_frame_view=view_id,
                      data_prior_event=handoff['from_event'] if handoff else '')
        body=element('g',transform=f'translate({tx} {ty}) scale({scale}) translate({-x0} {-y0})')
        for node in original:
            if node.get('id') in {'molecular-structures','electron-arrows','electron-markers','atom-ids'}:
                if side=='output' and node.get('id')=='electron-arrows':continue
                clone=copy.deepcopy(node)
                if node.get('id')!='electron-arrows':
                    # Select the actual side, rather than relying on a clipPath
                    # to hide another chemical state (QtSvg clipping is limited).
                    for child in list(clone):
                        if child.get('data-side')!=side:clone.remove(child)
                for item in clone.iter():
                    if item.get('id'):item.set('id',f'frame.{index}.'+item.get('id'))
                body.append(clone)
            elif node.get('id')=='reaction-flow':
                for child in node:
                    if child.get('data-separator-side')==side:body.append(copy.deepcopy(child))
        outer.append(body);root.append(outer)
    from .render import arrow_line
    for i,(x,y) in enumerate(positions[:-1]):
        nx,ny=positions[i+1]
        if ny==y:
            start=(x+tw+12,y+th/2) if nx>x else (x-12,y+th/2)
            end=(nx-12,ny+th/2) if nx>x else (nx+tw+12,ny+th/2)
        else:start,end=(x+tw/2,y+th+53),(nx+tw/2,ny-7)
        kind=transitions[i]
        if kind=='gap':
            label(root,'路径未衔接',(start[0]+end[0])/2,(start[1]+end[1])/2,size=12,text_anchor='middle',fill='#985936')
        else:
            arrow_line(root,start,end,resonance=kind=='resonance')
            event_id=frames[i][0]['event_id']
            title=arrow_titles.get(event_id) if kind=='reaction' and frames[i][2]=='input' else None
            if title:
                horizontal=ny==y
                label(root,str(title),(start[0]+end[0])/2+(0 if horizontal else 12),
                      (start[1]+end[1])/2-(12 if horizontal else 3),size=14,
                      text_anchor='middle' if horizontal else 'start',fill='#006b65',
                      **{'data-event-label':event_id})
    from .figure_keys import append_keys
    keys=append_keys(root,folds)
    from .local_rules import append_rules
    shown={p['event_id'] for p in events}
    append_rules(root,[meta['source_view_id'] if side=='input' else meta['target_view_id']
                       for _,meta,side,_ in frames], [e for e in document.events if e.event_id in shown])
    path=directory/'mechanism_chain.svg';write_svg(root,path);svg_to_png(path,path.with_suffix('.png'))
    return {'svg':path.name,'png':path.with_suffix('.png').name,'frames':len(frames),'events':len(events),'folded_columns':columns,
        'actual_views_only':True,'handoffs':handoffs,'abbreviation_key':keys,'frame_views':[
            m['source_view_id'] if side=='input' else m['target_view_id'] for _,m,side,_ in frames],
        'view_switches':[p['view_switch'] for p in events if p.get('view_switch')]}
