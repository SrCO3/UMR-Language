"""Display policy only. Never add/remove H atoms or alter the chemical ledger."""
RULE_VERSION="umr-hydrogen-depiction/0.3"


def carbonyl_sites(view):
    atoms={a.atom_id:a for a in view.atoms};sites=set()
    for b in view.bonds:
        if b.order==2 and {atoms[b.a].element,atoms[b.b].element}=={'C','O'}:
            c=b.a if atoms[b.a].element=='C' else b.b;sites.add(c)
    return sites


def document_hydrogens(document,event_reports):
    views=document.views;atoms={a.atom_id:a for v in views for a in v.atoms}
    protected=set()
    compiled={r['event_id']:r.get('compiled_actions') for r in event_reports}
    for event in document.events:
        ids=set()
        for value in event.bindings.values():ids.update(value if isinstance(value,list) else [value])
        actions=compiled.get(event.event_id) or event.actions or []
        for action in actions:
            a=action.model_dump(mode='json') if hasattr(action,'model_dump') else action
            ids.update(v for k,v in a.items() if k in {'a','b','c','target'})
            for end in ('source','destination'):ids.update(v for k,v in a.get(end,{}).items() if k in {'a','b','atom_id'})
        protected|={i for i in ids if i in atoms and atoms[i].element=='H'}
    return protected


def display_policy(view,required=()):
    atoms={a.atom_id:a for a in view.atoms};sites=carbonyl_sites(view)
    required=set(required);stereo_references=set()
    for obj in [*view.atoms,*view.bonds]:
        if obj.stereo.get('status')=='specified':
            stereo_references.update(obj.stereo.get('neighbor_order',[]));stereo_references.update(obj.stereo.get('reference_atoms',[]))
            if hasattr(obj,'atom_id'):stereo_references.add(obj.atom_id)
    hidden=[];visible=[];attached=[]
    for atom in view.atoms:
        if atom.element=='C' and atom.atom_id in sites and atom.hydrogens:
            attached.append({'parent_atom_id':atom.atom_id,'count':atom.hydrogens,'style':'separate_H_bond','reason':'directly_bonded_to_carbonyl_carbon','numbered_atom_ids':[]})
        if atom.element!='H':continue
        bonds=[b for b in view.bonds if atom.atom_id in {b.a,b.b}]
        parent=atoms[bonds[0].b if bonds[0].a==atom.atom_id else bonds[0].a] if len(bonds)==1 else None
        reason='ordinary_carbon_hydrogen'
        if atom.atom_id in required:reason='reaction_participant'
        elif atom.isotope:reason='isotope'
        elif atom.formal_charge or atom.radical_electrons or parent is None:reason='independent_hydrogen_species'
        elif atom.atom_id in stereo_references or parent.atom_id in stereo_references:reason='stereochemical_reference'
        elif parent.atom_id in sites:reason='directly_bonded_to_carbonyl_carbon'
        elif parent.element!='C':reason='heteroatom_or_metal_hydrogen'
        if reason=='ordinary_carbon_hydrogen':hidden.append(atom.atom_id)
        else:visible.append({'atom_id':atom.atom_id,'reason':reason})
    return {'rule_version':RULE_VERSION,'hidden_atom_ids':hidden,'explicit_visible':visible,'attached_hydrogen_labels':attached,
        'carbonyl_site_atom_ids':sorted(set(atoms)&sites),'all_hydrogen_atoms_retained_in_AST':True}


def draw_attached_hydrogens(view,coords,root,obstacles,config,policy):
    """Expand unnumbered attached H as depiction glyphs, never invented AST atoms."""
    import math
    from .geometry import add,mul,sub,segment_distance
    from .svg import element,label,text_box

    atoms={a.atom_id:a for a in view.atoms};reserved=list(obstacles);extra=[];glyphs=[]
    hidden=set(policy['hidden_atom_ids'])
    for entry in policy['attached_hydrogen_labels']:
        parent=entry['parent_atom_id'];center=coords[parent];count=entry['count']
        neighbors=[coords[b.b if b.a==parent else b.a] for b in view.bonds
            if parent in {b.a,b.b} and (b.b if b.a==parent else b.a) not in hidden]
        angles=sorted(math.atan2(p[1]-center[1],p[0]-center[0]) for p in neighbors)
        gaps=[(b-a,a) for a,b in zip(angles,angles[1:]+[angles[0]+2*math.pi])]
        gap,base=max(gaps);ideal=[base+gap*(j+1)/(count+1) for j in range(count)]
        candidates=[]
        for turn in (0,math.pi/18,-math.pi/18,math.pi/9,-math.pi/9):
            proposed=[];hits=0
            for angle in ideal:
                direction=(math.cos(angle+turn),math.sin(angle+turn))
                point=add(center,mul(direction,config.bond_pixels))
                start=add(center,mul(direction,13 if atoms[parent].formal_charge or atoms[parent].isotope else 0))
                end=add(point,mul(direction,-10));box=text_box('H',point[0]-config.font_pixels*.36,point[1]+config.font_pixels*.35,config.font_pixels)
                for o in reserved:
                    b=o['box']
                    if not (box[2]+2<b[0] or box[0]-2>b[2] or box[3]+2<b[1] or box[1]-2>b[3]):hits+=1
                    for k in range(3,11):
                        p=add(start,mul(sub(end,start),k/10))
                        if o['kind']=='bond':
                            hits+=any(segment_distance(p,a,b)<3 for a,b in o.get('segments',[]))
                        elif b[0]-2<=p[0]<=b[2]+2 and b[1]-2<=p[1]<=b[3]+2:hits+=1
                proposed.append((point,start,end,box))
            candidates.append((hits*10000+abs(turn),proposed,hits))
        _,chosen,hits=min(candidates,key=lambda c:c[0])
        for slot,(point,start,end,box) in enumerate(chosen):
            display_id=f'attached-h.{parent}.{slot+1}'
            attrs={'data-parent-atom-id':parent,'data-display-id':display_id,'data-display-only':'true'}
            root.append(element('path',d=f'M {start[0]} {start[1]} L {end[0]} {end[1]}',
                stroke='#000000',stroke_width=2,fill='none',**{'class':'attached-hydrogen-bond'},**attrs))
            label(root,'H',point[0]-config.font_pixels*.36,point[1]+config.font_pixels*.35,size=config.font_pixels,
                fill='#000000',**{'class':'attached-hydrogen'},**attrs)
            glyphs.append({'display_id':display_id,'parent_atom_id':parent,'slot':slot+1,
                'local_xy':list(point),'bond_local_xy':[list(start),list(end)],'display_only':True,
                'atom_id':None,'source':'parent.hydrogens','layout_overlap_count':hits})
            objects=[{'id':display_id+'.label','kind':'text','box':box,'atom_ids':[parent]},
                {'id':display_id+'.bond','kind':'bond','box':[min(start[0],end[0])-2,min(start[1],end[1])-2,max(start[0],end[0])+2,max(start[1],end[1])+2],
                 'atom_ids':[parent],'segments':[(start,end)],'width':2}]
            extra.extend(objects);reserved.extend(objects)
    policy['display_hydrogens']=glyphs
    return extra
