"""Sparse display of ledger electrons; marker centers are separate from curves."""
import math
from .palette import PURPLE, TEAL_TEXT
from .geometry import add,distance,mul,normal,sub,unit,segment_distance
from .svg import element


def place_electrons(view,requests,coords,parent,obstacles,config,side='input'):
    atoms={a.atom_id:a for a in view.atoms};by_atom={}
    for req in requests:
        if req['source']['type'] in {'lone_pair','radical','metal_owned_pair'}:
            by_atom.setdefault(req['source']['atom_ids'][0],[]).append(req)
    reserved=list(obstacles);objects=[];annotations=[]
    for i,a in atoms.items():
        if a.radical_electrons is None:
            if by_atom.get(i):raise ValueError('Cannot draw electron-slot arrows from an unknown metal occupation')
            annotations.append({'atom_id':i,'kind':'unknown_occupation','ledger_electrons':a.nonbonding_electrons,'ledger_radicals':None})
            continue
        selected=by_atom.get(i,[]);slots={}
        for req in selected:slots[(req['source']['type'],req['source'].get('slot',0))]=req
        for slot in range(a.radical_electrons):slots.setdefault(('radical',slot),None)
        # Spectator pairs may be omitted. Participating pairs are never omitted.
        pair_slots=sum(kind!='radical' for kind,slot in slots)
        available=max(0,((a.nonbonding_electrons or 0)-a.radical_electrons)//2)
        special=any(c.atom_id==i for c in getattr(view,"reactive_centers",[]))
        if (selected or special) and available>pair_slots and pair_slots<config.lone_pair_display_limit:
            for slot in range(min(available,config.lone_pair_display_limit)):slots.setdefault(('lone_pair',slot),None)
        for (kind,slot),req in slots.items():
            candidates=[];center=coords[i]
            for radius in (19,25,32,40):
                for k in range(16):
                    angle=-math.pi/2+k*math.pi/8+slot*.4;pos=add(center,(radius*math.cos(angle),radius*math.sin(angle)))
                    box=[pos[0]-6,pos[1]-6,pos[0]+6,pos[1]+6]
                    overlaps=0
                    for o in reserved:
                        b=o['box']
                        if box[2]<b[0] or box[0]>b[2] or box[3]<b[1] or box[1]>b[3]:continue
                        if o['kind']=='bond' and not o.get('filled') and all(segment_distance(pos,x,y)>7 for x,y in o.get('segments',[])):continue
                        overlaps+=1
                    # Prefer 90 degrees above the atom, then nearby free sectors.
                    deviation=min(k,16-k)*.45
                    approach=0
                    if req:
                        # Give a reacting pair an open exit towards its actual
                        # acceptor. The upper default is still used for idle pairs.
                        toward=unit(sub(req['target']['reference_xy'],center))
                        approach=12*(1-math.cos(angle)*toward[0]-math.sin(angle)*toward[1])
                    candidates.append((overlaps*10000+approach+deviation+(radius-19)*.3,pos,box))
            _,pos,box=min(candidates,key=lambda t:t[0]);tangent=normal(unit(sub(pos,center)))
            points=[pos] if kind=='radical' else [add(pos,mul(tangent,-2.8)),add(pos,mul(tangent,2.8))]
            for p in points:parent.append(element('circle',cx=p[0],cy=p[1],r=1.9,fill=PURPLE if kind=='radical' else TEAL_TEXT,data_atom_id=i,data_electron_kind=kind,data_side=side))
            obstacle={'id':f'{side}.electron.{i}.{kind}.{slot}','kind':'electron','atom_ids':[i],'box':box}
            objects.append(obstacle);reserved.append(obstacle)
            annotations.append({'atom_id':i,'kind':kind,'slot':slot,'position':pos,'ledger_electrons':a.nonbonding_electrons,'ledger_radicals':a.radical_electrons})
            for action in selected:
                if (action['source']['type'],action['source'].get('slot',0))==(kind,slot):
                    action['source'].update(positions=[pos],marker_position=list(pos),atom_reference_xy=list(center))
    return objects,annotations
