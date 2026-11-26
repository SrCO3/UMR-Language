"""Typed UMR action -> electron ports. No endpoint-state difference inference."""
import math
from collections import Counter
from .geometry import add,mul,normal,sub,unit


def requests_for(view,actions,coords,event_id,inactive=(),obstacles=(),declared_pair_donors=()):
    atoms={a.atom_id:a for a in view.atoms};orders={tuple(sorted((b.a,b.b))):b.order for b in view.bonds}
    used=Counter();requests=[];unsupported=[]
    def atom_port(i,kind,slot=0):
        if i not in coords:raise ValueError(f'Missing endpoint atom {i}')
        p=coords[i]
        charge=next((o for o in obstacles if o.get('kind')=='charge' and o.get('atom_ids')==[i]),None)
        if kind=='atom_acceptor' and charge:p=charge['position']
        radius=20 if kind in {'lone_pair','metal_owned_pair','radical'} else 15
        positions=[]
        own_labels=[o['box'] for o in obstacles if o.get('kind') in {'text','atom_id'} and o.get('atom_ids')==[i]]
        angles=[-math.pi/2+slot*.6+k*math.pi/3 for k in range(6)]
        if kind=='atom_acceptor' and not own_labels:
            radius=12 if charge else 5
            neighbors=[coords[j] for a,b in orders for j in ([b] if a==i else [a] if b==i else []) if j in coords]
            if not charge and len(neighbors)>=3:
                radius=10
                occupied=sorted(math.atan2(q[1]-p[1],q[0]-p[0]) for q in neighbors)
                gaps=sorted(((b-a,(b+a)/2) for a,b in zip(occupied,occupied[1:]+[occupied[0]+2*math.pi])),reverse=True)
                angles=[angle+offset for _,angle in gaps for offset in (0,.18,-.18)]
        for angle in angles:
            r=radius
            while r<54:
                point=add(p,(r*math.cos(angle),r*math.sin(angle)))
                if not any(b[0]-7<=point[0]<=b[2]+7 and b[1]-7<=point[1]<=b[3]+7 for b in own_labels):break
                r+=3
            positions.append(point)
        return {'type':kind,'atom_ids':[i],'slot':slot,'positions':positions,'reference_xy':list(p)}
    def bond_port(a,b,kind,pair='sigma',slot=None):
        if a not in coords or b not in coords:raise ValueError('Missing bond endpoint atom')
        p,q=coords[a],coords[b];mid=mul(add(p,q),.5);v=normal(unit(sub(q,p)))
        if kind=='bond_electrons' and tuple(sorted((a,b))) not in orders:raise ValueError('Source bond absent in actual input view')
        shift={'sigma':0,'pi1':6,'pi2':-6}.get(pair,0)
        pos=add(mid,mul(v,shift))
        if slot is not None:pos=add(pos,mul(unit(sub(q,p)),-3 if slot==0 else 3))
        return {'type':kind,'atom_ids':[a,b],'pair':pair,'slot':slot,'positions':[pos],
            'reference_xy':list(mid),'input_bond_order':orders.get(tuple(sorted((a,b))),0)}
    for index,action in enumerate(actions):
        action=action.model_dump(mode='json') if hasattr(action,'model_dump') else action
        meta={'arrow_id':f'{event_id}.arrow.{index+1:03d}','event_id':event_id,'action_id':f'{event_id}.action.{index+1:03d}',
              'action_index':index,'view_id':view.view_id,'raw_action':action,'coordinate_system':'final SVG pixels, top-left origin'}
        try:
            op=action['op']
            if op=='PAIR':
                a,b=action['a'],action['b'];slot=used[a];used[a]+=1
                conditional=a in declared_pair_donors and atoms[a].radical_electrons is None and used[a]==1 and atoms[a].nonbonding_electrons is not None and atoms[a].nonbonding_electrons>=2
                if not conditional and (atoms[a].nonbonding_electrons is None or atoms[a].radical_electrons is None or 2*used[a]>atoms[a].nonbonding_electrons-atoms[a].radical_electrons):raise ValueError('No ledger-backed unspent electron pair')
                kind='declared_metal_pair' if conditional else 'metal_owned_pair' if atoms[a].element in {'Pd','Zn'} else 'lone_pair'
                if conditional:meta['occupation_basis']='explicit conditional paired channel; orbital position unknown, no lone-pair dots invented'
                source=atom_port(a,kind,slot);target=atom_port(b,'atom_acceptor');electrons=2
            elif op in {'SPLIT','SHIFT'}:
                a,b=action['a'],action['b'];order=orders.get(tuple(sorted((a,b))),0);pair='pi1' if order>=2 else 'sigma'
                key=tuple(sorted((a,b)));debit=used[('bond',key)];used[('bond',key)]+=1
                if debit>=order:raise ValueError('Duplicate bond electron-pair consumption in drawing input')
                pair={1:'sigma',2:'pi1',3:'pi2'}[order-debit]
                source=bond_port(a,b,'bond_electrons',pair)
                target=atom_port(action['target'],'atom_acceptor') if op=='SPLIT' else bond_port(b,action['c'],'forming_bond','pi1' if orders.get(tuple(sorted((b,action['c']))),0) else 'sigma')
                electrons=2
            elif op=='FISHHOOK':
                src,dst=action['source'],action['destination'];electrons=1
                if src['kind']=='radical':
                    if atoms[src['atom_id']].radical_electrons<=src['slot']:raise ValueError('Radical slot not occupied in input ledger')
                    source=atom_port(src['atom_id'],'radical',src['slot'])
                elif src['kind']=='lone_pair':
                    a=atoms[src['atom_id']]
                    if a.radical_electrons is None or src['pair_index'] >= (a.nonbonding_electrons-a.radical_electrons)//2:raise ValueError('Unoccupied lone-pair electron slot')
                    source=atom_port(src['atom_id'],'lone_pair',src['pair_index'])
                    source['electron_slot']=src['slot']
                else:source=bond_port(src['a'],src['b'],'bond_electrons',src['pair'],src['slot'])
                target=atom_port(dst['atom_id'],'atom_acceptor') if dst['kind'] in {'atom','lone_pair'} else bond_port(dst['a'],dst['b'],'forming_bond',dst['pair'])
            else:raise ValueError(f'Unimplemented typed arrow: {op}')
            if (set(source['atom_ids'])|set(target['atom_ids']))&set(inactive):raise ValueError('Action uses declared inactive inventory; arrow retained as unsupported, not drawn as active chemistry')
            meta.update(source=source,target=target,electrons=electrons,arrow_shape='fishhook' if electrons==1 else 'full_head')
            requests.append(meta)
        except (ValueError,KeyError) as exc:unsupported.append({**meta,'status':'render_unsupported','reason':str(exc)})
    shared={}
    for req in requests:
        if req['target']['type']=='forming_bond':shared.setdefault(tuple(sorted(req['target']['atom_ids'])),[]).append(req)
    for ids,group in shared.items():
        if len(group)==1:
            req=group[0]
            if req['raw_action']['op']=='SHIFT' and not req['target']['input_bond_order']:
                # The midpoint of a not-yet-existing intermolecular bond can
                # sit on a spectator atom. Attach this single pair arrow near
                # the actual acceptor, while retaining both semantic atom IDs.
                acceptor=req['raw_action']['c'];port=atom_port(acceptor,'atom_acceptor')
                req['target'].update(positions=port['positions'],reference_xy=port['reference_xy'],placement_policy='new_bond_acceptor',acceptor_atom_id=acceptor)
            continue
        p,q=sorted([coords[i] for i in ids],key=lambda xy:(xy[0],xy[1]));axis=unit(sub(q,p))
        group.sort(key=lambda r:(sum((r['source']['reference_xy'][k]-p[k])*axis[k] for k in (0,1)),r['action_index']))
        for n,req in enumerate(group,1):
            fraction=n/(len(group)+1);point=add(p,mul(sub(q,p),fraction))
            req['target'].update(positions=[point],shared_bond_fraction=fraction,semantic_attachment_xy=list(point))
    return requests,unsupported
