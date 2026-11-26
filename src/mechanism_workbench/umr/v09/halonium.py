"""Explicit, bounded halonium opening; does not authorize ordinary vinyl SN2."""
import json
from collections import Counter
from ..ast import Pair, Split
from ..chemistry import atoms_by_id, order, bond_valences, FormalError, Unsupported

PROFILE='VINYL_HALONIUM_OPENING'
VERSION='umr-halonium/0.9-r1'
SPECS={PROFILE:{'event_type':'POLAR','roles':['carbon','partner','bridge','nucleophile'],
    'scope':'Explicit bridged C=C/Br+ or Cl+ and free Br-/Cl-; source-declared finite opening, not general vinyl SN2',
    'version':VERSION}}
def execute(source,event):
    from .ledger import validate_view
    from .execution import transact
    validate_view(source)
    roles=event.general_roles
    if set(roles)!={'carbon','partner','bridge','nucleophile'} or len(set(roles.values()))!=4:
        raise FormalError('Four distinct halonium roles required')
    a=atoms_by_id(source)
    if not set(roles.values())<=a.keys():raise FormalError('Missing halonium role atom')
    c,z,x,n=[roles[k] for k in ['carbon','partner','bridge','nucleophile']];b=bond_valences(source)
    if event.event_type!='POLAR' or not event.evidence:raise FormalError('POLAR and evidence required')
    if any((event.metal_template,event.rule_template,event.reactive_template,event.element_template,event.template_id,
            event.electron_transfers,event.electronic_outcome,event.metal_participation,event.released_fragments)):
        raise Unsupported('Mixed halonium/template declaration')
    if not (a[c].element==a[z].element=='C' and order(source,c,z)==2 and order(source,c,x)==order(source,z,x)==1):
        raise FormalError('Explicit three-membered vinyl halonium connectivity required')
    if not (a[x].element in {'Br','Cl'} and (b[x],a[x].formal_charge,a[x].nonbonding_electrons,a[x].radical_electrons)==(2,1,4,0)):
        raise Unsupported('Finite bridged chloronium/bromonium occupation required')
    if not (a[n].element in {'Br','Cl'} and (b[n],a[n].formal_charge,a[n].nonbonding_electrons,a[n].radical_electrons)==(0,-1,8,0)):
        raise Unsupported('Finite free halide nucleophile required')
    if any(a[i].formal_charge or a[i].radical_electrons or b[i]!=4 for i in (c,z)):
        raise FormalError('Neutral tetravalent carbon sites required')
    actions=[Pair(a=n,b=c),Split(a=c,b=x,target=x)]
    key=lambda x:json.dumps(x.model_dump(mode='json'),sort_keys=True)
    if event.actions is not None and Counter(map(key,event.actions))!=Counter(map(key,actions)):
        raise FormalError('Actions disagree with halonium opening roles')
    target=transact(source,actions,event.product_state_id,event.product_view_id)
    event.actions=actions
    return target,{'event_id':event.event_id,'event_type':'POLAR','general_event':PROFILE,
        'profile_version':VERSION,'atomic_commit':True,'formal_status':'passed','execution_mode':'electron_actions',
        'compiled_actions':[x.model_dump(mode='json') for x in actions],
        'chemical_feasibility':'unknown','stereoselectivity':'requires_external_evidence','human_approved':False}
