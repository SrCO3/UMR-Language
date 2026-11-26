"""One finite aqueous nitroxyl/hypochlorite proton-electron ledger.

This is an explicitly grouped redox/proton bookkeeping channel, not a claim
that two nitroxyls, water and hypochlorite collide in an elementary event.
"""
import json
from collections import Counter
from ..chemistry import atoms_by_id,order,bond_valences,FormalError,Unsupported
from ..v07.ast import Fishhook
from .ledger import validate_view
from .ast import RULE_VERSION
SPECS={'NITROXYL_HYPOCHLORITE_PCET':{'roles':['n1','o1','n2','o2','oxidant_o','halogen','water_o','hydrogen'],
       'source':'aqueous two-nitroxyl oxidation; all H/O/Cl retained; grouped formal channel'}}
def execute(v,e):
    from .execution import transact
    validate_view(v);r=e.general_roles;ats=atoms_by_id(v);vals=bond_valences(v)
    if set(r)!=set(SPECS[e.general_event]['roles']) or len(set(r.values()))!=8 or not set(r.values())<=ats.keys():raise FormalError('Eight distinct explicit PCET roles required')
    if e.event_type!='SET' or e.electron_transfers or e.electronic_outcome:raise FormalError('Registered grouped SET/PCET classification required')
    if not e.evidence or not e.mechanism_assumption or e.conditions.get('medium')!='aqueous':raise FormalError('Aqueous condition, evidence and grouped-channel assumption required')
    if v.electronic_state.character!='unknown' or v.electronic_state.coupling!='unknown':raise Unsupported('Defined spin coupling needs a separate PCET model')
    n1,o1,n2,o2,o,x,w,h=[r[k] for k in SPECS[e.general_event]['roles']]
    if [ats[i].element for i in [n1,o1,n2,o2,o,x,w,h]]!=['N','O','N','O','O','Cl','O','H']:raise FormalError('PCET role elements disagree')
    for n,ox in [(n1,o1),(n2,o2)]:
        if order(v,n,ox)!=1 or (vals[n],ats[n].formal_charge,ats[n].nonbonding_electrons,ats[n].radical_electrons)!=(3,0,2,0) or (vals[ox],ats[ox].formal_charge,ats[ox].nonbonding_electrons,ats[ox].radical_electrons)!=(1,0,5,1):raise Unsupported('Explicit neutral dialkyl nitroxyl profile required')
        if len([i for i in ats if order(v,n,i)==1 and ats[i].element=='C'])!=2:raise Unsupported('Nitroxyl N must have two carbon substituents')
    if order(v,o,x)!=1 or (vals[o],ats[o].formal_charge,ats[o].nonbonding_electrons,ats[o].radical_electrons)!=(1,-1,6,0) or vals[x]!=1:raise FormalError('Explicit hypochlorite anion required')
    hs=[i for i in ats if order(v,w,i)]
    if len(hs)!=2 or h not in hs or any(ats[i].element!='H' for i in hs) or ats[w].hydrogens or ats[w].formal_charge or ats[w].radical_electrons or ats[w].nonbonding_electrons!=4:raise FormalError('Water must have both hydrogens explicit')
    if order(v,o,h):raise FormalError('Transferred hydrogen already bound to oxidant')
    def fish(src,dst):return Fishhook.model_validate({'op':'FISHHOOK','source':src,'destination':dst})
    def lp(i,slot):return dict(kind='lone_pair',atom_id=i,pair_index=0,slot=slot)
    def pairdest(i):return dict(kind='lone_pair',atom_id=i)
    actions=[fish(dict(kind='radical',atom_id=i,slot=0),pairdest(o)) for i in [o1,o2]]
    for n,ox in [(n1,o1),(n2,o2)]:
        actions += [fish(lp(n,k),dict(kind='bond',a=n,b=ox,pair='pi1')) for k in [0,1]]
    actions += [fish(dict(kind='bond',a=o,b=x,pair='sigma',slot=k),pairdest(x)) for k in [0,1]]
    actions += [fish(lp(o,k),dict(kind='bond',a=o,b=h,pair='sigma')) for k in [0,1]]
    actions += [fish(dict(kind='bond',a=w,b=h,pair='sigma',slot=k),pairdest(w)) for k in [0,1]]
    key=lambda x:json.dumps(x.model_dump(mode='json'),sort_keys=True)
    if e.actions is not None and Counter(map(key,e.actions))!=Counter(map(key,actions)):raise FormalError('PCET actions differ from registered occupied-slot transaction')
    target=transact(v,actions,e.product_state_id,e.product_view_id);e.actions=actions
    return target,{'event_id':e.event_id,'event_type':'SET','general_event':e.general_event,'rule_version':RULE_VERSION,
        'execution_mode':'electron_actions','formal_status':'passed','atomic_commit':True,
        'compiled_actions':[a.model_dump(mode='json') for a in actions],'chemical_feasibility':'unknown',
        'mechanistic_resolution':'grouped_redox_proton_ledger_not_asserted_elementary',
        'retained_coproducts':{'chloride':[x],'hydroxides':[[o,h],[w,next(i for i in hs if i!=h)]]},'ledger':validate_view(target)}
