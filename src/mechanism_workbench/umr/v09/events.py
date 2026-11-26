"""Finite role-bound transactions; names of whole reactions never dispatch here."""
from ..ast import Pair,Split,Shift
from ..chemistry import FormalError,Unsupported,order,atoms_by_id
from ..v02.electrons import atom_to_bond,bond_to_atom,bond_to_bond
from .radical import SPECS as RADICALS,compile as compile_radical
from .metal_events import SPECS as METAL_SPECS,compile as compile_metal
from .eliminations import SPECS as ELIMINATIONS,compile as compile_elimination

REGISTRY={
 'EPISULFONE_SO2':{'roles':['s','a','b'],'hydrogen_transfer':'none','source':'three-membered episulfone topology; case evidence required'},
 'ANION_DECARBOXYLATION':{'roles':['r','carbon','oxo','oxygen'],'hydrogen_transfer':'none','source':'explicit carboxylate to carbanion + CO2 formal channel'},
 'RADICAL_DECARBOXYLATION':{'roles':['r','carbon','oxo','oxygen'],'hydrogen_transfer':'none','source':'explicit carboxyl radical to radical + CO2 channel'},
 'DIAZO_N2_EXTRUSION':{'roles':['carbon','proximal_n','terminal_n'],'hydrogen_transfer':'none','source':'C=N+=N- to declared singlet carbene + N2'},
 'GASSMAN_2_3':{'roles':['a','s','n','d','e'],'hydrogen_transfer':'none','source':'C-S-N-C-C onium ylide [2,3] topology, not full named mechanism'},
 'CARBENE_1_2_H_MIGRATION':{'roles':['center','donor','hydrogen'],'hydrogen_transfer':'hydrogen_migration','source':'source-declared paired carbene 1,2-H migration; kinetics and E/Z not inferred'},
 'ALKYLIDENECARBENE_1_2_H_MIGRATION':{'roles':['center','donor','hydrogen'],'hydrogen_transfer':'hydrogen_migration','rule_revision':'0.9-r7','source':'Kurti/Czako Strategic Applications, printed p402 / PDF456: paired alkylidenecarbene C=C to alkyne through explicit 1,2-H migration; not a ground-state or rate prediction'},
 'PHOSPHAZIDE_N2_EXTRUSION':{'roles':['phosphorus','organic_n','middle_n','terminal_n'],'hydrogen_transfer':'none','source':'source-declared cyclic P-N-N-N fragment; retained N2'},
}
REGISTRY.update({k:{'roles':v[0].split(),'source':v[1]} for k,v in RADICALS.items()})
REGISTRY.update({k:{'roles':v,'source':'finite ionic partition; original-case evidence required'} for k,v in METAL_SPECS.items()})
REGISTRY.update({k:{'roles':v,'source':'source-conditioned finite cyclic topology; geometry unverified'} for k,v in ELIMINATIONS.items()})
from .polar_metal import SPECS as POLAR_METAL_SPECS
REGISTRY.update({k:{'roles':list(v.elements),'source':v.scope,'execution_adapter':'polar_metal/0.9-r3'} for k,v in POLAR_METAL_SPECS.items()})
from .spin_events import SPECS as SPIN_SPECS
REGISTRY.update(SPIN_SPECS)
from .cluster_events import SPECS as CLUSTER_SPECS
REGISTRY.update(CLUSTER_SPECS)
from .pcet import SPECS as PCET_SPECS
REGISTRY.update(PCET_SPECS)
from .silyl import SPECS as SILYL_SPECS
REGISTRY.update(SILYL_SPECS)
from .textbook_events import SPECS as TEXTBOOK_SPECS, compile as compile_textbook
REGISTRY.update(TEXTBOOK_SPECS)
from .catalytic_profiles import SPECS as CATALYTIC_SPECS
REGISTRY.update(CATALYTIC_SPECS)
from .cross_coupling import SPECS as COUPLING_SPECS
REGISTRY.update(COUPLING_SPECS)
from .metathesis import SPECS as OM_SPECS
REGISTRY.update({name:{'roles':list(roles),'event_type':'POLAR',
    'execution_adapter':'umr-olefin-metathesis/1.0',
    'source':'Explicit catalyst convention and source-bound Chauvin coordination/electron transaction; kinetics unknown'}
    for name,roles in OM_SPECS.items()})

def event_type(name):
    if name in COUPLING_SPECS:return COUPLING_SPECS[name]['event_type']
    if name in TEXTBOOK_SPECS:return TEXTBOOK_SPECS[name]['event_type']
    if name in PCET_SPECS:return 'SET'
    if name in CLUSTER_SPECS:return 'COORDINATION'
    if name in SPIN_SPECS:return 'ISC'
    if name in {'DISSOCIATIVE_NA_REDUCTION','NA_RADICAL_REDUCTION','LI_RADICAL_REDUCTION','LI_CONJUGATED_REDUCTION'} or name.startswith('METAL_') and name in RADICALS:return 'SET'
    if name in RADICALS or name=='RADICAL_DECARBOXYLATION':return 'RADICAL'
    if name in ELIMINATIONS or name in {'EPISULFONE_SO2','PHOSPHAZIDE_N2_EXTRUSION','GASSMAN_2_3'}:return 'PERICYCLIC'
    return 'POLAR'
def compile(source,event):
    if event.general_event not in REGISTRY:raise Unsupported('General event profile not registered')
    spec=REGISTRY[event.general_event];roles=event.general_roles
    if set(roles)!=set(spec['roles']) or len(set(roles.values()))!=len(roles):raise FormalError('Exact distinct general event roles required')
    ids=[roles[x] for x in spec['roles']];ats=atoms_by_id(source)
    if not set(ids)<=ats.keys():raise FormalError('Missing role atom')
    if not event.evidence or not event.mechanism_assumption:raise FormalError('Event requires evidence and explicit mechanism assumption')
    expected=event_type(event.general_event)
    if event.event_type!=expected:raise FormalError('General profile/event classification mismatch; expected '+expected)
    if event.general_event in TEXTBOOK_SPECS:return compile_textbook(source,event)
    if event.general_event in RADICALS:return compile_radical(source,event)
    if event.general_event in METAL_SPECS:return compile_metal(source,event),[]
    if event.general_event in ELIMINATIONS:return compile_elimination(source,event)
    def require(a,b,n):
        if order(source,a,b)!=n:raise FormalError(f'Required bond {a}-{b} order {n}')
    t=event.general_event
    if t in {'CARBENE_1_2_H_MIGRATION','ALKYLIDENECARBENE_1_2_H_MIGRATION'}:
        c,d,h=ids
        if [ats[i].element for i in ids]!=['C','C','H']:raise FormalError('C-C-H migration roles required')
        center=next((x for x in source.reactive_centers if x.atom_id==c),None)
        if center is None or center.profile!='carbene_paired' or center.environment!='free':raise Unsupported('Explicit free paired-carbene occupation required')
        if source.electronic_state.coupling=='high_spin':raise Unsupported('requires_spin_transition')
        # Separate selector: old alkyl-carbene semantics remain unchanged.
        starting_order=2 if t=='ALKYLIDENECARBENE_1_2_H_MIGRATION' else 1
        for a,b,n in [(c,d,starting_order),(d,h,1),(c,h,0)]:require(a,b,n)
        return [Shift(a=h,b=d,c=c),Pair(a=c,b=h)],[]
    if t=='PHOSPHAZIDE_N2_EXTRUSION':
        p,n,m,t=ids
        if [ats[i].element for i in ids]!=['P','N','N','N']:raise FormalError('P-N-N-N roles required')
        for a,b,o in [(p,n,1),(n,m,1),(m,t,2),(t,p,1)]:require(a,b,o)
        return [Shift(a=n,b=m,c=t),Shift(a=t,b=p,c=n)],[m,t]
    if t=='EPISULFONE_SO2':
        s,a,b=ids
        if [ats[i].element for i in ids]!=['S','C','C']:raise Unsupported('Episulfone requires S-C-C ring')
        for x,y in [(s,a),(a,b),(b,s)]:require(x,y,1)
        oxy=[i for i in ats if ats[i].element=='O' and order(source,s,i)==2]
        if len(oxy)!=2:raise FormalError('SO2 release requires two preserved oxo atoms')
        return [Shift(a=s,b=a,c=b),Split(a=b,b=s,target=s)],[s,*oxy]
    if t in {'ANION_DECARBOXYLATION','RADICAL_DECARBOXYLATION'}:
        r,c,o,x=ids
        if [ats[i].element for i in [c,o,x]]!=['C','O','O']:raise Unsupported('Explicit carboxyl group required')
        require(r,c,1);require(c,o,2);require(c,x,1)
        if t.startswith('ANION'):
            if ats[x].formal_charge!=-1 or ats[x].radical_electrons!=0:raise FormalError('Carboxylate donor required')
            actions=[Split(a=c,b=r,target=r),Pair(a=x,b=c)]
        else:
            if ats[x].radical_electrons!=1 or ats[x].formal_charge:raise FormalError('Carboxyl radical required')
            actions=[atom_to_bond(x,x,c,'pi1'),bond_to_bond(r,c,'sigma',0,x,c,'pi1'),bond_to_atom(r,c,'sigma',1,r)]
        return actions,[c,o,x]
    if t=='DIAZO_N2_EXTRUSION':
        c,n,t=ids
        if [ats[i].element for i in ids]!=['C','N','N'] or (ats[n].formal_charge,ats[t].formal_charge)!=(1,-1):raise FormalError('Explicit diazo charge topology required')
        require(c,n,2);require(n,t,2)
        return [Split(a=n,b=c,target=c),Split(a=c,b=n,target=n),Pair(a=t,b=n)],[n,t]
    a,s,n,d,e=ids
    if ''.join(ats[i].element for i in ids)!='CSNCC' or (ats[a].formal_charge,ats[s].formal_charge)!=(-1,1):raise Unsupported('Exact C-S-N-C-C ylide profile required')
    for x,y,o in [(a,s,1),(s,n,1),(n,d,1),(d,e,2),(a,e,0)]:require(x,y,o)
    return [Pair(a=a,b=e),Shift(a=e,b=d,c=n),Split(a=n,b=s,target=s)],[]
