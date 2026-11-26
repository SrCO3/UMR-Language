"""Finite fluoride/silyl fragmentation and aryne capture, using native transactions."""
import json
from collections import Counter
from ..ast import Pair,Split
from ..chemistry import FormalError,Unsupported,atoms_by_id,bond_valences,order
from ..v04.ast import AryneRegion
from ..v04.ledger import check_isolated_ring
from ..v07.substitution import aromatic_ring
from .ledger import validate_view

VERSION='umr-silyl/0.9-r5'
ARYNE_SOURCE='https://doi.org/10.1039/D5SC03943F'
PG_SOURCE='https://doi.org/10.1039/C7RA07876E'
SPECS={
 'FLUORIDE_SILICON_ASSOCIATION':{'roles':['silicon','fluoride'],'source':ARYNE_SOURCE,'scope':'Neutral four-coordinate C/O silane plus free F-; explicit five-coordinate silicate'},
 'SILICON_ATE_CLEAVAGE':{'roles':['silicon','fluoride','leaving'],'source':ARYNE_SOURCE,'scope':'Fluorosilicate cleavage restricted to ortho aryl triflate, beta activated oxyethyl, alkynyl or O-silyl donors'},
 'ARYL_ANION_TRIFLATE_ELIMINATION':{'roles':['alpha','beta','oxygen'],'source':ARYNE_SOURCE,'scope':'Declared isolated six-carbon ortho triflate aryl anion to formal-triple aryne; ring in bindings'},
 'ARYNE_NEUTRAL_CAPTURE':{'roles':['aryne_a','aryne_b','nucleophile'],'source':'local:authored-neutral-aryne-capture-completion','scope':'Declared isolated ortho aryne plus neutral O-H or N-H paired donor; zwitterion retained; not a validated prediction of trap selectivity'},
 'SILYLETHYL_BETA_ELIMINATION':{'roles':['alpha','beta','oxygen'],'source':PG_SOURCE,'scope':'Explicit CH2(-)-CH2-O(acyl/acetal) to ethylene and oxygen anion'},
 'HEMIACETAL_ANION_FRAGMENTATION':{'roles':['oxygen','carbon','leaving_oxygen'],'source':'local:Greene4.oh_sem.deprotect','scope':'O(-)-CH2-O-R to formaldehyde and alkoxide; source-specific completion evidence needed'},
}

def execute(v,e):
    from .execution import transact
    from .ast import RULE_VERSION
    validate_view(v)
    spec=SPECS[e.general_event];r=e.general_roles;ats=atoms_by_id(v);val=bond_valences(v)
    def need(ok,message):
        if not ok:raise FormalError(message)
    need(e.event_type=='POLAR' and e.evidence and e.mechanism_assumption,'Polar profile requires explicit evidence and mechanism assumption')
    if any((e.metal_template,e.element_template,e.rule_template,e.template_id,e.reactive_template,e.electron_transfers,e.electronic_outcome,e.metal_participation,e.released_fragments)):
        raise Unsupported('Mixed silyl profile/other execution declaration')
    need(set(r)==set(spec['roles']) and len(set(r.values()))==len(r) and set(r.values())<=ats.keys(),'Exact distinct role atoms required')
    if e.general_event!='ARYL_ANION_TRIFLATE_ELIMINATION' and e.bindings:raise FormalError('Unregistered extra role bindings')
    ids=[r[k] for k in spec['roles']]
    def ns(i):return [j for j in ats if order(v,i,j)]
    def occ(i,element,b,n,q):return (ats[i].element,val[i],ats[i].nonbonding_electrons,ats[i].radical_electrons,ats[i].formal_charge)==(element,b,n,0,q)
    def hydrogens(i):return ats[i].hydrogens+sum(ats[j].element=='H' for j in ns(i))
    def triflate(o):
        ss=[j for j in ns(o) if ats[j].element=='S' and order(v,o,j)==1]
        if len(ss)!=1:return False
        s=ss[0];ox=[j for j in ns(s) if ats[j].element=='O' and order(v,s,j)==2]
        cs=[j for j in ns(s) if ats[j].element=='C' and order(v,s,j)==1]
        return (len(ns(s))==4 and len(ox)==2 and len(cs)==1 and len(ns(cs[0]))==4 and
            sum(ats[j].element=='F' and order(v,cs[0],j)==1 for j in ns(cs[0]))==3 and ats[s].formal_charge==0)
    def activated_oxygen(o,beta):
        other=[i for i in ns(o) if i!=beta]
        if not occ(o,'O',2,4,0) or len(other)!=1 or ats[other[0]].element!='C':return False
        c=other[0]
        return (any(ats[j].element=='O' and order(v,c,j)==2 for j in ns(c)) or
            hydrogens(c)==2 and any(j!=o and ats[j].element=='O' and order(v,c,j)==1 for j in ns(c)))
    created=None;consumed=None;fragment=[]
    if e.general_event=='FLUORIDE_SILICON_ASSOCIATION':
        si,f=ids
        need(occ(si,'Si',4,0,0) and len(ns(si))==4 and all(order(v,si,j)==1 and ats[j].element in {'C','O'} for j in ns(si)), 'Four explicit C/O single ligands on neutral silicon required')
        need(occ(f,'F',0,8,-1),'Free paired fluoride required; HF is not F-')
        actions=[Pair(a=f,b=si)]
    elif e.general_event=='SILICON_ATE_CLEAVAGE':
        si,f,x=ids
        need(occ(si,'Si',5,0,-1) and len(ns(si))==5 and all(order(v,si,j)==1 for j in ns(si)), 'Five-coordinate anionic silicon required')
        if any(j!=f and ats[j].element not in {'C','O'} for j in ns(si)):
            raise Unsupported('This cleavage profile retains C/O ligands plus one designated fluoride only')
        need(occ(f,'F',1,6,0) and order(v,si,f)==order(v,si,x)==1,'Designated F/Si/leaving bonds required')
        need(ats[x].radical_electrons==0 and ats[x].formal_charge==0,'Closed neutral leaving center required')
        scope=e.conditions.get('fragmentation_scope')
        if scope=='silyl_ether':ok=occ(x,'O',2,4,0)
        elif scope=='alkynyl':ok=occ(x,'C',4,0,0) and any(ats[j].element=='C' and order(v,x,j)==3 for j in ns(x))
        elif scope=='beta_heteroethyl':
            ok=occ(x,'C',4,0,0) and hydrogens(x)==2 and any(ats[b].element=='C' and hydrogens(b)==2 and order(v,x,b)==1 and any(activated_oxygen(o,b) for o in ns(b) if ats[o].element=='O') for b in ns(x) if b!=si)
        elif scope=='aryne_precursor':
            ring=aromatic_ring(v,x)
            ok=occ(x,'C',4,0,0) and ring is not None and any(j in ring and any(triflate(o) for o in ns(j) if ats[o].element=='O') for j in ns(x))
        else:raise Unsupported('Explicit registered fragmentation_scope required')
        if not ok:raise Unsupported('Substrate does not meet the declared silicon-cleavage activation scope')
        actions=[Split(a=si,b=x,target=x)]
    elif e.general_event=='ARYL_ANION_TRIFLATE_ELIMINATION':
        a,b,o=ids;ring=e.bindings.get('ring')
        need(set(e.bindings)=={'ring'} and isinstance(ring,list) and len(ring)==len(set(ring))==6 and set(ring)<=ats.keys(),'Six distinct explicit ring atoms required')
        need(a in ring and b in ring and o not in ring and all(ats[i].element=='C' and ats[i].radical_electrons==0 for i in ring),'Isolated closed six-carbon ring scope')
        check_isolated_ring(v,ring)
        need([order(v,i,ring[(j+1)%6]) for j,i in enumerate(ring)] in [[1,2,1,2,1,2],[2,1,2,1,2,1]],'Alternating ring required')
        need(occ(a,'C',3,2,-1) and occ(b,'C',4,0,0) and order(v,a,b)==2 and order(v,b,o)==1 and triflate(o),'Ortho aryl anion/triflate topology required')
        need(not v.aryne_regions,'Multiple aryne regions not supported in this generation profile')
        actions=[Pair(a=a,b=b),Split(a=b,b=o,target=o)]
        created=AryneRegion(region_id=e.event_id+'.aryne',ring_atom_ids=ring,aryne_atom_ids=[a,b],existence_basis='literature_analogy',reaction_context=spec['scope'],evidence=e.evidence)
    elif e.general_event=='ARYNE_NEUTRAL_CAPTURE':
        a,b,n=ids
        regions=[z for z in v.aryne_regions if set(z.aryne_atom_ids)=={a,b} and z.representation=='ortho_formal_triple']
        need(len(regions)==1 and order(v,a,b)==3 and not order(v,n,a) and not order(v,n,b),'Declared ortho formal-triple aryne required')
        need((occ(n,'O',2,4,0) or occ(n,'N',3,2,0)) and hydrogens(n)>=1,'Neutral paired O-H or N-H donor required')
        actions=[Pair(a=n,b=a),Split(a=a,b=b,target=b)];consumed=regions[0].region_id
    elif e.general_event=='SILYLETHYL_BETA_ELIMINATION':
        a,b,o=ids
        need(occ(a,'C',3,2,-1) and occ(b,'C',4,0,0) and hydrogens(a)==hydrogens(b)==2,'Explicit anionic CH2-CH2 region required')
        need(order(v,a,b)==order(v,b,o)==1 and activated_oxygen(o,b),'Activated ester/carbamate/acetal oxygen leaving group required')
        actions=[Pair(a=a,b=b),Split(a=b,b=o,target=o)];fragment=[a,b]
    else:
        o,c,x=ids
        need(occ(o,'O',1,6,-1) and occ(c,'C',4,0,0) and hydrogens(c)==2 and occ(x,'O',2,4,0),'Anionic formaldehyde hemiacetal occupation required')
        need(order(v,o,c)==order(v,c,x)==1,'O(-)-CH2-O connectivity required')
        actions=[Pair(a=o,b=c),Split(a=c,b=x,target=x)];fragment=[c,o]
    key=lambda a:json.dumps(a.model_dump(mode='json'),sort_keys=True)
    if e.actions is not None and Counter(map(key,e.actions))!=Counter(map(key,actions)):raise FormalError('Actions disagree with silyl profile')
    def after(w):
        if created:w.aryne_regions.append(created)
        if consumed:w.aryne_regions=[z for z in w.aryne_regions if z.region_id!=consumed]
    target=transact(v,actions,e.product_state_id,e.product_view_id,after=after)
    if fragment and any(bool(b.a in fragment)!=bool(b.b in fragment) for b in target.bonds):raise FormalError('Released fragment remains connected')
    e.actions=actions
    return target,{'event_id':e.event_id,'general_event':e.general_event,'event_type':'POLAR','rule_version':RULE_VERSION,'profile_version':VERSION,
        'atomic_commit':True,'execution_mode':'electron_actions','formal_status':'passed','compiled_actions':[a.model_dump(mode='json') for a in actions],
        'released_fragments':[fragment] if fragment else [],'scope':spec['scope'],'source':spec['source'],
        'chemical_feasibility':'unknown','human_approved':False}
