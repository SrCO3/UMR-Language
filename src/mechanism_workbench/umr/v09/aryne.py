"""Finite aryne profiles; exact transactions, no inferred stable intermediate."""
import json
from collections import Counter
from ..ast import Pair,Split,Shift
from ..chemistry import FormalError,Unsupported,atoms_by_id,bond_valences,order
from ..v04.ast import AryneRegion
from ..v04.ledger import check_isolated_ring
from .ledger import validate_view

RULE_VERSION='umr-aryne/0.9-r2'
SOURCE='https://authors.library.caltech.edu/records/k7yft-rxk49'
SPECS={
 'ARYNE_CONCERTED_ELIMINATION':{'event_type':'POLAR','roles':['ring','alpha_carbon','beta_carbon','hydrogen','leaving_group','base'],
   'scope':'explicit grouped elimination hypothesis; isolated substituted six-carbon ring, ortho H, Cl/Br, NH2-; no kinetic claim'},
 'ARYNE_NUCLEOPHILIC_ADDITION':{'event_type':'POLAR','roles':['aryne_a','aryne_b','nucleophile'],
   'scope':'declared ortho formal-triple six-carbon aryne; closed-shell NH2-, O-/S- with one bond; two regiochannels separate'},
 'ARYL_CARBA_PROTONATION':{'event_type':'POLAR','roles':['carbon','donor','hydrogen'],
   'scope':'resolved carbon-anion occupation + explicit NH3 H; preserve unknown spectator-metal spin'},
}

def execute(source,event):
    from .execution import transact
    validate_view(source)
    if event.event_type!='POLAR' or not event.evidence:raise FormalError('Aryne profile requires POLAR type and explicit evidence/assumption')
    if source.reactive_centers:raise Unsupported('Simultaneous special-center/aryne coupling is not registered')
    if source.electronic_state.character in {'open_shell','open_shell_singlet'}:raise Unsupported('Open-shell aryne path needs a separate electronic profile')
    spec=SPECS[event.reactive_template];b=event.bindings;atoms=atoms_by_id(source);val=bond_valences(source)
    if set(b)!=set(spec['roles']):raise FormalError('Exactly the registered aryne atom roles are required')
    if event.template_id or event.metal_template or event.metal_participation:
        raise Unsupported('Aryne polar profiles do not authorize an additional metal or unrelated template')
    scalar=[]
    for key in spec['roles']:
        x=b.get(key)
        if key=='ring':
            if not isinstance(x,list) or len(x)!=6 or len(set(x))!=6 or not set(x)<=set(atoms):raise FormalError('Six distinct ring atoms required')
        elif type(x) is not int or x not in atoms:raise FormalError('Missing aryne atom role: '+key)
        else:scalar.append(x)
    if len(scalar)!=len(set(scalar)):raise FormalError('Aryne scalar roles must be distinct')
    def require(ok,message):
        if not ok:raise FormalError(message)
    created=None;consumed=None
    if event.reactive_template=='ARYNE_CONCERTED_ELIMINATION':
        if event.mechanism_assumption!='concerted_elimination':raise FormalError('Explicit concerted_elimination hypothesis required; no mechanism inferred from net products')
        ring=b['ring'];a,z,h,x,q=(b[k] for k in ['alpha_carbon','beta_carbon','hydrogen','leaving_group','base'])
        check_isolated_ring(source,ring)
        require(a in ring and z in ring and h not in ring and x not in ring and q not in ring,'Aryne elimination role regions overlap')
        require(all(atoms[i].element=='C' and atoms[i].formal_charge==0 and atoms[i].radical_electrons==0 for i in ring),'Neutral closed-shell six-carbon input required')
        require([order(source,i,ring[(j+1)%6]) for j,i in enumerate(ring)] in [[1,2,1,2,1,2],[2,1,2,1,2,1]],'Alternating ring required')
        require(order(source,a,z)==2 and atoms[x].element in {'Cl','Br'} and order(source,a,x)==1,'Supported ortho aryl-halide sites required')
        require(atoms[h].element=='H' and order(source,z,h)==1,'Designated transferable H must be explicit')
        require((atoms[q].element,atoms[q].formal_charge,atoms[q].nonbonding_electrons,atoms[q].radical_electrons,val[q])==('N',-1,4,0,2),'Base requires NH2- occupation')
        require(atoms[q].hydrogens+sum(atoms[i].element=='H' for i in atoms if order(source,q,i))==2,'Base scope is NH2-')
        actions=[Pair(a=q,b=h),Shift(a=h,b=z,c=a),Split(a=a,b=x,target=x)]
        created=AryneRegion(region_id=event.event_id+'.aryne',ring_atom_ids=ring,aryne_atom_ids=[a,z],
            evidence=event.evidence,existence_basis='formal_model',reaction_context=spec['scope'])
    elif event.reactive_template=='ARYNE_NUCLEOPHILIC_ADDITION':
        a,z,q=(b[k] for k in ['aryne_a','aryne_b','nucleophile'])
        regions=[r for r in source.aryne_regions if set(r.aryne_atom_ids)=={a,z} and r.representation=='ortho_formal_triple']
        require(len(regions)==1 and order(source,a,z)==3,'One declared validated 1,2-aryne required')
        require(not order(source,q,a) and not order(source,q,z),'Nucleophile must not already bond to aryne sites')
        nu=atoms[q];nh=nu.hydrogens+sum(atoms[i].element=='H' for i in atoms if order(source,q,i))
        amide=(nu.element,nu.formal_charge,nu.nonbonding_electrons,nu.radical_electrons,val[q],nh)==('N',-1,4,0,2,2)
        oxothio=nu.element in {'O','S'} and (nu.formal_charge,nu.nonbonding_electrons,nu.radical_electrons,val[q])==(-1,6,0,1)
        if not (amide or oxothio):raise Unsupported('Aryne nucleophile profile currently limited to NH2-/monovalent O-/S-')
        actions=[Pair(a=q,b=a),Split(a=a,b=z,target=z)];consumed=regions[0].region_id
    else:
        c,n,h=(b[k] for k in ['carbon','donor','hydrogen'])
        ca,nu=atoms[c],atoms[n]
        require((ca.element,ca.formal_charge,ca.nonbonding_electrons,ca.radical_electrons,val[c])==('C',-1,2,0,3),
                'Resolved three-valent carbon anion required')
        require((nu.element,nu.formal_charge,nu.nonbonding_electrons,nu.radical_electrons,val[n])==('N',0,2,0,3),
                'Resolved neutral NH3 donor required')
        require(atoms[h].element=='H' and order(source,n,h)==1 and not order(source,c,h),'Designated H must be explicit on NH3')
        require(nu.hydrogens+sum(atoms[j].element=='H' for j in atoms if order(source,n,j))==3,'Only NH3 donor registered')
        actions=[Pair(a=c,b=h),Split(a=h,b=n,target=n)]
    key=lambda x:json.dumps(x.model_dump(mode='json'),sort_keys=True)
    if event.actions is not None and Counter(map(key,event.actions))!=Counter(map(key,actions)):
        raise FormalError('Explicit aryne arrows disagree with bound profile')
    def after(v):
        if consumed:v.aryne_regions=[r for r in v.aryne_regions if r.region_id!=consumed]
        if created:v.aryne_regions.append(created)
    target=transact(source,actions,event.product_state_id,event.product_view_id,after=after)
    event.actions=actions;event.reactive_rule_version=RULE_VERSION
    return target,{'event_id':event.event_id,'event_type':'POLAR','reactive_template':event.reactive_template,
        'rule_version':RULE_VERSION,'atomic_commit':True,'execution_mode':'electron_actions','formal_status':'passed',
        'compiled_actions':[a.model_dump(mode='json') for a in actions],'scope':spec['scope'],
        'source':SOURCE,'chemical_feasibility':'unknown','regioselectivity':'not_predicted','human_approved':False}
