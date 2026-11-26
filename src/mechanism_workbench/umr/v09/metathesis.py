"""Finite Chauvin olefin transactions with explicit initiation and pi binding.

Historical METALLACYCLE_* X2 semantics are untouched. New OM_* profiles use
the declared neutral-L Ru or X2 Mo/W convention and the shared executor.
"""
from ..ast import Shift
from ..chemistry import FormalError,Unsupported,atoms_by_id,bond_orders,identity,electron_count
from .ast import Ligand,MetalBinding
from .coordination import checks
from .ledger import validate_view

VERSION='umr-olefin-metathesis/1.0'
SPECS={
    'OM_DONOR_DISSOCIATION':('metal','donor'),
    'OM_OLEFIN_ASSOCIATION':('metal','a','b'),
    'OM_HOVEYDA_EXCHANGE_INITIATION':('metal','donor','a','b'),
    'OM_METALLACYCLE_FORMATION':('metal','alkylidene','a','b'),
    'OM_METALLACYCLE_CLEAVAGE':('metal','first','second','third'),
    'OM_PRODUCT_RELEASE':('metal','a','b'),
}


def _roles(view,event):
    ats=atoms_by_id(view);r=event.general_roles
    if set(r)!=set(SPECS[event.general_event]) or len(set(r.values()))!=len(r) or not set(r.values())<=ats.keys():
        raise FormalError('Metathesis needs exact distinct mapped roles')
    if event.event_type!='POLAR' or not event.evidence or not event.mechanism_assumption:
        raise FormalError('Metathesis is an evidenced metal-assisted event, not an unnamed pericyclic shortcut')
    m=r['metal'];el=ats[m].element
    if el not in {'Ru','Mo','W'}:raise Unsupported('Registered olefin metathesis metals are Ru/Mo/W in their finite ligand profiles')
    annotation=view.annotations.get('metathesis_representation',{})
    if annotation.get('metal_atom_id')!=m or annotation.get('carbene_convention')!=('neutral_L' if el=='Ru' else 'X2'):
        raise Unsupported('Explicit state-local metathesis carbene convention required; no automatic X2/L conversion')
    return r


def _component(view,root):
    ats=atoms_by_id(view);ids={root};stack=[root]
    from .elements import current_elements
    while stack:
        i=stack.pop()
        for a,b in bond_orders(view):
            if i not in (a,b):continue
            j=b if i==a else a
            if j not in ids and ats[j].element not in current_elements().metals:ids.add(j);stack.append(j)
    return ids


def bind_pi(view,m,a,b,event_id):
    ats=atoms_by_id(view);edges=bond_orders(view)
    if any(ats[i].element!='C' or ats[i].formal_charge or ats[i].radical_electrons!=0 for i in (a,b)) or edges.get(tuple(sorted((a,b))))!=2:
        raise Unsupported('Pi partner must be an actual neutral closed-shell C=C; alkynes/aromatics/hetero-olefins need distinct rules')
    if any(x.metal_atom_id==m and x.mode=='eta' for x in view.metal_bindings):raise FormalError('Metathesis active site already has a pi partner')
    if any(edges.get(tuple(sorted((m,i))),0) for i in (a,b)):raise FormalError('Do not disguise existing sigma ligands as an incoming pi partner')
    ids=_component(view,a)
    if b not in ids:raise FormalError('Disconnected alkene donor atoms')
    ls=[l for l in view.ligands if set(l.atom_ids)==ids]
    if len(ls)>1:raise FormalError('Ambiguous complete incoming ligand identity')
    if ls:
        ligand=ls[0];sigma=[x for x in view.metal_bindings if x.ligand_id==ligand.ligand_id and x.mode=='sigma_covalent']
        if sigma:
            if len(sigma)!=1 or ligand.ligand_type not in {'Ru_alkylidene','Schrock_alkylidene'}:
                raise Unsupported('Intramolecular pi binding requires one actual alkylidene ligand')
            ligand.ionic_class='L2' if ats[m].element=='Ru' else 'LX2'
        elif any(x.ligand_id==ligand.ligand_id for x in view.metal_bindings):raise Unsupported('Incoming ligand is bound elsewhere')
        ligand.donor_atom_ids=sorted(set(ligand.donor_atom_ids)|{a,b})
    else:
        ligand=Ligand(ligand_id='om.incoming.'+event_id,atom_ids=sorted(ids),donor_atom_ids=[a,b],
            ionic_class='L',ionic_charge=0,assignment_basis='Complete real organic component; one explicit C=C donor pool')
        view.ligands.append(ligand)
    view.metal_bindings.append(MetalBinding(binding_id='om.pi.'+event_id,metal_atom_id=m,ligand_id=ligand.ligand_id,
        donor_atom_ids=[a,b],mode='eta',hapticity=2,donated_electrons=2,ionic_allocation=0,pi_pool_slot=0))


def release_pi(view,m,a,b):
    matches=[x for x in view.metal_bindings if x.metal_atom_id==m and x.mode=='eta' and set(x.donor_atom_ids)=={a,b}]
    if len(matches)!=1:raise FormalError('Exactly one actual eta2 binding must be released')
    binding=matches[0];view.metal_bindings.remove(binding)
    ligand=next(l for l in view.ligands if l.ligand_id==binding.ligand_id)
    cov=[x for x in view.metal_bindings if x.ligand_id==ligand.ligand_id and x.mode=='sigma_covalent']
    if cov:
        ligand.donor_atom_ids=sorted({i for x in cov for i in x.donor_atom_ids})
        ligand.ionic_class='L' if atoms_by_id(view)[m].element=='Ru' else 'X2'


def release_donor(view,m,donor):
    matches=[x for x in view.metal_bindings if x.metal_atom_id==m and x.mode=='dative' and x.donor_atom_ids==[donor]]
    if len(matches)!=1:raise FormalError('Named dissociation donor is not uniquely coordinated')
    binding=matches[0];l=next(l for l in view.ligands if l.ligand_id==binding.ligand_id)
    a=atoms_by_id(view)[donor]
    if a.element not in {'P','N','O'} or a.formal_charge or a.radical_electrons!=0:
        raise Unsupported('Metathesis labile donor must have a resolved neutral P/N/O pair')
    view.metal_bindings.remove(binding)
    sigma=[x for x in view.metal_bindings if x.ligand_id==l.ligand_id and x.mode=='sigma_covalent']
    if sigma:
        if a.element!='O' or l.ligand_type!='Ru_alkylidene' or l.ionic_class!='L2':raise Unsupported('Only the registered Hoveyda ether chelate can open within an alkylidene ligand')
        l.ionic_class='L';l.donor_atom_ids=[i for x in sigma for i in x.donor_atom_ids]


def _reset(view,m):
    for c in view.metal_centers:
        if c.atom_id==m:c.oxidation_state=None;c.d_electrons=None;c.valence_electrons=None;c.relative_positions=[];c.vacant_site_ids=None
    for d in view.element_declarations:
        if d.atom_id==m:d.oxidation_state=None;d.coordination_number=None


def execute(source,event):
    validate_view(source);r=_roles(source,event);m=r['metal'];name=event.general_event
    before=checks(source)[1];actions=[]
    if name in {'OM_DONOR_DISSOCIATION','OM_OLEFIN_ASSOCIATION','OM_HOVEYDA_EXCHANGE_INITIATION','OM_PRODUCT_RELEASE'}:
        if event.actions:raise FormalError('Donor-owned coordination must not be counted as covalent electron transfer')
        target=source.model_copy(deep=True)
        if name in {'OM_DONOR_DISSOCIATION','OM_HOVEYDA_EXCHANGE_INITIATION'}:release_donor(target,m,r['donor'])
        if name in {'OM_OLEFIN_ASSOCIATION','OM_HOVEYDA_EXCHANGE_INITIATION'}:bind_pi(target,m,r['a'],r['b'],event.event_id)
        if name=='OM_PRODUCT_RELEASE':release_pi(target,m,r['a'],r['b'])
        target.state_id=event.product_state_id;target.view_id=event.product_view_id;_reset(target,m)
        validate_view(target)
    else:
        from .execution import transact
        from .metal_events import repartition
        edges=bond_orders(source)
        def bond(a,b,n):
            if edges.get(tuple(sorted((a,b))),0)!=n:raise FormalError(f'Metathesis expected bond {a}-{b} order {n}')
        if name=='OM_METALLACYCLE_FORMATION':
            c,a,b=r['alkylidene'],r['a'],r['b']
            for x,y,n in [(m,c,2),(a,b,2),(c,a,0),(m,b,0)]:bond(x,y,n)
            pi=[x for x in source.metal_bindings if x.metal_atom_id==m and x.mode=='eta' and set(x.donor_atom_ids)=={a,b}]
            if len(pi)!=1:raise FormalError('Chauvin formation requires a preceding explicit eta2 partner')
            actions=[Shift(a=m,b=c,c=a),Shift(a=a,b=b,c=m)]
            def post(v):
                release_pi(v,m,a,b);repartition(v);_reset(v,m)
            released={pi[0].binding_id}
        else:
            a,b,c=r['first'],r['second'],r['third']
            if any(atoms_by_id(source)[i].element!='C' for i in (a,b,c)):raise Unsupported('Metallacyclobutane has three mapped carbon atoms')
            for x,y in [(m,a),(a,b),(b,c),(c,m)]:bond(x,y,1)
            bond(a,c,0);bond(m,b,0)
            actions=[Shift(a=a,b=m,c=c),Shift(a=c,b=b,c=a)]
            def post(v):
                repartition(v)
                from ...libraries.organometallic.metathesis import decorate
                decorate(v,m,convention='neutral_L' if atoms_by_id(v)[m].element=='Ru' else 'X2')
                bind_pi(v,m,a,b,event.event_id);_reset(v,m)
            released=set()
        if event.actions is not None and sorted(a.model_dump_json() for a in event.actions)!=sorted(a.model_dump_json() for a in actions):
            raise FormalError('Declared metathesis arrows disagree with the registered grouped transaction')
        target=transact(source,actions,event.product_state_id,event.product_view_id,after=post,_released_binding_ids=released)
        event.actions=actions
    if identity(source)!=identity(target) or electron_count(source)!=electron_count(target) or sum(a.formal_charge for a in source.atoms)!=sum(a.formal_charge for a in target.atoms):
        raise FormalError('Metathesis material/electron/charge conservation failed')
    return target,{'event_id':event.event_id,'event_type':'POLAR','general_event':name,'extension_rule_version':VERSION,
        'formal_status':'passed','execution_mode':'electron_actions','atomic_commit':True,
        'compiled_actions':[a.model_dump(mode='json') for a in actions],
        'metal_counts_before':before,'metal_counts_after':checks(target)[1],
        'chemical_feasibility':'unknown','selectivity':'unknown','total_spin_inferred':False,
        'initiation_resolution':'grouped coordination interchange; not a computed rate law' if name=='OM_HOVEYDA_EXCHANGE_INITIATION' else 'registered finite event'}
