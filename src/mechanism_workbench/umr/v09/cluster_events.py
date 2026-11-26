"""Finite Co2(CO)6 / mu-eta2:eta2 alkyne donor-owned convention."""
from ..chemistry import atoms_by_id,order,FormalError,Unsupported,identity,electron_count
from ..v07.ledger import project_symbols,symbols
from .ast import Ligand,MetalBinding,RULE_VERSION
from .ledger import validate_view
from .coordination import checks
SPECS={'DICOBALT_ALKYNE_ASSOCIATION':{'roles':['metal_a','metal_b','alkyne_a','alkyne_b'],
    'source':'Co2(CO)6 with explicit two distinct alkyne pi donation pools; no covalent/eta resonance assertion'}}
def execute(v,e):
    validate_view(v);r=e.general_roles;ats=atoms_by_id(v)
    if set(r)!=set(SPECS[e.general_event]['roles']) or len(set(r.values()))!=4 or not set(r.values())<=ats.keys():raise FormalError('Four distinct cluster/alkyne roles required')
    if e.event_type!='COORDINATION' or e.actions or e.electron_transfers or e.electronic_outcome:raise FormalError('Donor-owned coordination cannot contain covalent or redox actions')
    if not e.evidence or not e.mechanism_assumption:raise FormalError('Explicit cluster convention and evidence required')
    m,n,a,b=[r[k] for k in SPECS[e.general_event]['roles']]
    if [ats[i].element for i in [m,n,a,b]]!=['Co','Co','C','C'] or order(v,m,n)!=1 or order(v,a,b)!=3:raise Unsupported('Finite Co-Co/alkyne triple-bond input required')
    counts={x['metal_atom_id']:x for x in checks(v)[1]};ligs={l.ligand_id:l for l in v.ligands}
    for i in [m,n]:
        bs=[x for x in v.metal_bindings if x.metal_atom_id==i]
        if counts[i]['oxidation_state']!=0 or len(bs)!=3:raise Unsupported('Each Co0 needs exactly three bound CO ligands; explicitly release the other CO first')
        for bind in bs:
            ds=bind.donor_atom_ids;ids=ligs[bind.ligand_id].atom_ids
            if bind.mode!='dative' or len(ds)!=1 or len(ids)!=2 or ats[ds[0]].element!='C' or not any(ats[j].element=='O' and order(v,ds[0],j)==3 for j in ids):raise Unsupported('Only complete terminal CO ligands in this cluster profile')
    ids={a};todo=[a];metals=set(counts)
    while todo:
        i=todo.pop()
        for edge in v.bonds:
            if i not in {edge.a,edge.b}:continue
            j=edge.b if edge.a==i else edge.a
            if j in metals:raise Unsupported('Alkyne already has covalent metal binding')
            if j not in ids:ids.add(j);todo.append(j)
    if b not in ids or any(set(bind.donor_atom_ids)&ids for bind in v.metal_bindings):raise FormalError('Require one free complete alkyne component')
    target=v.model_copy(deep=True);target.state_id=e.product_state_id;target.view_id=e.product_view_id
    q=[ats[i].formal_charge for i in ids]
    opaque=symbols(v)&ids
    if any(ats[i].formal_charge is None for i in ids-opaque):raise Unsupported('Active alkyne ligand charge not resolved')
    if sum(ats[i].formal_charge for i in ids-opaque):raise Unsupported('Finite locally neutral alkyne ligand only')
    lid=e.event_id+'.alkyne'
    target.ligands.append(Ligand(ligand_id=lid,atom_ids=sorted(ids),donor_atom_ids=[a,b],ionic_class='L',ionic_charge=None if opaque else 0,
        assignment_basis='two independent pi-pair pools; intact C#C is the donor-owned carrier convention, not measured bond order'))
    for slot,metal in enumerate([m,n]):
        target.metal_bindings.append(MetalBinding(binding_id=f'{lid}.{metal}',metal_atom_id=metal,ligand_id=lid,donor_atom_ids=[a,b],mode='eta',hapticity=2,ionic_allocation=0,donated_electrons=2,pi_pool_slot=slot))
        center=next(x for x in target.metal_centers if x.atom_id==metal);center.valence_electrons=None;center.vacant_site_ids=None;center.relative_positions=[]
    ledger=validate_view(target);s,t=project_symbols(v),project_symbols(target)
    if identity(v)!=identity(target) or electron_count(s)!=electron_count(t):raise FormalError('Coordination must preserve atom and electron ownership')
    return target,{'event_id':e.event_id,'event_type':e.event_type,'general_event':e.general_event,'execution_mode':'electron_actions',
        'rule_version':RULE_VERSION,'formal_status':'passed','atomic_commit':True,'compiled_actions':[],
        'coordination_change':{'metals':[m,n],'donors':[a,b],'pi_pool_slots':[0,1]},'electron_ownership_changed':False,
        'chemical_feasibility':'unknown','ledger':ledger,'scope':'association prefix only; no inference of subsequent cluster insertion/decomplexation'}
