"""Metal-to-metal SET of known total N; unknown spin is never replaced by zero."""
from collections import Counter
from ..chemistry import FormalError,Unsupported,atoms_by_id,bond_valences,electron_count,identity
from .elements import charge, current_elements
from .ledger import validate_view
from .coordination import checks

def reduce_carbon_radical(source,event):
    """One owned metal electron fills one explicit neutral C radical slot."""
    validate_view(source)
    if event.event_type!='SET' or event.actions or len(event.electron_transfers)!=1 or not event.evidence or not event.mechanism_assumption:
        raise FormalError('Metal/radical SET needs one explicit transfer and evidence')
    tx=event.electron_transfers[0];d,a=tx.donor_atom_id,tx.acceptor_atom_id;ats=atoms_by_id(source)
    if d==a or not {d,a}<=ats.keys():raise FormalError('Distinct existing SET endpoints required')
    counts={x['metal_atom_id']:x for x in checks(source)[1]}
    if d not in counts or (ats[d].element,counts[d]['oxidation_state']) not in {('Cu',1),('Cr',2),('Sm',2)}:
        raise Unsupported('Finite CuI/CrII or explicitly ionic SmII donor scope for carbon radical reduction')
    if ats[d].element=='Sm' and (bond_valences(source)[d] or any(b.metal_atom_id==d for b in source.metal_bindings) or event.conditions.get('electron_transfer_model')!='outer_sphere_formal'):
        raise Unsupported('SmII extension requires isolated ionic ownership and explicit outer_sphere_formal declaration')
    if (ats[a].element,ats[a].formal_charge,ats[a].nonbonding_electrons,ats[a].radical_electrons,bond_valences(source)[a])!=('C',0,1,1,3):
        raise Unsupported('Explicit neutral tricoordinate one-electron carbon radical required')
    if source.electronic_state.character!='unknown' or source.electronic_state.coupling!='unknown':
        raise Unsupported('Defined coupled spin requires a separate transition model')
    if any((b.direction!='inside' or b.domain in {'surface','bulk_metal','external_electron_reservoir'}) and {d,a}&set(b.atom_ids) for b in source.material_boundaries):
        raise Unsupported('Material/reservoir electron exchange has not been specified')
    if ats[d].nonbonding_electrons<1:raise FormalError('Donor owns no electron')
    out=source.model_copy(deep=True);out.state_id=event.product_state_id;out.view_id=event.product_view_id
    atoms=atoms_by_id(out);atoms[d].nonbonding_electrons-=1;atoms[d].radical_electrons=None
    atoms[a].nonbonding_electrons+=1;atoms[a].radical_electrons=0
    vals=bond_valences(out)
    for i in (a,d):atoms[i].formal_charge=charge(atoms[i],vals[i])
    for c in out.metal_centers:
        if c.atom_id==d:c.oxidation_state=None;c.d_electrons=None;c.valence_electrons=None;c.spin_state='unknown'
    for c in out.element_declarations:
        if c.atom_id in {a,d}:c.oxidation_state=None;c.coordination_number=None
    result=validate_view(out)
    from ..v07.ledger import project_symbols
    s,t=project_symbols(source),project_symbols(out)
    if identity(source)!=identity(out) or electron_count(s)!=electron_count(t) or sum(x.formal_charge for x in s.atoms)!=sum(x.formal_charge for x in t.atoms):
        raise FormalError('Metal/radical SET conservation failed')
    after={x['metal_atom_id']:x for x in checks(out)[1]}
    if after[d]['oxidation_state']!=counts[d]['oxidation_state']+1:raise FormalError('Donor oxidation inconsistent with one-electron transfer')
    from .ast import RULE_VERSION
    return out,{'event_id':event.event_id,'event_type':'SET','resolved_profile':'METAL_C_RADICAL_REDUCTION','rule_version':RULE_VERSION,
        'execution_mode':'electron_actions','atomic_commit':True,'formal_status':'passed','compiled_actions':[],
        'electron_transfers':[tx.model_dump(mode='json')],'acceptor_slot':{'atom_id':a,'before_unpaired':1,'after_unpaired':0},
        'spin_status':'unknown_global_coupling_preserved','chemical_feasibility':'unknown','ledger':result}

def execute_dissociative(source,event):
    from ..chemistry import order
    from ..v02.electrons import bond_to_atom
    validate_view(source);r=event.general_roles;ats=atoms_by_id(source)
    if set(r)!={'metal','carbon','halogen'} or len(set(r.values()))!=3 or not set(r.values())<=ats.keys():raise FormalError('Distinct metal/C/X roles required')
    m,c,x=[r[k] for k in ['metal','carbon','halogen']]
    counts={a['metal_atom_id']:a for a in checks(source)[1]}
    if m not in counts or (ats[m].element,counts[m]['oxidation_state']) not in {('Cu',0),('Cr',2),('Li',0),('Sm',2)}:raise Unsupported('Finite Cu0/CrII/Li0 or explicit ionic SmII dissociative ownership channel; no bulk surface model')
    if ats[m].element=='Sm' and (bond_valences(source)[m] or any(b.metal_atom_id==m for b in source.metal_bindings) or event.conditions.get('electron_transfer_model')!='outer_sphere_formal'):
        raise Unsupported('SmII extension requires isolated ionic ownership and explicit outer_sphere_formal declaration')
    if ats[c].element!='C' or ats[x].element not in {'Cl','Br','I'} or order(source,c,x)!=1:raise FormalError('Explicit intact C-halogen bond required')
    if event.event_type!='SET' or not event.evidence or not event.mechanism_assumption:raise FormalError('Declared SET channel and evidence required')
    if event.actions or len(event.electron_transfers)!=1:raise FormalError('Ownership SET has explicit transfer; do not assert an unknown orbital slot')
    tx=event.electron_transfers[0]
    if (tx.donor_atom_id,tx.acceptor_atom_id)!=(m,x):raise FormalError('Transferred electron must end in the leaving-halide pair')
    if ats[m].nonbonding_electrons<1 or any(ats[i].radical_electrons!=0 for i in [c,x]):raise FormalError('Donor electron and closed-shell C-X bond required')
    if source.electronic_state.coupling not in {'unknown','unimplemented'}:raise Unsupported('Defined global spin requires a specific transition model')
    if any(b.direction!='inside' for b in source.material_boundaries):raise Unsupported('Open reservoir SET requires external exchange accounting')
    out=source.model_copy(deep=True);out.state_id=event.product_state_id;out.view_id=event.product_view_id
    out.bonds=[b for b in out.bonds if {b.a,b.b}!={c,x}];a=atoms_by_id(out)
    a[m].nonbonding_electrons-=1;a[m].radical_electrons=None
    a[c].nonbonding_electrons+=1;a[c].radical_electrons=1
    a[x].nonbonding_electrons+=2
    vals=bond_valences(out)
    for i in [m,c,x]:a[i].formal_charge=charge(a[i],vals[i])
    for center in out.metal_centers:
        if center.atom_id==m:center.oxidation_state=None;center.d_electrons=None;center.valence_electrons=None;center.spin_state='unknown'
    for d in out.element_declarations:
        if d.atom_id in {m,c,x}:d.oxidation_state=None;d.coordination_number=None
    check=validate_view(out)
    from ..v07.ledger import project_symbols
    s,t=project_symbols(source),project_symbols(out)
    if identity(source)!=identity(out) or electron_count(s)!=electron_count(t) or sum(a.formal_charge for a in s.atoms)!=sum(a.formal_charge for a in t.atoms):raise FormalError('Dissociative SET conservation failed')
    return out,{'event_id':event.event_id,'event_type':'SET','general_event':event.general_event,'execution_mode':'electron_actions','atomic_commit':True,
        'formal_status':'passed','compiled_actions':[bond_to_atom(c,x,'sigma',0,c).model_dump(mode='json'),bond_to_atom(c,x,'sigma',1,x).model_dump(mode='json')],
        'electron_transfers':[tx.model_dump(mode='json')],'spin_status':'metal_coupling_unknown','chemical_feasibility':'unknown','ledger':check}

def execute(source,event):
    validate_view(source)
    if event.event_type!='SET' or not event.evidence or not event.mechanism_assumption:raise FormalError('SET needs classification, source and declared pathway')
    if event.actions or len(event.electron_transfers)!=1:raise Unsupported('Finite SET profile uses exactly one typed ownership transfer')
    tx=event.electron_transfers[0];ats=atoms_by_id(source);d,a=tx.donor_atom_id,tx.acceptor_atom_id
    if d==a or not {d,a}<=ats.keys():raise FormalError('Distinct existing SET endpoints required')
    if ats[d].element=='C' and ats[a].element=='Cu':return oxidize_carbon_radical(source,event)
    if ats[d].element in current_elements().metals and ats[a].element=='C':return reduce_carbon_radical(source,event)
    if any(ats[i].element not in current_elements().metals for i in (d,a)):raise Unsupported('Ownership SET supports two metals or a registered metal/C-radical reduction')
    if source.electronic_state.character!='unknown' or source.electronic_state.coupling!='unknown':raise Unsupported('Definite coupled spin requires an explicit supported transition model')
    if ats[d].nonbonding_electrons<1:raise FormalError('SET donor has no owned electron')
    if any(b.direction!='inside' for b in source.material_boundaries):raise Unsupported('Open reservoir SET needs an external exchange ledger')
    target=source.model_copy(deep=True);target.state_id=event.product_state_id;target.view_id=event.product_view_id
    out=atoms_by_id(target);out[d].nonbonding_electrons-=1;out[a].nonbonding_electrons+=1
    vals=bond_valences(target)
    for i in (d,a):
        out[i].formal_charge=charge(out[i],vals[i]);out[i].radical_electrons=None
    for m in target.metal_centers:
        if m.atom_id in {d,a}:m.oxidation_state=None;m.d_electrons=None;m.valence_electrons=None;m.spin_state='unknown'
    for decl in target.element_declarations:
        if decl.atom_id in {d,a}:decl.oxidation_state=None
    result=validate_view(target)
    from ..v07.ledger import project_symbols
    s,t=project_symbols(source),project_symbols(target)
    if identity(source)!=identity(target) or electron_count(s)!=electron_count(t) or sum(x.formal_charge for x in s.atoms)!=sum(x.formal_charge for x in t.atoms):raise FormalError('SET conservation failed')
    before={r['metal_atom_id']:r['oxidation_state'] for r in checks(source)[1]};after={r['metal_atom_id']:r['oxidation_state'] for r in checks(target)[1]}
    if after[d]-before[d]!=1 or after[a]-before[a]!=-1:raise FormalError('SET ownership and oxidation changes disagree')
    return target,{'event_id':event.event_id,'event_type':'SET','general_event':'METAL_SET','atomic_commit':True,
       'execution_mode':'electron_actions','formal_status':'passed','compiled_actions':[],
       'electron_transfers':[tx.model_dump(mode='json')],'spin_status':'unknown_coupling_not_inferred',
       'oxidation_before':before,'oxidation_after':after,'chemical_feasibility':'unknown','ledger':result}

def oxidize_carbon_radical(source,event):
    """Finite outer-sphere CuII oxidation; no guessed metal orbital slot."""
    tx=event.electron_transfers[0];d,a=tx.donor_atom_id,tx.acceptor_atom_id;ats=atoms_by_id(source)
    if event.conditions.get('electron_transfer_model')!='outer_sphere_formal':raise Unsupported('Explicit outer-sphere declaration required')
    if (ats[a].element,ats[a].formal_charge,bond_valences(source)[a])!=('Cu',2,0) or any(b.metal_atom_id==a for b in source.metal_bindings):raise Unsupported('Isolated ionic CuII acceptor only')
    if (ats[d].element,ats[d].formal_charge,ats[d].nonbonding_electrons,ats[d].radical_electrons,bond_valences(source)[d])!=('C',0,1,1,3):raise Unsupported('Neutral tricoordinate carbon radical donor required')
    if source.electronic_state.character!='unknown' or source.electronic_state.coupling!='unknown':raise Unsupported('Defined global coupling requires explicit transition model')
    out=source.model_copy(deep=True);out.state_id=event.product_state_id;out.view_id=event.product_view_id
    atoms=atoms_by_id(out);atoms[d].nonbonding_electrons-=1;atoms[d].radical_electrons=0
    atoms[a].nonbonding_electrons+=1;atoms[a].radical_electrons=None
    vals=bond_valences(out)
    for i in [d,a]:atoms[i].formal_charge=charge(atoms[i],vals[i])
    for c in out.metal_centers:
        if c.atom_id==a:c.oxidation_state=None;c.d_electrons=None;c.valence_electrons=None;c.spin_state='unknown'
    for c in out.element_declarations:
        if c.atom_id in {a,d}:c.oxidation_state=None;c.coordination_number=None
    check=validate_view(out)
    if identity(source)!=identity(out) or electron_count(source)!=electron_count(out) or sum(x.formal_charge for x in source.atoms)!=sum(x.formal_charge for x in out.atoms):raise FormalError('CuII radical oxidation conservation failed')
    return out,{'event_id':event.event_id,'event_type':'SET','resolved_profile':'CU_II_C_RADICAL_OXIDATION','execution_mode':'electron_actions','atomic_commit':True,'formal_status':'passed','compiled_actions':[],
        'electron_transfers':[tx.model_dump(mode='json')],'spin_status':'unknown_global_coupling_preserved','chemical_feasibility':'unknown','ledger':check}
