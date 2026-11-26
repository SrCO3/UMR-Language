import json
from ..chemistry import FormalError,Unsupported,bond_valences,atoms_by_id,electron_count
from ..v05.ledger import validate_ledger as shared_ledger
from ..v04.ledger import validate_view as reactive_checks
from ..v07.ledger import symbols,project_symbols,signature as old_signature,compare_views as old_compare
from .elements import LEGACY,current_elements,profile,charge
from .coordination import checks
from ..v08.elements import oxidation
from .local_rules import active_rules, active_reagents, check_local_rules, check_reagents, require_equal

def validate_view(view,state=None,*,_symbolic_ids=frozenset()):
    elements = current_elements()
    METALS, MAIN = elements.metals, elements.main
    boundary_ids=set();seen_atoms=set();known={a.atom_id for a in view.atoms}
    for boundary in view.material_boundaries:
        if boundary.boundary_id in boundary_ids or len(set(boundary.atom_ids))!=len(boundary.atom_ids) or not boundary.atom_ids or not set(boundary.atom_ids)<=known:
            raise FormalError('Duplicate/absent material boundary reference')
        if seen_atoms&set(boundary.atom_ids):raise FormalError('Atom belongs to conflicting material boundaries')
        boundary_ids.add(boundary.boundary_id);seen_atoms.update(boundary.atom_ids)
    if symbols(view):
        shadow=project_symbols(view);result=validate_view(shadow,_symbolic_ids=symbols(view))
        if state and (state.state_id!=view.state_id or state.electronic_state!=view.electronic_state):raise FormalError('Symbolic state/view reference mismatch')
        result.update(symbolic_boundary=True,total_ledger='unknown',electrons=None,total_charge=None,
          full_electron_total_verified=False,validation_scope='symbolic_unchanged_boundary',opaque_symbols=sorted(symbols(view)));return result
    vals=bond_valences(view);shadow=view.model_copy(deep=True)
    rules = active_rules(view.view_id)
    reagents = active_reagents(view.view_id)
    if reagents:check_reagents(view, reagents, _symbolic_ids)
    local = check_local_rules(view, rules, symbolic_ids=_symbolic_ids) if rules else {}
    for center in shadow.reactive_centers:
        if center.environment=='metal_bound':
            bs=[b for b in view.metal_bindings if center.atom_id in b.donor_atom_ids and b.mode=='dative']
            ls={l.ligand_id:l for l in view.ligands}
            if len(bs)!=1 or ls[bs[0].ligand_id].ligand_type not in {'Fischer_carbene','NHC'} or center.profile!='carbene_paired':
                raise Unsupported('Metal-bound low-valent center needs a registered donor-owned Fischer/NHC ligand')
            # Reuse only the local occupancy/evidence check on a private copy.
            center.environment='free'
    unknown=[];known_r=0
    for a in view.atoms:
        if a.nonbonding_electrons is None:raise Unsupported('Unresolved nonbonding electron ownership at '+str(a.atom_id))
        if a.formal_charge!=charge(a,vals[a.atom_id]):raise FormalError('Charge ownership mismatch at '+str(a.atom_id))
        r=a.radical_electrons
        if r is None:
            if a.element not in METALS:raise Unsupported('Unknown main-group local occupation')
            unknown.append(a.atom_id);atoms_by_id(shadow)[a.atom_id].radical_electrons=0
        else:
            known_r+=r
            if r>a.nonbonding_electrons or (a.nonbonding_electrons-r)%2:raise FormalError('Inconsistent single-electron slots')
            if r and a.formal_charge and a.element in {'C','N','O'} and a.atom_id not in local:
                if (a.element,vals[a.atom_id],a.nonbonding_electrons,r,a.formal_charge) not in {
                    ('C',2,1,1,1),('N',3,1,1,1),('O',2,3,1,1),('C',2,3,1,-1)}:
                    raise Unsupported('Charged-radical local profile not registered')
                cap=atoms_by_id(shadow)[a.atom_id];cap.nonbonding_electrons+=1;cap.radical_electrons=0;cap.formal_charge-=1
    # Unknown coupling and open-shell singlet are retained in authority. Only
    # local occupancy is checked on the private shadow; no ground state inferred.
    shadow.electronic_state.character='unknown';shadow.electronic_state.coupling='unknown'
    reserved,counts=checks(view,_symbolic_ids=_symbolic_ids)
    original=atoms_by_id(view)
    def prof(v,a):
        rule=local.get(a.atom_id)
        if rule is None:return profile(view,original[a.atom_id])
        metal = a.element in METALS
        if metal and a.hydrogens:raise FormalError('Metal hydrides need an explicit atom ID')
        if not metal and rule.occupied_sites is not None:
            raise FormalError(f'{rule.rule_id}: occupied_sites requires a metal binding model')
        if not metal and rule.oxidation_state is not None:
            require_equal(rule,'oxidation_state',rule.oxidation_state,oxidation(view,a.atom_id))
        return dict(atom_id=a.atom_id,element=rule.element,bond_order_sum=rule.bond_order_sum,
                    N=rule.nonbonding_electrons,r=rule.radical_electrons,formal_charge=rule.formal_charge,
                    shell_electrons=rule.nonbonding_electrons+2*rule.bond_order_sum,rule=rule.rule_id)
    special_radicals=[]
    for a in shadow.atoms:
        if a.element in {'Si','Sn','Se'} and a.radical_electrons and a.atom_id not in local:
            special_radicals.append(profile(view,original[a.atom_id]));a.element='O' if a.element=='Se' else 'C'
    local_ids = set(local) | {a.atom_id for a in view.atoms if a.element in elements.new_local}
    # Oxidation state is derived on the submitted graph. A charged radical's
    # private, closed-shell cap changes its formal charge (e.g. N in NO2).
    derive_os = oxidation if elements is LEGACY else lambda _, i: oxidation(view, i)
    report=shared_ledger(shadow,_profile=prof,_coordination=lambda _: (reserved,counts),_main=MAIN,_metals=METALS,_oxidation=derive_os,_local_atoms=local_ids)
    if elements is not LEGACY and view.element_declarations:
        # The original organic checker emits rows only for special centres.
        # Return the checked oxidation declaration for ordinary sites as well.
        checked = {row["atom_id"]: row for row in report["atom_checks"]}
        for declaration in view.element_declarations:
            i = declaration.atom_id
            atom = original[i]
            if atom.element not in {"H", "C", "N", "O", "F"}:
                continue
            if i not in checked:
                row = {"atom_id": i, "element": atom.element, "rule": "declared_Lewis_site/1.0"}
                report["atom_checks"].append(row)
                checked[i] = row
            checked[i].update(oxidation_state=oxidation(view, i), formal_charge=atom.formal_charge,
                              bond_order_sum=vals[i], N=atom.nonbonding_electrons, r=atom.radical_electrons,
                              shell_electrons=2*vals[i]+atom.nonbonding_electrons)
            if "q" in checked[i]:
                checked[i]["q"] = atom.formal_charge
    if special_radicals:report['atom_checks']=[x for x in report['atom_checks'] if x['atom_id'] not in {s['atom_id'] for s in special_radicals}]+special_radicals
    reactive_checks(shadow,_ledger=lambda _:report)
    character=view.electronic_state.character;coupling=view.electronic_state.coupling
    if character=='closed_shell' and known_r:raise FormalError('Closed-shell assertion contradicts local radicals')
    if character=='open_shell_singlet' and (unknown or known_r<2 or known_r%2):raise Unsupported('Singlet coupling cannot be verified with this local occupancy')
    if coupling=='closed_shell' and known_r:raise FormalError('Closed-shell coupling contradicts radicals')
    if coupling=='one_unpaired' and (unknown or known_r!=1):raise Unsupported('One-unpaired coupling not established by supplied occupation')
    if coupling=='unimplemented':raise Unsupported('Declared spin coupling is explicitly unsupported')
    if state:
        if state.state_id!=view.state_id or state.total_charge!=sum(a.formal_charge for a in view.atoms):raise FormalError('State charge/reference mismatch')
        if state.electronic_state!=view.electronic_state:raise FormalError('State/view electronic declarations differ')
        if state.multiplicity is not None and not unknown and state.multiplicity not in range(known_r+1,0,-2):raise FormalError('Multiplicity incompatible with local occupation')
        if coupling=='high_spin' and state.multiplicity is not None and not unknown and state.multiplicity!=known_r+1:raise FormalError('High-spin declaration conflicts with multiplicity')
    report.update(total_charge=sum(a.formal_charge for a in view.atoms),electrons=electron_count(view),
        unpaired_electrons=None if unknown else known_r,spin_status='unknown_metal_coupling' if unknown else 'occupation_compatible_not_wavefunction',rule_version=elements.version)
    if rules:report['local_rule_ids']=[r.rule_id for r in rules]
    if reagents:report['reagent_ids']=[r.reagent_id for r in reagents]
    return report
def graph_identity_view(view):
    # Published resonance-node identifiers depend on this serialisation. The 1.0
    # rule table changes acceptance, while the representation of a Lewis graph
    # and its declarations stays compatible with 0.9. Preserve existing node IDs.
    if view.element_rule_version == 'umr-elements/1.0':
        return view.model_copy(update={'element_rule_version':'umr-elements/0.9'})
    return view


def signature(v):return old_signature(graph_identity_view(v)),json.dumps({'ligands':[l.model_dump(mode='json') for l in v.ligands],'bindings':[b.model_dump(mode='json') for b in v.metal_bindings],'boundaries':[b.model_dump(mode='json') for b in v.material_boundaries]},sort_keys=True)
def compare_views(a,b):
    old_compare(graph_identity_view(a),graph_identity_view(b))
    if signature(a)!=signature(b):raise FormalError('v0.9 coordination/boundary differs from target')
    if active_rules(b.view_id) or active_reagents(b.view_id):validate_view(b)
