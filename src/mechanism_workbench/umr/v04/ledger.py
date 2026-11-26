"""Finite local occupation checks; no ground-state or quantum-wavefunction inference."""
import json

from ..chemistry import (FormalError, Unsupported, atoms_by_id, bond_valences, derived_charge,
                         electron_count, graph_signature, order)
from ..v02.ledger import validate_view as validate02
from ..v03.ledger import compare_views as compare03, signature as signature03, validate_view as validate03

# species, element, b, N, r; N includes every unpaired electron.
PROFILES = {
    "carbene_paired": ("carbene", "C", 2, 2, 0),
    "carbene_two_singly": ("carbene", "C", 2, 2, 2),
    "carbyne_one_singly": ("carbyne", "C", 1, 3, 1),
    "carbyne_three_singly": ("carbyne", "C", 1, 3, 3),
    "nitrene_paired": ("nitrene", "N", 1, 4, 0),
    "nitrene_two_singly": ("nitrene", "N", 1, 4, 2),
}


def occupation_profiles(view):
    valences = bond_valences(view)
    return {a.atom_id: name for a in view.atoms for name, (_, element, b, n, r) in PROFILES.items()
            if (a.element, valences[a.atom_id], a.nonbonding_electrons, a.radical_electrons) == (element,b,n,r)}


def electronic_signature(view):
    return json.dumps({"electronic_state":view.electronic_state.model_dump(mode="json"),
        "centers":sorted((c.model_dump(mode="json", exclude={"center_id", "evidence", "existence_basis", "reaction_context"}) for c in view.reactive_centers),key=lambda c:c["atom_id"]),
        "arynes":sorted((r.model_dump(mode="json", exclude={"region_id", "evidence", "existence_basis", "reaction_context"}) for r in view.aryne_regions),key=lambda r:r["ring_atom_ids"])},sort_keys=True)


def signature(view):
    return signature03(view), electronic_signature(view)


def compare_views(actual, expected):
    compare03(actual,expected)
    if electronic_signature(actual) != electronic_signature(expected):
        raise FormalError("Electronic state/local occupation/aryne declaration differs from target")


def check_isolated_ring(view, ring):
    """This release allows substituents, but not fused/bridged aryne rings."""
    members=set(ring)
    edges={frozenset((i,ring[(j+1)%len(ring)])) for j,i in enumerate(ring)}
    adjacent={a.atom_id:set() for a in view.atoms}
    for bond in view.bonds:
        adjacent[bond.a].add(bond.b);adjacent[bond.b].add(bond.a)
        if {bond.a,bond.b}<=members and frozenset((bond.a,bond.b)) not in edges:
            raise Unsupported("Fused/bridged/chorded aryne rings need a dedicated scope")
    seen=set()
    for start in set(adjacent)-members:
        if start in seen: continue
        todo=[start];attached=set()
        while todo:
            i=todo.pop()
            if i in seen: continue
            seen.add(i);attached|=adjacent[i]&members;todo+=list(adjacent[i]-members-seen)
        if len(attached)>1:
            raise Unsupported("Fused/bridged aryne systems require a dedicated template")


def validate_arynes(view):
    for region in view.aryne_regions:
        if region.representation != "ortho_formal_triple":
            raise Unsupported("Only six-membered 1,2-aryne formal-triple regions are implemented")
        ring=region.ring_atom_ids; a,b=region.aryne_atom_ids; atoms=atoms_by_id(view)
        check_isolated_ring(view,ring)
        if any(atoms[i].element!="C" for i in ring):
            raise Unsupported("Heteroarynes are outside the registered six-carbon-ring scope")
        if order(view,a,b)!=3:
            raise FormalError("Declared aryne requires the specified formal triple bond")
        triples=0; doubles=0
        for j,i in enumerate(ring):
            k=ring[(j+1)%6]; o=order(view,i,k)
            if not o: raise FormalError("Aryne ring ordering is not a connected cycle")
            triples+=o==3; doubles+=o==2
            if atoms[i].formal_charge or atoms[i].nonbonding_electrons or atoms[i].radical_electrons:
                raise FormalError("Formal-triple aryne ring requires neutral ordinary tetravalent carbon")
        if triples!=1 or doubles!=2 or not any({ring[j],ring[(j+1)%6]}=={a,b} for j in range(6)):
            raise FormalError("Registered aryne needs one triple and two residual double bonds in the ring")
        if atoms[a].hydrogens or atoms[b].hydrogens or any(x.element=="H" and (order(view,a,x.atom_id) or order(view,b,x.atom_id)) for x in view.atoms):
            raise FormalError("Aryne active carbons cannot retain hydrogen")


def unannotated_ring_triples(view):
    """Flag six-membered cyclic C#C candidates; do not infer their electronic state."""
    atoms=atoms_by_id(view);adj={i:set() for i in atoms}
    for b in view.bonds: adj[b.a].add(b.b);adj[b.b].add(b.a)
    declared={frozenset(r.aryne_atom_ids) for r in view.aryne_regions}
    for b in view.bonds:
        if b.order!=3 or frozenset((b.a,b.b)) in declared: continue
        paths=[[b.a]]
        while paths:
            path=paths.pop()
            if len(path)==6:
                if path[-1]==b.b and all(atoms[i].element=="C" for i in path): return True
                continue
            for j in adj[path[-1]]-set(path):
                if len(path)==1 and j==b.b: continue
                paths.append(path+[j])
    return False


def validate_ledger(view):
    """Transaction-local graph validation independent of stale descriptive labels."""
    if any(a.radical_electrons is None or a.nonbonding_electrons is None for a in view.atoms):
        raise Unsupported("Unknown local electronic occupancy; preserve input, do not infer spin from SMILES")
    profiles=occupation_profiles(view)
    if view.metal_centers or view.metal_bindings or view.ligands:
        if profiles or view.reactive_centers or view.aryne_regions:
            raise Unsupported("Metal-bound/reactive-center mixed systems require a dedicated adapter")
        return validate03(view)
    shadow=view.model_copy(deep=True); valences=bond_valences(view); checks=[]
    for atom in view.atoms:
        if atom.formal_charge != derived_charge(atom,valences[atom.atom_id]):
            raise FormalError(f"Charge ledger mismatch at atom {atom.atom_id}")
        if atom.nonbonding_electrons<atom.radical_electrons or (atom.nonbonding_electrons-atom.radical_electrons)%2:
            raise FormalError(f"N must contain r and N-r must be even at {atom.atom_id}")
        if atom.atom_id in profiles:
            checks.append({"rule":profiles[atom.atom_id]+"/0.4","atom_id":atom.atom_id,"b":valences[atom.atom_id],
                "N":atom.nonbonding_electrons,"r":atom.radical_electrons,"q":atom.formal_charge,
                "shell_electrons":atom.nonbonding_electrons+2*valences[atom.atom_id]})
            # Private validation cap, following v0.2/v0.3 shadow validation. Never
            # enters the authoritative state, global ledger, export or drawing.
            cap=atoms_by_id(shadow)[atom.atom_id]
            cap.element="C"; cap.hydrogens+=4-valences[atom.atom_id]
            cap.nonbonding_electrons=cap.radical_electrons=cap.formal_charge=0
    old=validate02(shadow)
    return {"total_charge":sum(a.formal_charge for a in view.atoms),"electrons":electron_count(view),
        "unpaired_electrons":sum(a.radical_electrons for a in view.atoms),
        "atom_checks":[c for c in old["atom_checks"] if c["atom_id"] not in profiles]+checks,
        "chemical_feasibility":"unknown"}


def validate_view(view,state=None,*,_ledger=validate_ledger):
    for item in [*view.reactive_centers,*view.aryne_regions]:
        if item.existence_basis=="unknown" or item.reaction_context.strip() in {"", "unknown"} or not item.evidence:
            raise Unsupported("Intermediate existence needs source/reaction context; a name tag is not evidence")
        if any(e.source.strip() in {"", "unknown"} or not e.note.strip() for e in item.evidence):
            raise FormalError("Intermediate evidence requires a traceable source and applicability note")
    for center in view.reactive_centers:
        if center.environment!="free":
            raise Unsupported("Metal-bound or unspecified carbene/carbyne/nitrene execution is unsupported")
    ledger=_ledger(view); profiles=occupation_profiles(view)
    if set(profiles)!={c.atom_id for c in view.reactive_centers}:
        raise FormalError("Low-valent occupations require exact reactive-center declarations")
    for c in view.reactive_centers:
        if c.profile=="unknown": raise Unsupported("Local profile is explicitly unknown")
        if profiles.get(c.atom_id)!=c.profile or PROFILES[c.profile][0]!=c.species:
            raise FormalError("Reactive-center tag contradicts structure/electron occupancy")
    validate_arynes(view)
    if unannotated_ring_triples(view):
        raise Unsupported("Cyclic C#C candidate requires explicit aryne region/electronic-state evidence; ordinary alkyne templates cannot stand in")
    r=ledger['unpaired_electrons'] if 'unpaired_electrons' in ledger else sum(a.radical_electrons for a in view.atoms)
    if r is None:raise Unsupported('Unknown global spin occupation; local ledger retained')
    character=view.electronic_state.character
    if character=="closed_shell" and r: raise FormalError("Closed-shell declaration conflicts with occupied single-electron slots")
    if character in {"open_shell","open_shell_singlet"} and not r:
        raise FormalError("Open-shell declaration has no singly occupied slots")
    if character=="open_shell_singlet":
        raise Unsupported("Open-shell singlet declaration retained; spin-coupled occupation model unsupported")
    if view.electronic_state.coupling=="closed_shell" and r:
        raise FormalError("Closed-shell coupling conflicts with occupancy")
    if state:
        if view.state_id!=state.state_id or state.total_charge!=ledger["total_charge"]:
            raise FormalError("State/view identity or total charge mismatch")
        if view.electronic_state!=state.electronic_state:
            raise FormalError("Views of a state must share its electronic-state declaration")
        m=state.multiplicity
        if m is not None and m not in range(r+1,0,-2):
            raise Unsupported("requires_spin_transition: total multiplicity incompatible with local occupancy")
        if m==1 and r:
            raise Unsupported("Open-shell singlet retained; coupling execution unsupported")
        if m is not None and r>1 and m!=r+1:
            raise Unsupported("Non-high-spin multi-center coupling retained but execution unsupported")
        if state.electronic_mode=="CLOSED_SHELL" and r:
            raise FormalError("State electronic_mode conflicts with occupancy")
        coupling=view.electronic_state.coupling
        if m is not None and coupling=="high_spin" and m!=r+1:
            raise FormalError("Declared high-spin coupling contradicts total multiplicity")
        if coupling=="one_unpaired" and (r!=1 or m not in {None,2}):
            raise FormalError("one_unpaired coupling needs exactly one unpaired electron")
    ledger["spin_status"]="declared_compatibility_only_not_ground_state"
    ledger["intermediate_existence"]="provenance_present_not_independently_verified"
    return ledger
