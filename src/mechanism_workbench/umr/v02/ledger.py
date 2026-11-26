from ..chemistry import (FormalError, Unsupported, atoms_by_id, bond_valences, derived_charge,
                         electron_count, validate_view as validate01)

# Explicit finite neutral local radical profiles: (bond valence, N, r).
RADICAL_PROFILES = {"C": (3, 1, 1), "O": (1, 5, 1), "N": (2, 3, 1), "H": (0, 1, 1),
                    "Cl": (0, 7, 1), "Br": (0, 7, 1), "I": (0, 7, 1), "S": (1, 5, 1)}


def radical_count(view):
    return sum(a.radical_electrons for a in view.atoms)


def validate_spin(view, multiplicity):
    if multiplicity is None:
        return "unknown"
    count = radical_count(view)
    if multiplicity not in range(count + 1, 0, -2):
        raise Unsupported("requires_spin_transition: supplied multiplicity is incompatible with local open-shell occupancy")
    return "compatible_local_occupancy_not_spin_dynamics"


def validate_view(view, state=None):
    valences = bond_valences(view)
    # A shadow validates the unchanged closed-shell atoms using the original rules.
    # Radical atoms are independently checked below, never relabelled in the real graph.
    shadow = view.model_copy(deep=True)
    radical_ids, checks = set(), []
    for atom in view.atoms:
        n, r = atom.nonbonding_electrons, atom.radical_electrons
        if n is None:
            raise Unsupported(f"Unresolved electron ledger at atom {atom.atom_id}")
        if n < r or (n-r) % 2:
            raise FormalError(f"N >= r and even N-r required at atom {atom.atom_id}")
        charge = derived_charge(atom, valences[atom.atom_id])
        if atom.formal_charge != charge:
            raise FormalError(f"Charge ledger mismatch at atom {atom.atom_id}")
        if not r and atom.element=="C" and charge==1 and (valences[atom.atom_id],n)!=(3,0):
            raise Unsupported("Carbon-cation profile outside classical three-coordinate carbenium; not a universal invalidity judgment")
        if r:
            if (valences[atom.atom_id], n, r) != RADICAL_PROFILES.get(atom.element) or charge != 0:
                raise Unsupported(f"Radical profile not implemented at atom {atom.atom_id}; radical ions are retained but unsupported")
            radical_ids.add(atom.atom_id)
            checks.append({"rule": f"neutral-{atom.element}-radical/0.2", "atom_id": atom.atom_id,
                           "shell_electrons": n+2*valences[atom.atom_id], "N": n, "r": r, "q": charge})
            substitute = atoms_by_id(shadow)[atom.atom_id]
            substitute.nonbonding_electrons += -1 if atom.element in {"C", "H"} else 1
            substitute.radical_electrons = 0
            substitute.formal_charge = derived_charge(substitute, valences[atom.atom_id])
    original_checks = validate01(shadow)
    charge = sum(a.formal_charge for a in view.atoms)
    if state:
        if state.state_id != view.state_id or state.total_charge != charge:
            raise FormalError("State/view identity or total charge mismatch")
        if radical_ids and state.electronic_mode != "RADICAL":
            raise FormalError("Open-shell state requires explicit electronic_mode=RADICAL")
        spin = validate_spin(view, state.multiplicity)
    else:
        spin = "unknown"
    return {"total_charge": charge, "electrons": electron_count(view), "unpaired_electrons": radical_count(view),
            "atom_checks": [c for c in original_checks["atom_checks"] if c["atom_id"] not in radical_ids] + checks,
            "spin_status": spin, "chemical_feasibility": "unknown"}
