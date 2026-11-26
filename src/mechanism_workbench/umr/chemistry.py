"""Integer Lewis electron accounting. No editable charge delta exists."""

from collections import Counter

from .ast import Atom, ChemicalState, LewisView

VALENCE = {"H": 1, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7,
           "P": 5, "S": 6, "Cl": 7, "Br": 7, "I": 7, "Li": 1, "Na": 1, "K": 1}


class FormalError(ValueError):
    pass


class Unsupported(FormalError):
    pass


def atoms_by_id(view: LewisView) -> dict[int, Atom]:
    return {a.atom_id: a for a in view.atoms}


def bond_orders(view: LewisView) -> dict[tuple[int, int], int]:
    return {tuple(sorted((b.a, b.b))): b.order for b in view.bonds}


def order(view: LewisView, a: int, b: int) -> int:
    return bond_orders(view).get(tuple(sorted((a, b))), 0)


def bond_valences(view: LewisView) -> Counter:
    result = Counter({a.atom_id: a.hydrogens for a in view.atoms})
    for bond in view.bonds:
        result[bond.a] += bond.order
        result[bond.b] += bond.order
    return result


def derived_charge(atom: Atom, valence: int) -> int:
    if atom.element not in VALENCE or atom.nonbonding_electrons is None:
        raise Unsupported(f"Unknown or unsupported electron ledger at atom {atom.atom_id}: {atom.element}")
    return VALENCE[atom.element] - atom.nonbonding_electrons - valence


def carbocation_rule(atom: Atom, valence: int) -> dict:
    """Classical two-centre carbenium centres, not nonclassical carbonium ions."""
    occupancy = atom.nonbonding_electrons + 2 * valence
    valid = (atom.element == "C" and atom.formal_charge == 1 and valence == 3
             and atom.nonbonding_electrons == 0 and atom.radical_electrons == 0 and occupancy == 6)
    return {"rule": "classical-carbenium-sextet/0.1", "atom_id": atom.atom_id,
            "valid": valid, "shell_electrons": occupancy,
            "chemical_stability": "unknown", "orbital_geometry": "not_verified"}


def validate_view(view: LewisView, state: ChemicalState | None = None) -> dict:
    valences = bond_valences(view)
    checks = []
    if state and state.multiplicity != 1:
        raise Unsupported("Only closed-shell singlet states are implemented")
    for atom in view.atoms:
        if atom.radical_electrons or (atom.nonbonding_electrons is not None and atom.nonbonding_electrons % 2):
            raise Unsupported(f"Single-electron/radical species at atom {atom.atom_id}")
        q = derived_charge(atom, valences[atom.atom_id])
        if atom.formal_charge != q:
            raise FormalError(f"Charge ledger mismatch at atom {atom.atom_id}: stored={atom.formal_charge}, computed={q}")
        shell = atom.nonbonding_electrons + 2 * valences[atom.atom_id]
        if atom.element in {"Li", "Na", "K"}:
            if valences[atom.atom_id] != 0 or q != 1 or atom.nonbonding_electrons != 0:
                raise Unsupported("Only disconnected closed-shell alkali spectator cations are supported")
        elif atom.element == "H":
            if atom.hydrogens or shell not in {0, 2} or valences[atom.atom_id] > 1:
                raise FormalError(f"Invalid hydrogen valence at atom {atom.atom_id}")
        elif atom.element == "C" and q == 1:
            report = carbocation_rule(atom, valences[atom.atom_id])
            if not report["valid"]:
                raise FormalError(f"Carbon cation {atom.atom_id} must have bond-order valence 3, no lone electrons and six shell electrons")
            checks.append(report)
        elif atom.element == "B" and shell == 6 and q == 0:
            checks.append({"rule": "trivalent-boron/0.1", "atom_id": atom.atom_id, "shell_electrons": 6})
        elif shell > 8:
            if atom.element in {"P", "S", "Cl", "Br", "I"}:
                raise Unsupported("Expanded-valence/coordination species need a future adapter")
            raise FormalError(f"Octet exceeded at atom {atom.atom_id}: {shell}")
        elif shell != 8:
            raise Unsupported(f"Closed-shell under-octet centre outside the implemented rules at atom {atom.atom_id}")
    charge = sum(a.formal_charge for a in view.atoms)
    if state and (view.state_id != state.state_id or charge != state.total_charge):
        raise FormalError("View/state identity or total charge mismatch")
    return {"total_charge": charge, "electrons": electron_count(view), "atom_checks": checks,
            "chemical_feasibility": "unknown"}


def electron_count(view: LewisView) -> int:
    return sum(a.nonbonding_electrons for a in view.atoms) + 2 * sum(b.order for b in view.bonds) + 2 * sum(a.hydrogens for a in view.atoms)


def identity(view: LewisView) -> dict:
    return {a.atom_id: (a.element, a.isotope, a.hydrogens, a.label) for a in view.atoms}


def hydrogen_connections(view: LewisView) -> set:
    atoms = atoms_by_id(view)
    return {tuple(sorted((b.a, b.b))) for b in view.bonds
            if atoms[b.a].element == "H" or atoms[b.b].element == "H"}


def graph_signature(view: LewisView) -> tuple:
    return (tuple(sorted((a.atom_id, a.element, a.isotope, a.hydrogens, a.formal_charge,
                          a.nonbonding_electrons, a.radical_electrons, a.label) for a in view.atoms)),
            tuple(sorted((min(b.a, b.b), max(b.a, b.b), b.order) for b in view.bonds)))


def compare_views(actual: LewisView, expected: LewisView) -> None:
    if identity(actual) != identity(expected):
        raise FormalError("Atom identity/mapping/isotope/attached hydrogen mismatch (including a lost H)")
    if graph_signature(actual) != graph_signature(expected):
        raise FormalError("Mapped Lewis graph differs from the asserted target view")
    for atom in actual.atoms:
        target = atoms_by_id(expected)[atom.atom_id]
        if target.stereo.get("status") == "specified" and atom.stereo != target.stereo:
            raise Unsupported("Target asserts stereochemistry that the executor has not established")
    actual_bonds = {tuple(sorted((b.a, b.b))): b for b in actual.bonds}
    for bond in expected.bonds:
        if bond.stereo.get("status") == "specified" and actual_bonds[tuple(sorted((bond.a, bond.b)))].stereo != bond.stereo:
            raise Unsupported("Target asserts bond stereochemistry that the executor has not established")
