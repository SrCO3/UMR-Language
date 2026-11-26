"""UMR 1.0 element rules. Oxidation states and Lewis occupations are distinct."""
import csv
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType

from ..chemistry import FormalError, Unsupported, atoms_by_id, bond_valences
from ..v07.elements import oxidation
from ..v09 import elements as legacy
from .ast import LewisView


@dataclass(frozen=True)
class Reference:
    number: int
    block: str
    group: int
    common: frozenset[int]
    states: frozenset[int]


def read_references():
    result = {}
    with Path(__file__).with_name("oxidation-states.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            common = frozenset(int(v) for key in ("mendeleev_main", "pubchem", "rsc_common")
                               for v in row[key].split())
            extended = frozenset(map(int, row["mendeleev_extended"].split()))
            result[row["symbol"]] = Reference(int(row["Z"]), row["block"], int(row["group"]),
                                               common, common | extended | {0})
    return MappingProxyType(result)


REFERENCES = read_references()
LOCAL = frozenset({"Ge", "As", "Sb", "Te", "Xe"})
NEW_METALS = frozenset(REFERENCES.keys() - legacy.VALENCE.keys() - LOCAL)


def ownership_basis(ref):
    return ref.group - 10 if ref.block == "p" else ref.group


METALS = MappingProxyType({**legacy.METALS, **{s: ownership_basis(REFERENCES[s]) for s in NEW_METALS}})
VALENCE = MappingProxyType({**legacy.VALENCE, **METALS, **{s: ownership_basis(REFERENCES[s]) for s in LOCAL}})
OXIDATION = MappingProxyType({
    s: frozenset(legacy.OXIDATION.get(s, ())) | (REFERENCES[s].states if s in REFERENCES else frozenset())
    for s in METALS
})


def metal_counts(symbol, state):
    if symbol not in REFERENCES:  # Pre-existing lanthanide rules.
        return legacy.legacy_metal_counts(symbol, state)
    residual = METALS[symbol] - state
    # A formal group-minus-OS count above ten cannot be assigned to one d shell.
    d_electrons = residual if REFERENCES[symbol].block == "d" and 0 <= residual <= 10 else None
    return d_electrons, residual


# Each triple is (bond-order sum including H, nonbonding electrons, unpaired electrons).
# Charge is independently reconstructed from electron ownership by the shared ledger.
OCTET = {(b, 8 - 2 * b, 0) for b in range(5)}
PNICTOGEN = OCTET | {(4, 2, 0), (5, 0, 0), (6, 0, 0), (0, 5, 3)}
LOCAL_OCCUPATIONS = {
    "Ge": OCTET | {(2, 2, 0), (3, 0, 0), (3, 1, 1), (5, 0, 0), (6, 0, 0), (0, 4, 2)},
    "As": PNICTOGEN,
    "Sb": PNICTOGEN,
    "Te": OCTET | {(1, 5, 1), (4, 2, 0), (5, 2, 0), (6, 0, 0), (0, 6, 2)},
}
# (number of terminal fluorides, number of terminal double-bonded oxygens).
XENON_FORMS = {(0, 0), (2, 0), (4, 0), (6, 0), (2, 1), (4, 1), (2, 2), (0, 3), (0, 4)}

# Named closed-shell structures supplement the inherited Si, Sn and P rules.
# Each key gives element, bond-order sum, nonbonding electrons and formal charge.
# Nonempty entries require equivalent terminal halogens joined by single bonds.
HALIDE_OR_ION_PROFILES = {
    ("Si", 0, 8, -4): ("silicide_ionic", set()),
    ("Sn", 0, 8, -4): ("stannide_ionic", set()),
    ("Si", 2, 2, 0): ("singlet_dihalosilylene", {"F", "Cl"}),
    ("Sn", 2, 2, 0): ("tin_II_difluoride", {"F"}),
    ("P", 6, 0, -1): ("hexahalophosphate_V", {"F", "Cl"}),
    ("Sn", 6, 0, -2): ("hexahalostannate_IV", {"F", "Cl"}),
}


def supplemental_profile(view, atom):
    if atom.element not in {"Si", "Sn", "P"} or atom.hydrogens or atom.radical_electrons != 0:
        return None
    sums = bond_valences(view)
    i = atom.atom_id
    b, n, q = sums[i], atom.nonbonding_electrons, atom.formal_charge
    profile = HALIDE_OR_ION_PROFILES.get((atom.element, b, n, q))
    if profile is None:
        return None
    name, allowed = profile
    atoms = atoms_by_id(view)
    neighbours = [(atoms[e.b if e.a == i else e.a], e.order) for e in view.bonds if i in (e.a, e.b)]
    if allowed and (len(neighbours) != b or len({a.element for a, _ in neighbours}) != 1 or
                    any(a.element not in allowed or order != 1 or sums[a.atom_id] != 1 or
                        a.formal_charge != 0 for a, order in neighbours)):
        return None
    return {"atom_id": i, "element": atom.element, "rule": name + "/1.0",
            "bond_order_sum": b, "N": n, "r": 0, "formal_charge": q, "shell_electrons": 2*b+n,
            "coordination_number": len(neighbours), "oxidation_state": oxidation(view, i)}


def local_profile(view, atom):
    el, idx = atom.element, atom.atom_id
    bond_sums = bond_valences(view)
    b = bond_sums[idx]
    n, radical = atom.nonbonding_electrons, atom.radical_electrons
    if n is None or atom.formal_charge != VALENCE[el] - n - b:
        raise FormalError(f"{el}{idx}: inconsistent local electron ownership")
    atoms = atoms_by_id(view)
    neighbours = [(atoms[e.b if e.a == idx else e.a], e.order) for e in view.bonds if idx in (e.a, e.b)]
    if el == "Xe":
        ligands = [(a.element, order) for a, order in neighbours]
        forms = (ligands.count(("F", 1)), ligands.count(("O", 2)))
        if (atom.hydrogens or sum(forms) != len(ligands) or forms not in XENON_FORMS or
                atom.formal_charge != 0 or radical != 0 or n != 8 - b or
                any(bond_sums[a.atom_id] != order or a.formal_charge != 0 for a, order in neighbours)):
            raise Unsupported(f"Xe{idx}: expected elemental Xe or a registered terminal fluoride/oxo structure")
    elif (b, n, radical) not in LOCAL_OCCUPATIONS[el]:
        raise Unsupported(f"{el}{idx}: local occupation (bond order {b}, nonbonding {n}, unpaired {radical}) is unregistered")
    elif b >= 5 and any(order != 1 and not (a.element in {"O", "S"} and order == 2 and bond_sums[a.atom_id] == 2)
                        for a, order in neighbours):
        raise Unsupported(f"{el}{idx}: hypervalent profile requires single bonds or terminal oxo/thioxo double bonds")
    state = oxidation(view, idx)
    if state is not None and state not in REFERENCES[el].states:
        raise Unsupported(f"{el}{idx}: oxidation state {state} lies outside the registered set")
    return {"atom_id": idx, "element": el, "rule": f"{el}_local/1.0", "bond_order_sum": b, "N": n,
            "r": radical, "formal_charge": atom.formal_charge, "shell_electrons": 2*b+n,
            "oxidation_state": state, "coordination_number": len(neighbours)+atom.hydrogens}


RULES = legacy.ElementRules("umr-elements/1.0", VALENCE, METALS, frozenset(legacy.MAIN) | LOCAL,
                            OXIDATION, metal_counts, LewisView, local_profile, LOCAL, supplemental_profile)


def validate_view(view, state=None):
    """Check an isolated Lewis view using the same rules as a UMR 1.0 replay."""
    from ..v09.ledger import validate_view as validate
    with legacy.use_elements(RULES):
        if not isinstance(view, LewisView):
            view = LewisView.model_validate(view)
        return validate(view, state)
