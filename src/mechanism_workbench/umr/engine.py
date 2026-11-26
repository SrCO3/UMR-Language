from __future__ import annotations

import math
from collections import Counter

from .ast import Action, Bond, LewisView, Pair, ReactionEvent, Shift, Split
from .chemistry import (FormalError, Unsupported, atoms_by_id, bond_orders, bond_valences,
                        derived_charge, electron_count, identity, order, validate_view)


def action_key(action: Action) -> tuple:
    if isinstance(action, Pair):
        return ("PAIR", action.a, action.b)
    if isinstance(action, Split):
        return ("SPLIT", action.a, action.b, action.target)
    return ("SHIFT", action.a, action.b, action.c)


def apply_actions(source: LewisView, actions: list[Action], state_id: str, view_id: str,
                  preserve_stereo: bool = False, _validator=validate_view, _charge=derived_charge) -> LewisView:
    """Read only the source; debit its resources; commit exactly once."""
    _validator(source)
    atoms = atoms_by_id(source)
    orders = bond_orders(source)
    bond_delta, lone_delta, bond_debits, lone_debits = Counter(), Counter(), Counter(), Counter()
    keys = [action_key(action) for action in actions]
    if len(keys) != len(set(keys)):
        raise FormalError("Duplicate electron consumption: repeated identical arrow")
    touched = set()

    def edge(a, b):
        return tuple(sorted((a, b)))

    for action in actions:
        ids = [action.a, action.b] + ([action.c] if isinstance(action, Shift) else [])
        if len(ids) != len(set(ids)) or any(x not in atoms for x in ids):
            raise FormalError("Action has repeated or missing atom references")
        touched.update(ids)
        if isinstance(action, Pair):
            lone_debits[action.a] += 2
            lone_delta[action.a] -= 2
            bond_delta[edge(action.a, action.b)] += 1
        else:
            source_bond = edge(action.a, action.b)
            if source_bond not in orders:
                raise FormalError(f"No source bond for {action_key(action)}")
            bond_debits[source_bond] += 1
            bond_delta[source_bond] -= 1
            if isinstance(action, Split):
                lone_delta[action.b] += 2
            else:
                bond_delta[edge(action.b, action.c)] += 1
    for atom_id, debit in lone_debits.items():
        if debit > atoms[atom_id].nonbonding_electrons - atoms[atom_id].radical_electrons:
            raise FormalError(f"Duplicate/excess lone-pair consumption at atom {atom_id}; cannot spend electrons received in this group")
    for key, debit in bond_debits.items():
        if debit > orders[key]:
            raise FormalError(f"Duplicate/excess source bond electron consumption at {key}")
    result = source.model_copy(deep=True)
    result.state_id, result.view_id = state_id, view_id
    result.bonds = []
    old_bonds = {edge(b.a, b.b): b for b in source.bonds}
    for key in sorted(set(orders) | set(bond_delta)):
        value = orders.get(key, 0) + bond_delta[key]
        if value < 0 or value > 3:
            raise FormalError(f"Invalid resulting bond order {value} at {key}")
        if value:
            previous = old_bonds.get(key)
            if previous and previous.order == value:
                bond = previous.model_copy(deep=True)
            else:
                bond = Bond(a=key[0], b=key[1], order=value)
            if not preserve_stereo and any(x in touched for x in key):
                bond.stereo = {"status": "unknown", "reason": "reaction geometry not established"}
            result.bonds.append(bond)
    for atom in result.atoms:
        atom.nonbonding_electrons += lone_delta[atom.atom_id]
        if not preserve_stereo and atom.atom_id in touched:
            atom.stereo = {"status": "unknown", "reason": "reaction stereochemistry not established"}
    valences = bond_valences(result)
    for atom in result.atoms:
        atom.formal_charge = _charge(atom, valences[atom.atom_id])
    from .stereochemistry import reassess_centers
    reassess_centers(source,result,touched,resonance=preserve_stereo)
    _validator(result)
    if identity(source) != identity(result):
        raise FormalError("Atom identity or attached H changed")
    if electron_count(source) != electron_count(result):
        raise FormalError("Electron conservation failed")
    if sum(a.formal_charge for a in source.atoms) != sum(a.formal_charge for a in result.atoms):
        raise FormalError("Total charge conservation failed")
    return result


def _dihedral(coordinates: dict, ids: list[int]) -> float | None:
    try:
        points = [[float(x) for x in coordinates[str(i)]] for i in ids]
    except (KeyError, TypeError, ValueError):
        return None
    if any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
        return None
    def sub(a, b): return [x - y for x, y in zip(a, b)]
    def dot(a, b): return sum(x * y for x, y in zip(a, b))
    def cross(a, b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    b0, b1, b2 = sub(points[0], points[1]), sub(points[2], points[1]), sub(points[3], points[2])
    length = math.sqrt(dot(b1, b1))
    if length < 1e-10:
        return None
    axis = [x / length for x in b1]
    v = sub(b0, [dot(b0, axis) * x for x in axis])
    w = sub(b2, [dot(b2, axis) * x for x in axis])
    if dot(v, v) < 1e-10 or dot(w, w) < 1e-10:
        return None
    return math.degrees(math.atan2(dot(cross(axis, v), w), dot(v, w)))


def check_event_rule(source: LewisView, event: ReactionEvent) -> dict:
    if event.event_type in {"RADICAL", "SET", "COORDINATION", "PERICYCLIC"}:
        raise Unsupported(f"{event.event_type} specialized execution is not implemented")
    atoms = atoms_by_id(source)
    for value in event.bindings.values():
        if value not in atoms:
            raise FormalError("Binding references a missing atom")
    details = {"rule_version": "umr-event-rules/0.1", "chemical_feasibility": "unknown",
               "uncertainties": ["Formal electron/graph consistency is not a kinetic or selectivity prediction."]}

    def roles(*names):
        if any(name not in event.bindings for name in names):
            raise FormalError(f"{event.event_type} requires bindings: {names}")
        ids = [event.bindings[name] for name in names]
        if len(ids) != len(set(ids)):
            raise FormalError("Distinct chemical roles must bind distinct atoms")
        return ids

    def exact(*expected):
        if Counter(action_key(x) for x in event.actions) != Counter(expected):
            raise FormalError(f"Actions do not match {event.event_type} combination rule")

    def explicit_hydrogen(h):
        if atoms[h].element != "H":
            raise FormalError("Transferred hydrogen must have its own explicit atom_id")

    def bonded(a, b, value=1):
        if order(source, a, b) != value:
            raise FormalError(f"Required bond ({a},{b}) order {value} is absent")

    if event.event_type == "HETEROLYSIS":
        if len(event.actions) != 1 or not isinstance(event.actions[0], Split):
            raise FormalError("HETEROLYSIS requires exactly one SPLIT")
    elif event.event_type == "PROTON_TRANSFER":
        base, h, acid = roles("base", "hydrogen", "acid")
        explicit_hydrogen(h)
        bonded(h, acid)
        exact(("PAIR", base, h), ("SPLIT", h, acid, acid))
    elif event.event_type == "CARBONYL_ADDITION":
        nu, carbon, oxygen = roles("nucleophile", "carbonyl_carbon", "oxygen")
        if atoms[carbon].element != "C" or atoms[oxygen].element != "O":
            raise FormalError("Carbonyl bindings must identify C and O")
        bonded(carbon, oxygen, 2)
        exact(("PAIR", nu, carbon), ("SPLIT", carbon, oxygen, oxygen))
    elif event.event_type == "SN2":
        nu, carbon, leaving = roles("nucleophile", "electrophile", "leaving_group")
        if atoms[carbon].element != "C" or bond_valences(source)[carbon] != 4:
            raise FormalError("SN2 requires a tetravalent carbon electrophile")
        bonded(carbon, leaving)
        exact(("PAIR", nu, carbon), ("SPLIT", carbon, leaving, leaving))
        details["stereochemical_result"] = "unknown: no geometry/inversion solver in v0.1"
    elif event.event_type.startswith("E1") or event.event_type == "E2":
        alpha, beta, h, leaving, base = roles("alpha_carbon", "beta_carbon", "hydrogen", "leaving_group", "base")
        if atoms[alpha].element != "C" or atoms[beta].element != "C":
            raise FormalError("Elimination alpha and beta atoms must be carbon")
        explicit_hydrogen(h)
        bonded(alpha, beta)
        if event.event_type == "E1CB_ELIMINATION":
            bonded(h, base)
            bonded(h, beta, 0)
            if atoms[beta].formal_charge != -1 or atoms[beta].nonbonding_electrons < 2:
                raise FormalError("E1cB elimination requires the bound beta carbanion Lewis view")
        else:
            bonded(h, beta)
        if event.event_type == "E1_ELIMINATION":
            bonded(alpha, leaving, 0)
            if atoms[alpha].formal_charge != 1:
                raise FormalError("E1 elimination requires its alpha carbocation intermediate")
        else:
            bonded(alpha, leaving)
        if event.event_type == "E1_IONIZATION":
            exact(("SPLIT", alpha, leaving, leaving))
        elif event.event_type == "E1_ELIMINATION":
            exact(("PAIR", base, h), ("SHIFT", h, beta, alpha))
        elif event.event_type == "E1CB_DEPROTONATION":
            exact(("PAIR", base, h), ("SPLIT", h, beta, beta))
        elif event.event_type == "E1CB_ELIMINATION":
            exact(("PAIR", beta, alpha), ("SPLIT", alpha, leaving, leaving))
        else:
            exact(("PAIR", base, h), ("SHIFT", h, beta, alpha), ("SPLIT", alpha, leaving, leaving))
            angle = _dihedral(event.conditions.get("coordinates_3d", {}), [h, beta, alpha, leaving])
            details["geometry"] = {"status": "unknown" if angle is None else
                                   "antiperiplanar_in_supplied_coordinates" if abs(abs(angle) - 180) <= 20 else
                                   "not_antiperiplanar_in_supplied_coordinates", "dihedral_degrees": angle,
                                   "is_kinetic_feasibility_proof": False}
            if angle is None:
                details["uncertainties"].append("E2 conformational requirement is unverified; 2D connectivity is insufficient.")
    return details


def execute_event(source: LewisView, event: ReactionEvent) -> tuple[LewisView, dict]:
    if event.reactant_view_id != source.view_id:
        raise FormalError("ReactionEvent reactant_view_id differs from the supplied view")
    if event.product_state_id == source.state_id:
        raise FormalError("ReactionEvent must create a new ChemicalState")
    report = check_event_rule(source, event)
    result = apply_actions(source, event.actions, event.product_state_id, event.product_view_id)
    report.update({"event_id": event.event_id, "event_type": event.event_type,
                   "atomic_commit": True, "formal_status": "passed",
                   "ledger": validate_view(result)})
    return result, report
