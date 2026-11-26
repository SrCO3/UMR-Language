"""Formal execution of a participant-specified elementary electron transaction."""
from ..chemistry import FormalError, Unsupported, atoms_by_id, bond_valences, electron_count, identity
from ..v06.stereo import apply_event
from ..v07.ledger import project_symbols, symbols
from ..v07.transaction import touched
from .elements import charge, current_elements
from .ledger import validate_view


PROFILE = "AUTHORED_ELECTRON_ACTIONS"


def ownership_transfer(source, event):
    transfers = event.electron_transfers
    if event.event_type != "SET" or event.actions or len(transfers) != 1:
        raise FormalError("Ownership SET needs one explicit transfer and no simultaneous bond arrows")
    tx = transfers[0]
    ids = (tx.donor_atom_id, tx.acceptor_atom_id)
    atoms = atoms_by_id(source)
    if len(set(ids)) != 2 or not set(ids) <= atoms.keys() or set(ids) & symbols(source):
        raise FormalError("SET needs two distinct explicit atoms in the material inventory")
    if not event.mechanism_assumption:
        raise FormalError("State the assumed electron donor, acceptor and electronic model")
    after = event.conditions.get("radicals_after")
    if not isinstance(after, dict) or set(after) != {"donor","acceptor"}:
        raise FormalError("Give conditions.radicals_after with donor and acceptor occupations")
    target = source.model_copy(deep=True)
    target.state_id, target.view_id = event.product_state_id, event.product_view_id
    out = atoms_by_id(target)
    valences = bond_valences(target)
    for i, role, delta in zip(ids, ("donor","acceptor"), (-1,1)):
        a = out[i]
        if a.nonbonding_electrons is None or a.nonbonding_electrons + delta < 0:
            raise FormalError("The specified atom does not own the requested electron")
        unpaired = after[role]
        if a.radical_electrons is None:
            if a.element not in current_elements().metals or unpaired is not None:
                raise Unsupported("Unresolved metal spin must remain unknown in an ownership transaction")
        elif type(unpaired) is not int or unpaired < 0 or abs(unpaired-a.radical_electrons) != 1:
            raise FormalError("Each resolved endpoint must gain or lose one unpaired electron")
        a.nonbonding_electrons += delta
        a.radical_electrons = unpaired
        a.formal_charge = charge(a, valences[i])
    for center in target.metal_centers:
        if center.atom_id in ids:
            center.oxidation_state = center.d_electrons = center.valence_electrons = None
            center.spin_state = "unknown"
    for item in target.element_declarations:
        if item.atom_id in ids:
            item.oxidation_state = item.coordination_number = None
    validate_view(target)
    return target


def execute(source, event, transact):
    if event.execution_level != "electron_actions" or event.exception_ids or any((event.rule_template,
            event.template_id, event.element_template, event.reactive_template, event.metal_template,
            event.net_template, event.foundation_template, event.electronic_outcome, event.released_fragments)):
        raise FormalError("Supply one authored electron transaction without a second profile or declared product")
    validate_view(source)
    if event.electron_transfers:
        target = ownership_transfer(source, event)
    else:
        actions = event.actions or []
        if not actions or event.event_type not in {"POLAR","RADICAL","PERICYCLIC","SET"}:
            raise FormalError("Specify the complete arrows for one polar, radical, pericyclic or single-electron event")
        role_atoms = {i for value in event.bindings.values()
                      for i in (value if isinstance(value,list) else [value])}
        if not role_atoms <= atoms_by_id(source).keys():
            raise FormalError("A role refers to an atom outside the current inventory")
        single = event.event_type in {"RADICAL","SET"}
        if any((a.op == "FISHHOOK") != single for a in actions):
            raise FormalError("Use fishhook arrows for radical or SET events and pair arrows for closed-shell events")
        if event.event_type == "SET" and (len(actions) != 1 or actions[0].destination.kind != "atom"):
            raise FormalError("An elementary SET moves one electron to one acceptor")
        # A concerted transaction can join different fragments, but its arrows must share reaction sites.
        links = [touched([a]) for a in actions]
        connected = set(links[0])
        while True:
            joined = connected.union(*(part for part in links if part & connected))
            if joined == connected:
                break
            connected = joined
        if connected != touched(actions):
            raise FormalError("Separate unrelated arrow groups into separate elementary events")
        target = transact(source, actions, event.product_state_id, event.product_view_id)
    before, after = project_symbols(source), project_symbols(target)
    if (identity(source) != identity(target) or electron_count(before) != electron_count(after)
            or sum(a.formal_charge for a in before.atoms) != sum(a.formal_charge for a in after.atoms)):
        raise FormalError("Authored transaction violates atom, charge or electron conservation")
    stereo = apply_event(source, target, event)
    validate_view(target)
    from .ast import RULE_VERSION
    return target, {"event_id":event.event_id, "event_type":event.event_type, "rule_version":RULE_VERSION,
                    "execution_mode":"electron_actions", "atomic_commit":True, "formal_status":"passed",
                    "chemical_feasibility":"unknown", "authored":True, "stereochemistry":stereo,
                    "compiled_actions":[a.model_dump(mode="json") for a in event.actions or []]}
