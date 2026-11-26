"""Experimental, explicitly selected electrophilic vinyl iododestannylation."""
from collections import Counter
from copy import deepcopy
from rdkit import Chem

from ..ast import Shift
from ..chemistry import FormalError, Unsupported, atoms_by_id, order

PROFILE = "PILOT_VINYL_TIN_I2_RETENTIVE_V1"


def execute(source, event):
    from .execution import transact
    from .ledger import validate_view
    validate_view(source)
    roles = event.general_roles
    names = ("tin", "vinyl", "partner", "iodo_carbon", "iodo_tin")
    if set(roles) != set(names) or len(set(roles.values())) != 5:
        raise FormalError("Expected five distinct Sn/vinyl/I2 atom roles")
    if event.event_type != "POLAR" or not event.evidence or any((
            event.template_id, event.rule_template, event.reactive_template, event.element_template,
            event.metal_template, event.electron_transfers, event.electronic_outcome,
            event.metal_participation, event.released_fragments)):
        raise Unsupported("Expected a single evidenced polar vinyl tin profile")
    sn, c, z, i, j = (roles[key] for key in names)
    atoms = atoms_by_id(source)
    if not set(roles.values()) <= atoms.keys():
        raise FormalError("Unknown vinyl tin role atom")
    if [atoms[a].element for a in (sn,c,z,i,j)] != ["Sn","C","C","I","I"]:
        raise FormalError("Wrong vinyl tin role elements")
    if any((atoms[a].formal_charge, atoms[a].radical_electrons) != (0,0) for a in roles.values()):
        raise Unsupported("Only neutral, locally paired centres are supported")
    from .smiles import to_rdkit
    carrier = to_rdkit(source)
    Chem.SanitizeMol(carrier)
    if any(a.GetAtomMapNum() in {c,z} and a.GetIsAromatic() for a in carrier.GetAtoms()):
        raise Unsupported("Aromatic ipso substitution is a separate channel")
    def links(a):
        return {b.b if b.a == a else b.a for b in source.bonds if a in {b.a, b.b}}
    if (order(source,c,z) != 2 or order(source,sn,c) != 1 or len(links(sn)) != 4
            or any(atoms[a].element != "C" or order(source,sn,a) != 1 for a in links(sn))
            or order(source,i,j) != 1 or links(i) != {j} or links(j) != {i}
            or atoms[i].isotope != atoms[j].isotope or atoms[sn].stereo.get("status") == "specified"
            or any(atoms[a].hydrogens for a in (sn,i,j))):
        raise Unsupported("Requires a tetraorganotin vinyl group and a separate I2 molecule")
    actions = [Shift(a=sn,b=c,c=i), Shift(a=i,b=j,c=sn)]
    if event.actions is not None and Counter(a.model_dump_json() for a in event.actions) != Counter(a.model_dump_json() for a in actions):
        raise FormalError("Supplied arrows disagree with vinyl tin roles")
    original = next(b for b in source.bonds if {b.a,b.b} == {c,z})

    def geometry(target):
        if original.stereo.get("status") != "specified":
            return
        new = next(b for b in target.bonds if {b.a,b.b} == {c,z})
        new.stereo = deepcopy(original.stereo)
        refs = [i if a == sn else a for a in original.stereo["reference_atoms"]]
        new.stereo["reference_atoms"] = refs if new.a == original.a else refs[::-1]
        new.stereo.pop("cip", None)

    target = transact(source, actions, event.product_state_id, event.product_view_id, after=geometry)
    event.actions = actions
    return target, {"event_id": event.event_id, "event_type":"POLAR", "general_event":PROFILE,
                    "execution_mode":"electron_actions", "atomic_commit":True, "formal_status":"passed",
                    "compiled_actions":[a.model_dump() for a in actions], "chemical_feasibility":"unknown",
                    "geometry_basis":"retentive channel; transfer mapped vinyl geometry only if specified on input",
                    "mechanism_assumption":"grouped polar exchange; elementary timing unestablished", "human_approved":False}
