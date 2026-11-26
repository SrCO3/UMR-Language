"""Read the earlier 1.0 covalent As representation at the input boundary."""
from dataclasses import replace

from ..chemistry import FormalError, Unsupported, bond_valences
from ..v07.elements import oxidation
from ..v09.ast import LewisView as PreviousView
from ..v09.coordination import checks
from ..v09.elements import use_elements
from .elements import METALS, PNICTOGEN, REFERENCES, RULES


def localize(view, ids):
    centers = [c for c in view.metal_centers if c.atom_id in ids]
    if len(centers) != len(ids):
        raise FormalError("Duplicate legacy arsenic center")
    bindings = [b for b in view.metal_bindings if b.metal_atom_id in ids]
    labels = {b.ligand_id for b in bindings}
    bound_centers = {label: {b.metal_atom_id for b in bindings if b.ligand_id == label} for label in labels}
    ligands = [l for l in view.ligands if l.ligand_id in labels]
    members = ids.union(*(l.atom_ids for l in ligands))
    if (any(b.mode != "sigma_covalent" for b in bindings)
            or any(len(attached) > 1 for attached in bound_centers.values())
            or any(c.geometry != "unknown" or c.site_ids or c.vacant_site_ids or c.relative_positions for c in centers)
            or any(b.metal_atom_id not in ids and (b.ligand_id in labels or ids.intersection(b.donor_atom_ids))
                   for b in view.metal_bindings)
            or any(ids.intersection(l.internal_metal_ids) for l in view.ligands)
            or any(a.atom_id in members and a.element in METALS for a in view.atoms)
            or any((b.a in members) != (b.b in members) for b in view.bonds)):
        raise Unsupported("Legacy As coordination requires explicit Lewis bonds or a local intermediate rule; "
                          "donor, bridging and geometry assignments must be resolved")
    # Recheck the old ligand charges, metal counts and declarations before removing
    # redundant sigma annotations. Only this private input copy uses the old partition.
    prior = PreviousView(view_id=view.view_id, state_id=view.state_id,
                         atoms=[a for a in view.atoms if a.atom_id in members],
                         bonds=[b for b in view.bonds if b.a in members],
                         metal_centers=centers, ligands=ligands, metal_bindings=bindings)
    old = replace(RULES, view_type=PreviousView, metals={"As":5}, oxidation={"As":REFERENCES["As"].states},
                  metal_counts=lambda symbol, state: (None, 5-state))
    with use_elements(old):
        checks(prior)
    sums = bond_valences(view)
    by_id = {c.atom_id:c for c in centers}
    updated = []
    for atom in view.atoms:
        if atom.atom_id not in ids:
            updated.append(atom)
            continue
        center = by_id[atom.atom_id]
        os = oxidation(view, atom.atom_id)
        if center.oxidation_state is not None and center.oxidation_state != os:
            raise FormalError(f"As{atom.atom_id}: recorded oxidation state disagrees with local bonds")
        allowed = {r for b,n,r in PNICTOGEN if (b,n) == (sums[atom.atom_id],atom.nonbonding_electrons)}
        if atom.radical_electrons is not None:
            allowed &= {atom.radical_electrons}
        if center.spin_state != "unknown":
            allowed = {r for r in allowed if (r == 0) == (center.spin_state == "closed_shell")}
        if len(allowed) != 1:
            raise Unsupported(f"As{atom.atom_id}: provide an explicit local Lewis occupation for this structure")
        updated.append(atom.model_copy(update={"radical_electrons":allowed.pop()}))
    # Retain the previous representation as provenance, including evidence and unresolved spin.
    history = {"metal_centers":[c.model_dump(exclude_defaults=True) for c in centers],
               "ligands":[l.model_dump(exclude_defaults=True) for l in ligands],
               "metal_bindings":[b.model_dump(exclude_defaults=True) for b in bindings],
               "original_unpaired":{str(a.atom_id):a.radical_electrons for a in view.atoms if a.atom_id in ids}}
    if "legacy_arsenic_coordination" in view.annotations:
        raise FormalError("Duplicate legacy arsenic coordination annotation")
    view.annotations = {**view.annotations, "legacy_arsenic_coordination":history}
    view.atoms = updated
    view.metal_centers = [c for c in view.metal_centers if c.atom_id not in ids]
    view.metal_bindings = [b for b in view.metal_bindings if b.metal_atom_id not in ids]
    view.ligands = [l for l in view.ligands if l.ligand_id not in labels]
