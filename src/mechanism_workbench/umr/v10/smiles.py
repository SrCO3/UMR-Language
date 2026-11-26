"""Mapped SMILES import with the UMR 1.0 electron-ownership basis."""
from ..v07.elements import oxidation
from ..v09.elements import use_elements
from ..v09.smiles import import_mapped_smiles as import_carrier
from .elements import LOCAL, RULES


def import_mapped_smiles(smiles, *, metadata=None, **kwargs):
    with use_elements(RULES):
        state, view = import_carrier(smiles, metadata=metadata, **kwargs)
    if metadata is None:
        local_ids = {atom.atom_id for atom in view.atoms if atom.element in LOCAL}
        for declaration in view.element_declarations:
            if declaration.atom_id in local_ids:
                declaration.oxidation_state = oxidation(view, declaration.atom_id)
    return state, view
