"""SMILES is a graph carrier; spin/occupancy declarations travel in a checked sidecar."""
import hashlib

from ..chemistry import FormalError, bond_valences, identity, bond_orders
from ..smiles import import_mapped_smiles as import01, export_mapped_smiles
from .ast import ChemicalState, ElectronicState, LewisView


def sidecar(state,view,smiles):
    return {"schema_version":"umr-smiles/0.4","smiles_sha256":hashlib.sha256(smiles.encode()).hexdigest(),
        "represented_view_id":view.view_id,
        "state":state.model_dump(mode="json"),"view":view.model_dump(mode="json"),
        "lossless_without_sidecar":False,"bare_smiles_electronic_state":"ambiguous"}


def import_mapped_smiles(smiles,state_id="s0",view_id="v0",labels=None,metadata=None):
    raw_state,raw_view=import01(smiles,state_id,view_id,labels,_preserve_lewis=True)
    if metadata is not None:
        if metadata.get("schema_version")!="umr-smiles/0.4" or metadata.get("smiles_sha256")!=hashlib.sha256(smiles.encode()).hexdigest():
            raise FormalError("SMILES sidecar version/hash mismatch")
        state=ChemicalState.model_validate(metadata["state"]);view=LewisView.model_validate(metadata["view"])
        # Unknown labels and stereochemical evidence live in the sidecar; compare
        # every structural field that SMILES can carry, ignoring its spin guesses.
        def structural(v):
            return {a.atom_id:(a.element,a.isotope,a.hydrogens,a.formal_charge,a.nonbonding_electrons) for a in v.atoms},bond_orders(v)
        if structural(raw_view)!=structural(view): raise FormalError("Sidecar contradicts mapped SMILES graph/charge/H")
        # A resonance contributor need not be the state's primary view. Retain
        # both references; never relabel a contributor as primary on import.
        if view.state_id!=state.state_id or view.view_id!=metadata.get("represented_view_id") or view.electronic_state!=state.electronic_state:
            raise FormalError("Sidecar state/view identity or electronic-state mismatch")
        return state,view
    data=raw_view.model_dump(mode="json");val=bond_valences(raw_view);ambiguous=[]
    for a in data["atoms"]:
        if a["formal_charge"]==0 and ((a["element"]=="C" and (val[a["atom_id"]],a["nonbonding_electrons"]) in {(2,2),(1,3)}) or
                                      (a["element"]=="N" and (val[a["atom_id"]],a["nonbonding_electrons"])==(1,4))):
            a["radical_electrons"]=None;ambiguous.append(a["atom_id"])
    data["annotations"].update(electronic_state_source="bare SMILES; RDKit radical assignment is not a spin measurement",
        ambiguous_occupation_atom_ids=ambiguous)
    view=LewisView.model_validate(data)
    state=ChemicalState(state_id=state_id,primary_view_id=view_id,total_charge=raw_state.total_charge,
        multiplicity=None,electronic_mode="UNKNOWN",electronic_state=ElectronicState())
    return state,view
