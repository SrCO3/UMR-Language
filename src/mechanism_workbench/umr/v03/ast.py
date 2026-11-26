from typing import Literal

from pydantic import Field, model_validator

from ..ast import AtomID, Evidence, Identifier, LewisView as View01, Node
from ..v02.ast import (AnalysisRecord as Analysis02, ChemicalState, MechanismDocument as Document02,
                       ReactionEvent as Event02, ResonanceTransform)


class Ligand(Node):
    ligand_id: Identifier
    atom_ids: list[AtomID] = Field(min_length=1)
    donor_atom_ids: list[AtomID] = Field(min_length=1)
    ionic_class: Literal["L", "X", "unknown"] = "unknown"
    ionic_charge: int | None = None
    identity_note: str = "Declared fragment identity; unbound entries do not assert an isolated species."
    evidence: list[Evidence] = Field(default_factory=list)


class MetalBinding(Node):
    binding_id: Identifier
    metal_atom_id: AtomID
    ligand_id: Identifier
    donor_atom_ids: list[AtomID] = Field(min_length=1)
    mode: Literal["sigma_covalent", "dative", "eta"]
    hapticity: int = Field(default=1, ge=1)
    site_id: Identifier | None = None


class SiteRelation(Node):
    first_binding_id: Identifier
    second_binding_id: Identifier
    relation: Literal["cis", "trans", "unknown"] = "unknown"
    evidence: list[Evidence] = Field(default_factory=list)


class MetalCenter(Node):
    atom_id: AtomID
    oxidation_state: int | None = None
    d_electrons: int | None = None
    valence_electrons: int | None = None
    counting_method: Literal["ionic"] = "ionic"
    spin_state: Literal["closed_shell", "unknown"] = "unknown"
    geometry: Literal["unknown", "linear", "square_planar", "tetrahedral", "octahedral"] = "unknown"
    site_ids: list[Identifier] = Field(default_factory=list)
    vacant_site_ids: list[Identifier] | None = None
    relative_positions: list[SiteRelation] = Field(default_factory=list)
    evidence: list[Evidence] = Field(default_factory=list)


class LewisView(View01):
    metal_centers: list[MetalCenter] = Field(default_factory=list)
    ligands: list[Ligand] = Field(default_factory=list)
    metal_bindings: list[MetalBinding] = Field(default_factory=list)
    electron_ownership_convention: Literal["donor_owned_dative/1"] = "donor_owned_dative/1"

    @model_validator(mode="after")
    def unique_metal_ids(self):
        for values, key in ((self.metal_centers, "atom_id"), (self.ligands, "ligand_id"),
                            (self.metal_bindings, "binding_id")):
            ids = [getattr(v, key) for v in values]
            if len(ids) != len(set(ids)):
                raise ValueError(f"Duplicate {key}")
        return self


class ReactionEvent(Event02):
    metal_participation: bool = False
    metal_template: str | None = None
    metal_rule_version: Literal["umr-metal/0.3"] = "umr-metal/0.3"
    ligand_roles: dict[str, Identifier] = Field(default_factory=dict)
    binding_roles: dict[str, Identifier] = Field(default_factory=dict)
    site_roles: dict[str, Identifier] = Field(default_factory=dict)


class AnalysisRecord(Analysis02):
    operation: Literal["RESONANCE_SCAN", "SITE_ASSESSMENT", "PKA_QUERY", "REPRESENTATION_SENSITIVITY",
                       "METAL_ELECTRON_COUNT", "METAL_REPRESENTATION_CHECK"]
    rule_version: str = "umr/0.3"


class MechanismDocument(Document02):
    schema_version: Literal["umr/0.3"] = "umr/0.3"
    views: list[LewisView] = Field(min_length=1)
    events: list[ReactionEvent] = Field(default_factory=list)
    analyses: list[AnalysisRecord] = Field(default_factory=list)
