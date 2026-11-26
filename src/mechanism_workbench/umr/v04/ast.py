from typing import Annotated, Literal

from pydantic import Field, model_validator

from ..ast import Atom as Atom01, AtomID, Evidence, Identifier, Node
from ..v03.ast import (AnalysisRecord, ChemicalState as State03, LewisView as View03,
                       MechanismDocument as Document03, ReactionEvent as Event03, ResonanceTransform)


class Atom(Atom01):
    # Bare SMILES does not resolve the paired/unpaired occupation of low-valent centers.
    radical_electrons: Annotated[int, Field(ge=0)] | None = 0


class ElectronicState(Node):
    character: Literal["closed_shell", "open_shell", "open_shell_singlet", "unknown"] = "unknown"
    coupling: Literal["closed_shell", "high_spin", "one_unpaired", "unknown", "unimplemented"] = "unknown"
    label: str = "unknown"
    evidence: list[Evidence] = Field(default_factory=list)


class ReactiveCenter(Node):
    center_id: Identifier
    atom_id: AtomID
    species: Literal["carbene", "carbyne", "nitrene"]
    profile: Literal["carbene_paired", "carbene_two_singly", "carbyne_one_singly",
                     "carbyne_three_singly", "nitrene_paired", "nitrene_two_singly", "unknown"]
    environment: Literal["free", "metal_bound", "unknown"] = "free"
    existence_basis: Literal["formal_model", "literature_analogy", "source_supported", "unknown"] = "unknown"
    reaction_context: str = "unknown"
    evidence: list[Evidence] = Field(default_factory=list)


class AryneRegion(Node):
    region_id: Identifier
    ring_atom_ids: list[AtomID] = Field(min_length=6, max_length=6)
    aryne_atom_ids: list[AtomID] = Field(min_length=2, max_length=2)
    representation: Literal["ortho_formal_triple", "diradical", "unknown"] = "ortho_formal_triple"
    existence_basis: Literal["formal_model", "literature_analogy", "source_supported", "unknown"] = "unknown"
    reaction_context: str = "unknown"
    evidence: list[Evidence] = Field(default_factory=list)


class ChemicalState(State03):
    electronic_mode: Literal["CLOSED_SHELL", "RADICAL", "UNKNOWN"] = "UNKNOWN"
    electronic_state: ElectronicState = Field(default_factory=ElectronicState)


class LewisView(View03):
    atoms: list[Atom] = Field(min_length=1)
    electronic_state: ElectronicState = Field(default_factory=ElectronicState)
    reactive_centers: list[ReactiveCenter] = Field(default_factory=list)
    aryne_regions: list[AryneRegion] = Field(default_factory=list)

    @model_validator(mode="after")
    def local_references(self):
        ids = {a.atom_id for a in self.atoms}
        centers = [c.atom_id for c in self.reactive_centers]
        if len(centers) != len(set(centers)) or not set(centers) <= ids:
            raise ValueError("Reactive centers require unique, existing atom IDs")
        for sequence, key in ((self.reactive_centers, "center_id"), (self.aryne_regions, "region_id")):
            names = [getattr(x, key) for x in sequence]
            if len(names) != len(set(names)):
                raise ValueError(f"Duplicate {key}")
        for region in self.aryne_regions:
            if len(set(region.ring_atom_ids)) != 6 or len(set(region.aryne_atom_ids)) != 2:
                raise ValueError("Aryne ring/active atom IDs must be distinct")
            if not set(region.aryne_atom_ids) <= set(region.ring_atom_ids) <= ids:
                raise ValueError("Aryne region refers to missing atoms")
        return self


class ReactionEvent(Event03):
    # Independent of event_type=POLAR/RADICAL/etc.; not a new reaction-family enum.
    reactive_template: str | None = None
    reactive_rule_version: Literal["umr-reactive/0.4"] = "umr-reactive/0.4"


class MechanismDocument(Document03):
    schema_version: Literal["umr/0.4"] = "umr/0.4"
    states: list[ChemicalState] = Field(min_length=1)
    views: list[LewisView] = Field(min_length=1)
    events: list[ReactionEvent] = Field(default_factory=list)
