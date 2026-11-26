from typing import Annotated, Literal

from pydantic import Field

from ..ast import (AnalysisRecord as Analysis01, AtomID, ChemicalState as State01,
                   Identifier, LewisView, MechanismDocument as Document01,
                   Node, Pair, ReactionEvent as Event01, ResonanceTransform as Resonance01,
                   Shift, Split)


class RadicalSource(Node):
    kind: Literal["radical"] = "radical"
    atom_id: AtomID
    slot: Annotated[int, Field(ge=0)] = 0


class BondSource(Node):
    kind: Literal["bond"] = "bond"
    a: AtomID
    b: AtomID
    pair: Literal["sigma", "pi1", "pi2"]
    slot: Literal[0, 1]


class AtomDestination(Node):
    kind: Literal["atom"] = "atom"
    atom_id: AtomID


class BondDestination(Node):
    kind: Literal["bond"] = "bond"
    a: AtomID
    b: AtomID
    pair: Literal["sigma", "pi1", "pi2"]


class Fishhook(Node):
    op: Literal["FISHHOOK"] = "FISHHOOK"
    source: Annotated[RadicalSource | BondSource, Field(discriminator="kind")]
    destination: Annotated[AtomDestination | BondDestination, Field(discriminator="kind")]


Action = Annotated[Pair | Split | Shift | Fishhook, Field(discriminator="op")]


class ChemicalState(State01):
    multiplicity: Annotated[int, Field(ge=1)] | None = None
    electronic_mode: Literal["CLOSED_SHELL", "RADICAL"] = "CLOSED_SHELL"


class ReactionEvent(Event01):
    event_type: Literal["POLAR", "PROTON_TRANSFER", "HETEROLYSIS", "CARBONYL_ADDITION", "SN2",
                        "E1_IONIZATION", "E1_ELIMINATION", "E1CB_DEPROTONATION", "E1CB_ELIMINATION",
                        "E2", "RADICAL", "PERICYCLIC", "SET", "ISC", "COORDINATION"]
    actions: list[Action] | None = None
    bindings: dict[str, AtomID | list[AtomID]] = Field(default_factory=dict)
    template_id: str | None = None
    template_version: str | None = None
    region_atom_ids: list[AtomID] = Field(default_factory=list)
    activation: Literal["thermal", "photochemical", "hν", "unknown"] = "unknown"
    direction: Literal["forward", "reverse"] = "forward"
    facial_mode: str = "unknown"
    rotation_mode: str = "unknown"
    mechanism_assumption: str | None = None
    mechanism_evidence: Literal["assumed_by_template", "declared_only", "supported", "unknown"] = "unknown"
    model_assumption: str | None = None
    radical_polarity: Literal["nucleophilic", "electrophilic", "ambiphilic", "unknown"] = "unknown"
    spin_correlation: Literal["singlet", "triplet", "unknown"] = "unknown"
    product_multiplicity: Annotated[int, Field(ge=1)] | None = None
    stereo_constraints: list[dict] = Field(default_factory=list)


class ResonanceTransform(Resonance01):
    actions: list[Action] = Field(min_length=1)
    resonance_mode: Literal["POLAR", "RADICAL"] = "POLAR"
    rule_version: str = "umr-resonance/0.2"


class AnalysisRecord(Analysis01):
    status: Literal["complete", "not_found", "truncated", "passed", "failed", "unknown",
                    "unsupported", "unsupported_for_resonance"] = "unknown"
    rule_version: str = "umr/0.2"


class MechanismDocument(Document01):
    schema_version: Literal["umr/0.2"] = "umr/0.2"
    states: list[ChemicalState] = Field(min_length=1)
    views: list[LewisView] = Field(min_length=1)
    events: list[ReactionEvent] = Field(default_factory=list)
    resonance_transforms: list[ResonanceTransform] = Field(default_factory=list)
    analyses: list[AnalysisRecord] = Field(default_factory=list)
