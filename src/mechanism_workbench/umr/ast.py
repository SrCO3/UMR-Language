from __future__ import annotations

from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, JsonValue, model_validator

Identifier = Annotated[str, Field(pattern=r"^[A-Za-z_][A-Za-z0-9_.:~\-]*$")]
AtomID = Annotated[int, Field(gt=0)]


class Node(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)


class Evidence(Node):
    level: Literal["observed", "inferred", "proposed_completion"] = "inferred"
    source: str = "unknown"
    locations: list[dict[str, JsonValue]] = Field(default_factory=list)
    note: str = ""


class Atom(Node):
    atom_id: AtomID
    element: str
    isotope: Annotated[int, Field(ge=0)] = 0
    formal_charge: int | None
    nonbonding_electrons: Annotated[int, Field(ge=0)] | None
    hydrogens: Annotated[int, Field(ge=0)] = 0  # Attached, unnumbered, nontransferable H.
    radical_electrons: Annotated[int, Field(ge=0)] = 0
    stereo: dict[str, JsonValue] = Field(default_factory=lambda: {"status": "unknown"})
    label: str | None = None


class Bond(Node):
    a: AtomID
    b: AtomID
    order: Annotated[int, Field(ge=1, le=3)]
    stereo: dict[str, JsonValue] = Field(default_factory=lambda: {"status": "unknown"})

    @model_validator(mode="after")
    def distinct(self):
        if self.a == self.b:
            raise ValueError("Self bonds are not allowed")
        return self


class ChemicalState(Node):
    kind: Literal["ChemicalState"] = "ChemicalState"
    state_id: Identifier
    primary_view_id: Identifier
    total_charge: int
    multiplicity: Annotated[int, Field(ge=1)] = 1
    conditions: dict[str, JsonValue] = Field(default_factory=dict)
    annotations: dict[str, JsonValue] = Field(default_factory=dict)


class LewisView(Node):
    kind: Literal["LewisView"] = "LewisView"
    view_id: Identifier
    state_id: Identifier
    atoms: list[Atom] = Field(min_length=1)
    bonds: list[Bond] = Field(default_factory=list)
    annotations: dict[str, JsonValue] = Field(default_factory=dict)

    @model_validator(mode="after")
    def unique_graph(self):
        ids = [a.atom_id for a in self.atoms]
        if len(ids) != len(set(ids)):
            raise ValueError("Duplicate atom_id")
        keys = [tuple(sorted((b.a, b.b))) for b in self.bonds]
        if len(keys) != len(set(keys)):
            raise ValueError("Duplicate bond")
        if any(b.a not in ids or b.b not in ids for b in self.bonds):
            raise ValueError("Bond references missing atom_id")
        return self


class Pair(Node):
    op: Literal["PAIR"] = "PAIR"
    a: AtomID
    b: AtomID


class Split(Node):
    op: Literal["SPLIT"] = "SPLIT"
    a: AtomID
    b: AtomID
    target: AtomID

    @model_validator(mode="after")
    def target_is_b(self):
        if self.target != self.b:
            raise ValueError("SPLIT(a,b->b) requires target == b; reverse a,b for the other direction")
        return self


class Shift(Node):
    op: Literal["SHIFT"] = "SHIFT"
    a: AtomID
    b: AtomID
    c: AtomID


Action = Annotated[Pair | Split | Shift, Field(discriminator="op")]
EventType = Literal[
    "POLAR", "PROTON_TRANSFER", "HETEROLYSIS", "CARBONYL_ADDITION", "SN2",
    "E1_IONIZATION", "E1_ELIMINATION", "E1CB_DEPROTONATION", "E1CB_ELIMINATION", "E2",
    "RADICAL", "SET", "COORDINATION", "PERICYCLIC",
]


class ReactionEvent(Node):
    kind: Literal["ReactionEvent"] = "ReactionEvent"
    event_id: Identifier
    event_type: EventType
    reactant_view_id: Identifier
    product_state_id: Identifier
    product_view_id: Identifier
    actions: list[Action] = Field(min_length=1)
    bindings: dict[str, AtomID] = Field(default_factory=dict)
    depends_on_event: Identifier | None = None
    conditions: dict[str, JsonValue] = Field(default_factory=dict)
    explanation: str = ""
    evidence: list[Evidence] = Field(default_factory=list)


class ResonanceTransform(Node):
    kind: Literal["ResonanceTransform"] = "ResonanceTransform"
    transform_id: Identifier
    state_id: Identifier
    source_view_id: Identifier
    target_view_id: Identifier
    actions: list[Action] = Field(min_length=1)
    rule: str
    rule_version: str = "umr-resonance/0.1"
    evidence: list[Evidence] = Field(default_factory=list)


class AnalysisRecord(Node):
    kind: Literal["AnalysisRecord"] = "AnalysisRecord"
    analysis_id: Identifier
    operation: Literal["RESONANCE_SCAN", "SITE_ASSESSMENT", "PKA_QUERY", "REPRESENTATION_SENSITIVITY"]
    state_id: Identifier | None = None
    view_id: Identifier | None = None
    status: Literal["complete", "not_found", "truncated", "passed", "failed", "unknown", "unsupported"] = "unknown"
    rule_version: str = "umr/0.1"
    input: dict[str, JsonValue] = Field(default_factory=dict)
    output: dict[str, JsonValue] = Field(default_factory=dict)
    uncertainties: list[str] = Field(default_factory=list)


class MechanismDocument(Node):
    schema_version: Literal["umr/0.1"] = "umr/0.1"
    initial_state_id: Identifier
    states: list[ChemicalState] = Field(min_length=1)
    views: list[LewisView] = Field(min_length=1)
    events: list[ReactionEvent] = Field(default_factory=list)
    resonance_transforms: list[ResonanceTransform] = Field(default_factory=list)
    analyses: list[AnalysisRecord] = Field(default_factory=list)
    metadata: dict[str, JsonValue] = Field(default_factory=dict)

    @model_validator(mode="after")
    def unique_ids(self):
        for sequence, key in ((self.states, "state_id"), (self.views, "view_id"),
                              (self.events, "event_id"), (self.resonance_transforms, "transform_id"),
                              (self.analyses, "analysis_id")):
            ids = [getattr(x, key) for x in sequence]
            if len(ids) != len(set(ids)):
                raise ValueError(f"Duplicate {key}")
        return self
