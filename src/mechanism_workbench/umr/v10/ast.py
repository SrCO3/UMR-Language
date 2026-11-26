from typing import Literal

from pydantic import Field, model_validator

from ..ast import AtomID, Evidence, Node
from ..v09.ast import LewisView as View09, MechanismDocument as Document09

RULE_VERSION = "umr-rules/1.0"


class LewisView(View09):
    element_rule_version: Literal["umr-elements/0.7", "umr-elements/0.8", "umr-elements/0.9",
                                  "umr-elements/1.0"] = "umr-elements/1.0"

    @model_validator(mode="before")
    @classmethod
    def accept_legacy_view(cls, value):
        if isinstance(value, View09) and not isinstance(value, cls):
            return value.model_dump(exclude_unset=True)
        return value

    @model_validator(mode="after")
    def read_previous_arsenic(self):
        if self.metal_centers:
            centers = {c.atom_id for c in self.metal_centers}
            arsenic = {a.atom_id for a in self.atoms if a.element == "As" and a.atom_id in centers}
            if arsenic:
                from .arsenic import localize
                localize(self, arsenic)
        return self


class MethylenationObservation(Node):
    observation_id: str
    view_id: str
    profile: Literal["NYSTED_CARBONYL_METHYLENATION"]
    carbonyl: AtomID
    carbonyl_oxygen: AtomID
    methylene: AtomID
    zinc_ends: list[AtomID] = Field(min_length=2,max_length=2)
    titanium: AtomID
    unassigned_atom_ids: list[AtomID] = Field(min_length=1)
    evidence: list[Evidence] = Field(min_length=1)


class MechanismDocument(Document09):
    schema_version: Literal["umr/1.0"] = "umr/1.0"
    views: list[LewisView] = Field(min_length=1)
    observations: list[MethylenationObservation] = Field(default_factory=list)

    @model_validator(mode="after")
    def check_observations(self):
        ids = [o.observation_id for o in self.observations]
        if len(ids) != len(set(ids)) or set(ids) & {e.event_id for e in self.events}:
            raise ValueError("Observation IDs must be unique and distinct from event IDs")
        available = {v.view_id for v in self.views} | {e.product_view_id for e in self.events}
        if any(o.view_id not in available for o in self.observations):
            raise ValueError("Observation refers to an unknown source view")
        return self
