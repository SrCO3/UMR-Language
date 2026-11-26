from typing import Literal
from pydantic import Field
from ..umr.ast import Node


class GroupAbbreviation(Node):
    label: str = Field(min_length=1,max_length=16,pattern=r"^[\w′'()\-]+$")
    smiles: str = Field(min_length=1,max_length=500)
    source: str = Field(min_length=1,max_length=500)
    anchor_index: int = Field(default=0,ge=0)
    ligand: bool = False
    atom_ids: list[int] = Field(default_factory=list)


class RenderConfig(Node):
    version: Literal["umr-drawing/0.1","umr-drawing/0.2"] = "umr-drawing/0.2"
    bond_pixels: int = Field(default=42,ge=24,le=80)
    font_pixels: int = Field(default=17,ge=12,le=28)
    show_atom_ids: bool = False
    max_resonance_panels: int = Field(default=6,ge=0,le=32)
    curve_candidates: int = Field(default=72,ge=1,le=256)
    joint_beam: int = Field(default=8,ge=1,le=32)
    max_relayouts: int = Field(default=1,ge=0,le=3)
    max_panel_pixels: int = Field(default=4000000,ge=100000,le=12000000)
    collision_padding: float = Field(default=2.5,ge=0,le=10)
    routing_seconds: float = Field(default=4.0,gt=0,le=60)
    diagnostics: bool = False
    lone_pair_display_limit: int = Field(default=1,ge=1,le=4)
    electron_endpoint_gap: float = Field(default=5.0,ge=2,le=12)
    material_display: Literal["participants","all"] = "participants"
    use_scaffold_templates: bool = True
    abbreviate_groups: bool = True
    site_layout: bool = True
    spin_display: Literal["on_demand","details"] = "on_demand"
    abbreviations: list[GroupAbbreviation] = Field(default_factory=list)


class ReviewLocation(Node):
    event_id: str | None = None
    view_id: str | None = None
    action_id: str | None = None
    arrow_id: str | None = None
    atom_ids: list[int] = Field(default_factory=list)
    explanation: str = ""
    suggestion: str = ""


class SpotReview(Node):
    schema_version: Literal["umr-spot-review/0.1"] = "umr-spot-review/0.1"
    sample_id: str
    bundle_revision: int = Field(ge=1)
    source_sha256: str
    umr_sha256: str
    render_sha256: str
    render_version: Literal["umr-drawing/0.1","umr-drawing/0.2"] = "umr-drawing/0.2"
    reviewer: str = ""
    reviewed_at: str | None = None
    fidelity: Literal["unreviewed","approved","issue","undetermined"] = "unreviewed"
    chemistry: Literal["unreviewed","approved","issue","undetermined"] = "unreviewed"
    drawing: Literal["unreviewed","approved","issue","undetermined"] = "unreviewed"
    locations: list[ReviewLocation] = Field(default_factory=list)
    notes: str = ""
    explicit_human_submission: bool = False
    human_approved: bool = False
    submission_kind: Literal["human","synthetic_test"] = "human"
