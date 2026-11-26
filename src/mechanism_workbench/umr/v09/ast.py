from typing import Literal
from pydantic import Field, field_validator, model_validator
from ..ast import Node,AtomID,Evidence
from ..v08.ast import (MechanismDocument as Doc08,LewisView as View08,ReactionEvent as Event08,
    ChemicalState,AnalysisRecord,ResonanceTransform,SourceStereo,EvidenceException)
from ..v05.ast import Ligand as Ligand05,MetalBinding as Binding05,MetalCenter as MetalCenter05
from ..v04.ast import ElectronicState
from .local_rules import ValenceRule, IntermediateRule, ReagentDeclaration

RULE_VERSION='umr-rules/0.9-r11'
class MetalCenter(MetalCenter05):
    geometry: Literal['unknown','linear','square_planar','tetrahedral','octahedral',
                      'trigonal_bipyramidal','square_pyramidal']='unknown'

    @model_validator(mode='before')
    @classmethod
    def accept_legacy_center(cls, value):
        return value.model_dump() if isinstance(value,MetalCenter05) and not isinstance(value,cls) else value

class Ligand(Ligand05):
    ionic_class: Literal['L','X','X2','X3','L2','LX','L2X','LX2','L2X2','unknown']='unknown'
    ligand_type: Literal['generic','CO','NHC','oxo','imido','hydride','Fischer_carbene','Schrock_alkylidene','Ru_alkylidene','Schrock_alkylidyne','metallacyclobutadiene','Cp','allyl']='generic'
    internal_metal_ids:list[AtomID]=Field(default_factory=list)
    remote_ionic_charge:int=0
class MetalBinding(Binding05):
    ionic_allocation:int|None=None
    donated_electrons:int|None=Field(default=None,ge=0,le=12)
    pi_pool_slot:int|None=Field(default=None,ge=0,le=1)
class MaterialBoundary(Node):
    boundary_id:str
    domain:Literal['gas','bulk_metal','surface','external_electron_reservoir']
    atom_ids:list[AtomID]
    direction:Literal['inside','in','out']='inside'
    evidence:list[Evidence]=Field(min_length=1)
    model_note:str
class ElectronOwnershipTransfer(Node):
    donor_atom_id:AtomID
    acceptor_atom_id:AtomID
    electron_count:Literal[1]=1
class LewisView(View08):
    element_rule_version:Literal['umr-elements/0.7','umr-elements/0.8','umr-elements/0.9']='umr-elements/0.9'
    metal_centers:list[MetalCenter]=Field(default_factory=list)
    ligands:list[Ligand]=Field(default_factory=list)
    metal_bindings:list[MetalBinding]=Field(default_factory=list)
    material_boundaries:list[MaterialBoundary]=Field(default_factory=list)
class ReactionEvent(Event08):
    reactive_rule_version:Literal['umr-reactive/0.4','umr-aryne/0.9-r2']='umr-reactive/0.4'
    general_event:str|None=None
    general_roles:dict[str,int]=Field(default_factory=dict)
    released_fragments:list[list[AtomID]]=Field(default_factory=list)
    electronic_outcome:ElectronicState|None=None
    electron_transfers:list[ElectronOwnershipTransfer]=Field(default_factory=list)
    phase_id:str|None=None
class ReactionPhase(Node):
    phase_id:str
    kind:Literal['reaction','workup','quench']
    event_ids:list[str]=Field(default_factory=list)
    activate_atom_ids:list[AtomID]=Field(default_factory=list)
    conditions:dict=Field(default_factory=dict)
    evidence:list[Evidence]=Field(min_length=1)
    status:Literal['specified','missing_details']='specified'
class EvidenceException(EvidenceException):
    target_view:LewisView|None=None
class CycleRecord(Node):
    cycle_id:str
    kind:Literal['catalytic','radical_propagation']
    event_ids:list[str]=Field(min_length=1)
    activation_event_ids:list[str]=Field(default_factory=list)
    termination_event_ids:list[str]=Field(default_factory=list)
    initial_view_id:str
    final_view_id:str
    persistent_core_ids:list[AtomID]=Field(default_factory=list)
    carrier_mapping:dict[int,int]
    product_atom_groups:list[list[AtomID]]=Field(min_length=1)
    exchangeable_atom_ids:list[AtomID]=Field(default_factory=list)
    evidence:list[Evidence]=Field(min_length=1)
    closure_claim:Literal['candidate']='candidate'
    @field_validator('carrier_mapping', mode='before')
    @classmethod
    def json_atom_keys(cls, value):
        if not isinstance(value,dict):return value
        result={}
        for key,target in value.items():
            if isinstance(key,str):
                if not key.isascii() or not key.isdecimal() or str(int(key))!=key:
                    raise ValueError('Carrier keys must be canonical integer atom IDs')
                key=int(key)
            if key in result:raise ValueError('Duplicate carrier atom ID')
            result[key]=target
        return result
class MechanismDocument(Doc08):
    schema_version:Literal['umr/0.9']='umr/0.9'
    views:list[LewisView]=Field(min_length=1)
    events:list[ReactionEvent]=Field(default_factory=list)
    cycles:list[CycleRecord]=Field(default_factory=list)
    phases:list[ReactionPhase]=Field(default_factory=list)
    evidence_exceptions:list[EvidenceException]=Field(default_factory=list)
    valence_rules:list[ValenceRule]=Field(default_factory=list, exclude_if=lambda rows: not rows)
    intermediate_rules:list[IntermediateRule]=Field(default_factory=list, exclude_if=lambda rows: not rows)
    reagents:list[ReagentDeclaration]=Field(default_factory=list, exclude_if=lambda rows: not rows)
    @model_validator(mode='after')
    def new_references(self):
        eids={e.event_id for e in self.events}
        if len({c.cycle_id for c in self.cycles})!=len(self.cycles):raise ValueError('Duplicate cycle ID')
        for c in self.cycles:
            if not set(c.event_ids+c.activation_event_ids+c.termination_event_ids)<=eids:raise ValueError('Cycle references absent event')
        rules = [*self.valence_rules, *self.intermediate_rules]
        views = {v.view_id for v in self.views} | {e.product_view_id for e in self.events}
        if any(r.view_id not in views for r in rules):
            raise ValueError('Local rule references an undeclared Lewis view')
        if len({r.rule_id for r in rules}) != len(rules):
            raise ValueError('Duplicate local rule ID')
        if len({(r.view_id, r.atom_id) for r in self.valence_rules}) != len(self.valence_rules):
            raise ValueError('A view can declare only one local valence per atom')
        if len({r.reagent_id for r in self.reagents}) != len(self.reagents):
            raise ValueError('Duplicate reagent ID')
        phases = {p.phase_id: set(p.activate_atom_ids) for p in self.phases}
        for reagent in self.reagents:
            if reagent.view_id not in views:
                raise ValueError('Reagent references an undeclared Lewis view')
            if reagent.phase_id is not None and not set(reagent.atom_ids) <= phases.get(reagent.phase_id, set()):
                raise ValueError('Reagent atoms must be introduced by the declared phase')
        return self
def migrate(doc):
    if doc.schema_version!='umr/0.8':raise ValueError('Only explicit v0.8 copy migration is supported')
    d=doc.model_dump(mode='json');d['schema_version']='umr/0.9'
    d['metadata'].update(human_approved=False,migration={'from':'umr/0.8','to':'umr/0.9','reviews_inherited':False})
    return MechanismDocument.model_validate(d)
