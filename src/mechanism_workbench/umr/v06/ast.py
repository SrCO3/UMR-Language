from typing import Literal
from pydantic import Field
from ..ast import AtomID, Node, Bond
from ..v05.ast import (Atom, ChemicalState, LewisView, ResonanceTransform,
                      ReactionEvent as Event05, MechanismDocument as Document05,
                      AnalysisRecord as Analysis05)


class StereoEffect(Node):
    atom_id: AtomID
    effect: Literal['retention','inversion','destroyed','created_unspecified','unknown']
    basis: str


class ReactionEvent(Event05):
    foundation_template: str | None = None
    template_roles: dict[str,int] = Field(default_factory=dict)
    stereo_effects: list[StereoEffect] = Field(default_factory=list)
    phase: Literal['reaction','acid_treatment','base_treatment','workup','unknown'] = 'unknown'


class AnalysisRecord(Analysis05):
    operation: Literal['RESONANCE_SCAN','SITE_ASSESSMENT','PKA_QUERY','REPRESENTATION_SENSITIVITY',
                       'METAL_ELECTRON_COUNT','METAL_REPRESENTATION_CHECK','ELEMENT_RULE_CHECK',
                       'ELEMENT_REPRESENTATION_CHECK','PHYSICAL_PARAMETER_QUERY',
                       'OXIDATION_STATE','STEREOCHEMISTRY']
    rule_version: str = 'umr/0.6'


class MechanismDocument(Document05):
    schema_version: Literal['umr/0.6'] = 'umr/0.6'
    events: list[ReactionEvent] = Field(default_factory=list)
    analyses: list[AnalysisRecord] = Field(default_factory=list)


def migrate(document):
    """Explicit copy migration. Never rewrite an existing artifact or its review."""
    data=document.model_dump(mode='json')
    if data['schema_version']!='umr/0.5':
        raise ValueError('Automatic copy migration is restricted to UMR 0.5')
    data['schema_version']='umr/0.6'
    data['metadata']={**data.get('metadata',{}),'migration':{'from':'umr/0.5','to':'umr/0.6',
        'reviews_inherited':False,'requires_reexecution':True},'human_approved':False}
    template=data['metadata'].get('foundation_template',{}).get('template_id')
    for e in data['events']:
        e['foundation_template']=template
        e['template_roles']=data['metadata'].get('role_bindings',{})
    return MechanismDocument.model_validate(data)
