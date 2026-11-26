from typing import Literal
from pydantic import Field, model_validator
from ..v07.ast import (ReactionEvent as Event07, MechanismDocument as Document07, LewisView as View07,
                      ChemicalState, EvidenceException as Exception07, ResonanceTransform, Fishhook)
from ..ast import Node, Evidence, AtomID
from ...analysis.models import AnalysisRecord

RULE_VERSION = 'umr-rules/0.8-r1'


class LewisView(View07):
    element_rule_version: Literal['umr-elements/0.7','umr-elements/0.8'] = 'umr-elements/0.8'


class EvidenceException(Exception07):
    rule_version: str = RULE_VERSION
    target_view: LewisView | None = None


class SourceStereo(Node):
    atom_ids: list[AtomID] = Field(min_length=1)
    assertion: dict
    evidence: list[Evidence] = Field(min_length=1)
    basis: Literal['source_record'] = 'source_record'
    execution_derived: Literal[False] = False


class ReactionEvent(Event07):
    execution_level: Literal['electron_actions','net_transformation'] = 'electron_actions'
    net_template: str | None = None
    net_roles: dict[str,int] = Field(default_factory=dict)
    source_stereo: list[SourceStereo] = Field(default_factory=list)

    @model_validator(mode='after')
    def level(self):
        if self.execution_level == 'net_transformation':
            if not self.net_template or self.actions is not None or not self.evidence or not self.conditions:
                raise ValueError('Net event requires a registered template, evidence and conditions; no asserted microscopic arrows')
        elif self.net_template:
            raise ValueError('net_template needs net_transformation level')
        return self


class MechanismDocument(Document07):
    schema_version: Literal['umr/0.8'] = 'umr/0.8'
    events: list[ReactionEvent] = Field(default_factory=list)
    views: list[LewisView] = Field(min_length=1)
    analyses: list[AnalysisRecord] = Field(default_factory=list)
    evidence_exceptions: list[EvidenceException] = Field(default_factory=list)


def migrate(document):
    if document.schema_version != 'umr/0.7':
        raise ValueError('Only explicit 0.7 to 0.8 copy migration is supported')
    data = document.model_dump(mode='json'); data['schema_version'] = 'umr/0.8'
    data['metadata'] = {**data.get('metadata',{}), 'human_approved':False,
        'migration':{'from':'umr/0.7','to':'umr/0.8','reviews_inherited':False,'requires_reexecution':True}}
    # Do not silently change template IDs or an author's source stereochemistry.
    return MechanismDocument.model_validate(data)
