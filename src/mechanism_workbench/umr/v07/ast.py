from typing import Literal, Annotated
from pydantic import Field, model_validator
from ..ast import Node, Identifier, Evidence, AtomID
from ..ast import Pair, Split, Shift
from ..v02.ast import (Fishhook as Fishhook02, RadicalSource, BondSource, AtomDestination, BondDestination)
from ..v06.ast import (ReactionEvent as Event06, MechanismDocument as Document06,
                      ChemicalState, LewisView as View06, ResonanceTransform, AnalysisRecord)

RULE_VERSION = 'umr-rules/0.7-r1'


class LewisView(View06):
    element_rule_version: Literal['umr-elements/0.7'] = 'umr-elements/0.7'


class LonePairSource(Node):
    kind: Literal['lone_pair'] = 'lone_pair'
    atom_id: AtomID
    pair_index: int = Field(default=0,ge=0)
    slot: Literal[0,1]


class LonePairDestination(Node):
    kind: Literal['lone_pair'] = 'lone_pair'
    atom_id: AtomID


class Fishhook(Fishhook02):
    source: Annotated[RadicalSource | BondSource | LonePairSource, Field(discriminator='kind')]
    destination: Annotated[AtomDestination | BondDestination | LonePairDestination, Field(discriminator='kind')]


Action=Annotated[Pair | Split | Shift | Fishhook,Field(discriminator='op')]


class EvidenceException(Node):
    schema_version: Literal['umr-exception/0.1'] = 'umr-exception/0.1'
    exception_id: Identifier
    event_id: Identifier
    mode: Literal['rule_waiver','manual_bridge','external_event']
    rule_ids: list[str] = Field(default_factory=list)
    rule_version: str = RULE_VERSION
    explanation: str = Field(min_length=1)
    evidence: list[Evidence] = Field(min_length=1)
    conditions: dict = Field(default_factory=dict)
    domain: Literal['chemistry','nuclear','other'] = 'chemistry'
    target_state: ChemicalState | None = None
    target_view: LewisView | None = None
    correspondence: list[dict] = Field(default_factory=list)

    @model_validator(mode='after')
    def complete(self):
        if self.mode=='rule_waiver':
            if not self.rule_ids or self.target_state or self.target_view:
                raise ValueError('A waiver needs rule IDs and must not inject a product state')
        elif not self.target_state or not self.target_view:
            raise ValueError('A bridge requires a complete target state and view')
        if self.target_state and (self.target_state.state_id!=self.target_view.state_id or
                                  self.target_state.primary_view_id!=self.target_view.view_id):
            raise ValueError('Bridge state/view references disagree')
        return self


class ReactionEvent(Event06):
    actions: list[Action] | None = None
    exception_ids: list[Identifier] = Field(default_factory=list)
    rule_template: str | None = None
    expansion: dict = Field(default_factory=dict)


class MechanismDocument(Document06):
    schema_version: Literal['umr/0.7'] = 'umr/0.7'
    record_id: Identifier = 'document'
    revision: int = Field(default=1,ge=1)
    events: list[ReactionEvent] = Field(default_factory=list)
    views: list[LewisView] = Field(min_length=1)
    evidence_exceptions: list[EvidenceException] = Field(default_factory=list)
    symbolic_mode: Literal['disabled','unchanged_boundary'] = 'disabled'

    @model_validator(mode='after')
    def references(self):
        for seq,key in [(self.states,'state_id'),(self.views,'view_id'),(self.events,'event_id'),
                        (self.evidence_exceptions,'exception_id'),(self.resonance_transforms,'transform_id')]:
            names=[getattr(x,key) for x in seq]
            if len(names)!=len(set(names)):raise ValueError('Duplicate '+key)
        events={e.event_id:e for e in self.events}; exceptions={x.exception_id:x for x in self.evidence_exceptions}
        for e in self.events:
            if len(e.exception_ids)!=len(set(e.exception_ids)):raise ValueError('Duplicate exception reference')
            for key in e.exception_ids:
                if key not in exceptions or exceptions[key].event_id!=e.event_id:raise ValueError('Exception/event reference mismatch')
        if any(x.event_id not in events or x.exception_id not in events[x.event_id].exception_ids for x in self.evidence_exceptions):
            raise ValueError('Unbound evidence exception')
        return self


def migrate(document):
    if document.schema_version!='umr/0.6':raise ValueError('Only explicit 0.6 to 0.7 copy migration is supported')
    data=document.model_dump(mode='json');data['schema_version']='umr/0.7'
    for view in data['views']:view['element_rule_version']='umr-elements/0.7'
    data['metadata']={**data.get('metadata',{}),'human_approved':False,
        'migration':{'from':'umr/0.6','to':'umr/0.7','requires_reexecution':True,'reviews_inherited':False}}
    return MechanismDocument.model_validate(data)
