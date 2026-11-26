"""Reuse v0.7 orchestration, approvals and BFS scans; separate net path accounting."""
from ..v07.document import run_document as run07, RunResult07
from .ast import MechanismDocument, RULE_VERSION
from .execution import execute_event


class RunResult08(RunResult07):
    def to_dict(self):
        data=super().to_dict(); data['rule_version']=RULE_VERSION
        data.update(net_transformation_steps=sum(e.get('execution_mode')=='net_transformation' for e in self.events),
                    full_electron_path_verified=self.status=='passed',
                    net_path_closed=self.execution_complete)
        return data


def run_document(document,*,approval_store=None,max_views=16,max_expansions=64):
    from . import resonance as backend
    result=run07(document,approval_store=approval_store,max_views=max_views,max_expansions=max_expansions,
                 _document_class=MechanismDocument,_executor=execute_event,_backend=backend)
    dependencies={}
    events={e.event_id:e for e in result.document.events}
    views={v.view_id:v for v in result.document.views}
    for report in result.events:
        e=events[report['event_id']]; source=views[e.reactant_view_id].state_id
        inherited=dependencies.get(source,[])
        dependencies[e.product_state_id]=inherited+([e.event_id] if report.get('execution_mode')=='net_transformation' else [])
        report['depends_on_net_events']=dependencies[e.product_state_id]
        if dependencies[e.product_state_id]:report['whole_path_strictly_verified']=False
    has_net=any(r.get('execution_mode')=='net_transformation' for r in result.events)
    if has_net:
        result.strict_validation['full_electron_path_status']='not_verified_contains_net_event'
        if result.strict_validation['status']=='passed':result.strict_validation['status']='net_ledger_passed'
        if result.status=='passed':result.status='completed_with_net_events'
    result.__class__=RunResult08
    return result
