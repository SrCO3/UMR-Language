from ..chemistry import FormalError
from ..v09 import document as legacy
from ..v09.execution import _execute_event as existing, execute_event as checked_event
from ..v09.local_rules import use_rules
from ..v09.elements import use_elements
from .ast import MechanismDocument, RULE_VERSION
from .elements import RULES
from . import metal
from .observations import evaluate


def dispatch(view, event, *, waivers=frozenset()):
    if event.general_event in metal.SPECS:
        if waivers:
            raise FormalError("Metal coordination and electron accounting require the structural checks")
        return metal.execute(view,event)
    return existing(view,event,waivers=waivers)


def execute(view, event, **kwargs):
    with use_elements(RULES):
        target, receipt = checked_event(view,event,_dispatcher=dispatch,**kwargs)
        receipt["element_rule_version"] = RULES.version
        return target, receipt


class Result(legacy.RunResult09):
    def to_dict(self):
        receipt = super().to_dict()
        receipt["rule_version"] = RULE_VERSION
        receipt["element_rule_version"] = RULES.version
        receipt["observations"] = self.observation_reports
        receipt["observed_net_steps"] = len(self.observation_reports)
        receipt["full_electron_path_verified"] = bool(self.events) and not self.document.observations and receipt["full_electron_path_verified"]
        return receipt


def run_document(document, *, approval_store=None, max_views=16, max_expansions=64,
                 propose_quench=True, resonance_cache=None):
    with use_rules(document):
        result = legacy._run_document(document,approval_store=approval_store,max_views=max_views,
            max_expansions=max_expansions,propose_quench=propose_quench and not document.observations,resonance_cache=resonance_cache,
            document_class=MechanismDocument,event_executor=execute)
        result.__class__ = Result
        result.observation_reports = []
        if result.status=="passed" and document.observations:
            views = {v.view_id:v for v in result.document.views}
            for observation in document.observations:
                try:
                    result.observation_reports.append(evaluate(views[observation.view_id],observation))
                except (FormalError, ValueError) as exc:
                    result.status = "failed"
                    result.execution_complete = False
                    result.workflow_completion = "incomplete_observation"
                    message = f"{observation.observation_id}: {exc}"
                    result.errors.append(message)
                    result.strict_validation = {**result.strict_validation,"status":"failed","errors":[message]}
                    break
            else:
                result.status = "completed_with_observed_net"
                result.workflow_completion = "observed_organic_endpoint"
        return result
