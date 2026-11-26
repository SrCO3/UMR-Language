"""Reuse the existing event orchestrator, with explicit human trust boundaries."""
from types import SimpleNamespace
from ..chemistry import FormalError,Unsupported,atoms_by_id
from ..v02.document import run_document as shared_run,RunResult02
from ..v05.execution import product_spin as product_spin05,make_state as make_state05
from .ast import RULE_VERSION,MechanismDocument
from .ledger import validate_view,signature,compare_views,symbols
from .resonance import resonance_scan,validate_resonance
from .execution import execute_event,prepare_source,WAIVABLE


class RunResult07(RunResult02):
    def to_dict(self):
        data=super().to_dict()
        data.update(strict_validation=self.strict_validation,exception_decisions=self.exception_decisions,
            exception_usage=self.exception_usage,execution_complete=self.execution_complete,
            strict_formal_pass=self.status=='passed',human_approved=False,
            manual_bridge_steps=sum(e.get('execution_mode') in {'manual_bridge','external_event'} for e in self.events),
            electron_executed_steps=sum(e.get('execution_mode','electron_actions')=='electron_actions' for e in self.events),
            rule_version=RULE_VERSION)
        return data


def _run(document,decisions,max_views,max_expansions,*,_executor=execute_event,_backend=None):
    validate_view=globals()['validate_view'] if _backend is None else _backend.validate_view
    resonance_scan=globals()['resonance_scan'] if _backend is None else _backend.resonance_scan
    validate_resonance=globals()['validate_resonance'] if _backend is None else _backend.validate_resonance
    signature=globals()['signature'] if _backend is None else _backend.signature
    compare_views=globals()['compare_views'] if _backend is None else _backend.compare_views
    prepare=prepare_source if _backend is None else getattr(_backend,'prepare_source',prepare_source)
    doc=document.model_copy(deep=True);specs={x.exception_id:x for x in doc.evidence_exceptions}
    dependencies={};usage=[];bridges={}
    if doc.symbolic_mode=='unchanged_boundary':
        for v in doc.views:v.annotations['symbolic_mode']='unchanged_boundary'
    def execute(source,event):
        original=None;target=report=None
        def attempt(waivers=frozenset()):
            trial=event.model_copy(deep=True)
            target,report=_executor(source,trial,waivers=waivers)
            from .ast import Fishhook
            if trial.actions is not None:
                trial.actions=[Fishhook.model_validate(a.model_dump(mode='json')) if a.op=='FISHHOOK' else a for a in trial.actions]
            for key in type(event).model_fields:setattr(event,key,getattr(trial,key))
            return target,report
        try:target,report=attempt()
        except (FormalError,ValueError,KeyError) as exc:original=exc
        approved=[specs[i] for i in event.exception_ids if decisions.get(i,{}).get('status')=='approved']
        own=[]
        if original:
            strict={'status':'unsupported' if isinstance(original,Unsupported) else 'failed','reason':str(original),
                    **getattr(original,'report',{})}
            for spec in approved:
                if spec.mode=='rule_waiver':
                    if set(spec.rule_ids)-WAIVABLE:raise FormalError('Exception requests a non-waivable data/ledger rule')
                    if strict.get('rule_id') not in spec.rule_ids:continue
                    target,report=attempt(frozenset(spec.rule_ids))
                else:
                    target=spec.target_view.model_copy(deep=True);state=spec.target_state.model_copy(deep=True)
                    if (target.state_id,target.view_id)!=(event.product_state_id,event.product_view_id):raise FormalError('Bridge target differs from event target references')
                    if doc.symbolic_mode=='unchanged_boundary':target.annotations['symbolic_mode']='unchanged_boundary'
                    validate_view(target,state)
                    if spec.mode=='external_event' and spec.domain=='nuclear':
                        before=atoms_by_id(source);after=atoms_by_id(target)
                        if any((before[i].element,before[i].isotope)!=(after[i].element,after[i].isotope) for i in before.keys()&after.keys()):
                            raise FormalError('Nuclear transmutation needs new atom IDs with explicit correspondence')
                        if not spec.correspondence:raise FormalError('External nuclear event needs an explicit before/after correspondence')
                        for pair in spec.correspondence:
                            if pair.get('before_atom_id') not in before or pair.get('after_atom_id') not in after:
                                raise FormalError('External correspondence references absent species')
                    bridges[event.event_id]=state
                    report={'event_id':event.event_id,'event_type':event.event_type,'execution_mode':spec.mode,
                        'domain':spec.domain,'formal_status':'not_executed','compiled_actions':[],
                        'target_state_validation':validate_view(target,state),'atomic_commit':False,
                        'conservation_verified_for_transition':False,'nuclear_computation_performed':False,
                        'correspondence':spec.correspondence}
                own=[spec.exception_id];usage.append({'event_id':event.event_id,'exception_id':spec.exception_id,
                    'mode':spec.mode,'strict_check':strict,'decision':decisions[spec.exception_id]['decision']})
                report['original_strict_check']=strict
                break
            if target is None:raise original
        inherited=dependencies.get(source.state_id,[])
        dependencies[event.product_state_id]=sorted(set(inherited+own))
        report['depends_on_exceptions']=dependencies[event.product_state_id]
        report['whole_path_strictly_verified']=not dependencies[event.product_state_id]
        return target,report
    def spin(origin,event,target):
        if event.event_id in bridges:return bridges[event.event_id].multiplicity,{'status':'human_declared_not_computed','spin_transition_executed':False}
        return getattr(_backend,'product_spin',product_spin05)(origin,event,target)
    def state(origin,event,target,multiplicity,report):
        return bridges[event.event_id] if event.event_id in bridges else make_state05(origin,event,target,multiplicity,report)
    backend=SimpleNamespace(validate_view=validate_view,execute_event=execute,product_spin=spin,make_state=state,
        signature=signature,compare_views=compare_views,resonance_scan=resonance_scan,validate_resonance=validate_resonance,
        prepare_source=prepare)
    result=shared_run(doc,max_views=max_views,max_expansions=max_expansions,_backend=backend)
    return result,usage


def run_document(document,*,approval_store=None,max_views=16,max_expansions=64,
                 _document_class=MechanismDocument,_executor=execute_event,_backend=None):
    document=_document_class.model_validate(document.model_dump(mode='json'))
    strict,_=_run(document,{},max_views,max_expansions,_executor=_executor,_backend=_backend)
    decisions=approval_store.decisions(document) if approval_store else {x.exception_id:{'status':'not_approved','decision':None} for x in document.evidence_exceptions}
    result,usage=_run(document,decisions,max_views,max_expansions,_executor=_executor,_backend=_backend) if any(v['status']=='approved' for v in decisions.values()) else (strict,[])
    complete=result.status=='passed'
    status='completed_with_exceptions' if complete and usage else 'passed_symbolic' if complete and any(symbols(v) for v in result.document.views) else result.status
    output=RunResult07(status,result.document,result.events,result.errors,result.failed_event_report)
    strict_status='passed_symbolic' if strict.status=='passed' and any(symbols(v) for v in strict.document.views) else strict.status
    output.strict_validation={'status':strict_status,'errors':strict.errors,'failed_event_report':strict.failed_event_report,
        'event_reports':strict.events,'source_of_result':'independent run without approvals'}
    output.exception_decisions=decisions;output.exception_usage=usage;output.execution_complete=complete
    output.input_document=document.model_copy(deep=True)
    return output
