from ..v07.document import run_document as shared,RunResult07
from ..chemistry import FormalError
from .ast import MechanismDocument,RULE_VERSION
from .execution import execute_event
from .cycles import validate_cycle
from .local_rules import active_rules, active_reagents, rules_by_view, use_rules
class RunResult09(RunResult07):
    def to_dict(self):
        d=super().to_dict()
        if self.document.reagents:
            d['reagents']=[r.model_dump(exclude_none=True) for r in self.document.reagents]
        if self.document.valence_rules or self.document.intermediate_rules:
            d['local_rules']={'valence':[r.model_dump() for r in self.document.valence_rules],
                              'intermediate':[r.model_dump() for r in self.document.intermediate_rules]}
        transitions=[r['event_id'] for r in self.events if r.get('execution_mode')=='declared_electronic_transition']
        assumptions=[r['event_id'] for r in self.events if r.get('occupation_assumption')]
        grouped=[r['event_id'] for r in self.events if r.get('mechanistic_resolution')=='grouped_redox_proton_ledger_not_asserted_elementary']
        d.update(rule_version=RULE_VERSION,cycles=self.cycle_reports,
            phases=self.phase_reports,workflow_completion=self.workflow_completion,
            inferred_quench=getattr(self,'inferred_quench',None),
            full_electron_path_verified=self.status=='passed' and self.workflow_completion=='declared_sequence_complete' and not transitions and not assumptions and not grouped,
            declared_electronic_transition_steps=len(transitions),reaction_event_steps=len(self.events)-len(transitions),
            declared_electronic_transition_ids=transitions,occupation_assumption_event_ids=assumptions,
            grouped_ledger_event_ids=grouped,
            cycle_closure_independent_of_path_status=True)
        stereo=getattr(self,'stereo_path',{})
        d.update(stereo_path=stereo,stereochemistry_fully_specified=not bool(stereo.get('uncertain_events',[])),
            stereo_uncertainty_event_ids=sorted({r['event_id'] for r in stereo.get('uncertain_events',[])}),
            unknown_intermediate_stereo_blocks_formal_execution=False)
        return d
def run_document(document,*,approval_store=None,max_views=16,max_expansions=64,propose_quench=True,resonance_cache=None):
    with use_rules(document):
        return _run_document(document,approval_store=approval_store,max_views=max_views,max_expansions=max_expansions,
                             propose_quench=propose_quench,resonance_cache=resonance_cache)

def _run_document(document,*,approval_store,max_views,max_expansions,propose_quench,resonance_cache,
                  document_class=MechanismDocument,event_executor=execute_event):
    from . import resonance
    from .phases import PhasePlan
    backend=resonance
    if resonance_cache is not None:
        from types import SimpleNamespace
        from copy import deepcopy
        from ...libraries.core.models import digest
        entries=resonance_cache.setdefault('entries',{})
        local_key = {vid:[r.model_dump() for r in rules] for vid,rules in rules_by_view(document).items()}
        def cached_scan(view,**kw):
            key=digest([document.schema_version,view.model_dump(mode='json'),local_key,
                {k:v.model_dump(mode='json') if hasattr(v,'model_dump') else sorted(v) if isinstance(v,set) else v for k,v in kw.items()}])
            if key in entries:
                resonance_cache['hits']=resonance_cache.get('hits',0)+1
                return deepcopy(entries[key])
            resonance_cache['misses']=resonance_cache.get('misses',0)+1
            answer=resonance.resonance_scan(view,**kw)
            if len(entries)>=512:
                entries.pop(next(iter(entries)))
                resonance_cache['evictions']=resonance_cache.get('evictions',0)+1
            entries[key]=deepcopy(answer)
            return answer
        backend=SimpleNamespace(**{k:getattr(resonance,k) for k in
            ('validate_view','signature','compare_views','validate_resonance','prepare_source','product_spin')},resonance_scan=cached_scan)
    phases=PhasePlan(document)
    def staged(source,event,**kw):
        trial=phases.before(source,event);target,report=event_executor(source,trial,**kw)
        from .elements import current_elements
        target=current_elements().view_type.model_validate(target.model_dump(mode='json'))
        event.actions=trial.actions
        report.update(phase_id=event.phase_id,effective_conditions=trial.conditions)
        used=[r.rule_id for vid in (source.view_id,target.view_id) for r in active_rules(vid)]
        if used:report['local_rule_ids']=used
        reagents=[r.reagent_id for vid in (source.view_id,target.view_id) for r in active_reagents(vid)]
        if reagents:report['reagent_ids']=reagents
        return target,report
    r=shared(document,approval_store=approval_store,max_views=max_views,max_expansions=max_expansions,
             _document_class=document_class,_executor=staged,_backend=backend)
    r.__class__=RunResult09;r.cycle_reports=[]
    dependencies={};events={e.event_id:e for e in r.document.events};views={v.view_id:v for v in r.document.views}
    for report in r.events:
        e=events[report['event_id']];prior=dependencies.get(views[e.reactant_view_id].state_id,[])
        dependencies[e.product_state_id]=prior+([e.event_id] if report.get('execution_mode')=='net_transformation' else [])
        report['depends_on_net_events']=dependencies[e.product_state_id]
    if any(dependencies.values()) and r.status=='passed':r.status='completed_with_net_events'
    for cycle in r.document.cycles:
        try:out=validate_cycle(r.document,r.events,cycle)
        except (FormalError,KeyError,ValueError) as exc:out={'status':'invalid_cycle','compressible':False,'reason':str(exc),'reward':0}
        r.cycle_reports.append({'cycle_id':cycle.cycle_id,**out})
    r.phase_reports=phases.report(r.events)
    r.workflow_completion='pending_workup' if any(p['status']=='pending_workup' for p in r.phase_reports) else 'declared_sequence_complete' if r.execution_complete else 'incomplete_reaction'
    r.inferred_quench=None
    from .stereo_path import analyze_path,VERSION as STEREO_VERSION
    from .ast import AnalysisRecord
    from ...libraries.core.models import digest
    r.stereo_path=analyze_path(r.document)
    analysis_id='stereo.path.'+digest([STEREO_VERSION,r.document.record_id,r.document.revision,
        [e.model_dump(mode='json') for e in r.document.events],[v.model_dump(mode='json') for v in r.document.views]])[:20]
    if not any(a.analysis_id==analysis_id for a in r.document.analyses):
        r.document.analyses.append(AnalysisRecord(analysis_id=analysis_id,operation='STEREOCHEMISTRY',status='complete',
            rule_version=STEREO_VERSION,output=r.stereo_path,
            uncertainties=['Intermediate geometry or selectivity may remain unknown; formal ledger success is independent.'] if r.stereo_path['uncertain_events'] else []))
    if propose_quench:
        from .quench import propose
        r.inferred_quench=propose(r)
    return r
