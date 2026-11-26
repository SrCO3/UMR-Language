from dataclasses import dataclass

from ..chemistry import FormalError, Unsupported, compare_views, graph_signature
from ..document import RunResult
from .ast import ChemicalState
from .execution import execute_event, product_spin
from .ledger import radical_count, validate_view
from .resonance import resonance_scan, validate_resonance


@dataclass
class RunResult02(RunResult):
    failed_event_report: dict | None = None

    def to_dict(self):
        data=super().to_dict()
        data["formal_execution"]=self.status
        data["failed_event_report"]=self.failed_event_report
        data["resonance_analysis"]=[{"analysis_id":a.analysis_id,"state_id":a.state_id,"status":a.status}
                                    for a in self.document.analyses if a.operation=="RESONANCE_SCAN"]
        return data


def run_document(document,*,max_views=32,max_expansions=128,_backend=None):
    # v0.4 reuses the established transaction/dependency/commit orchestrator.
    # Defaults bind the original v0.2 functions; historical semantics are unchanged.
    validate_view=globals()["validate_view"] if _backend is None else _backend.validate_view
    execute_event=globals()["execute_event"] if _backend is None else _backend.execute_event
    product_spin=globals()["product_spin"] if _backend is None else _backend.product_spin
    resonance_scan=globals()["resonance_scan"] if _backend is None else _backend.resonance_scan
    validate_resonance=globals()["validate_resonance"] if _backend is None else _backend.validate_resonance
    graph_signature=globals()["graph_signature"] if _backend is None else _backend.signature
    compare_views=globals()["compare_views"] if _backend is None else _backend.compare_views
    result=RunResult02("passed",document.model_copy(deep=True))
    doc=result.document
    states={s.state_id:s for s in doc.states}; views={v.view_id:v for v in doc.views}
    transforms={r.transform_id:r for r in doc.resonance_transforms}
    analyses={a.analysis_id:a for a in doc.analyses}
    established=set(); active=set(); completed={}; pending=list(doc.resonance_transforms)
    current_event=None
    def scan(view,label,seeds=None):
        found,links,analysis=resonance_scan(view,state=states[view.state_id],seeds=seeds,max_views=max_views,
                                          max_expansions=max_expansions,analysis_id=f"auto.scan.{label}")
        aliases={}
        for value in found:
            # Deduplicate newly discovered contributors against explicitly named views.
            known=next((v for v in views.values() if v.state_id==value.state_id and graph_signature(v)==graph_signature(value)),None)
            if known:
                compare_views(value,known); aliases[value.view_id]=known.view_id
                established.add(known.view_id)
            else:
                views[value.view_id]=value; established.add(value.view_id)
        for link in links:
            link.source_view_id=aliases.get(link.source_view_id,link.source_view_id)
            link.target_view_id=aliases.get(link.target_view_id,link.target_view_id)
            transforms[link.transform_id]=link
        analyses[analysis.analysis_id]=analysis
        return analysis
    def resolve():
        progress=True
        while progress:
            progress=False
            for transform in pending[:]:
                if transform.source_view_id in established:
                    if transform.target_view_id not in views:
                        raise FormalError("Missing resonance target view")
                    if transform.resonance_mode=="RADICAL" and states[transform.state_id].electronic_mode!="RADICAL":
                        raise FormalError("Radical resonance needs an explicitly declared radical electronic state")
                    validate_resonance(views[transform.source_view_id],views[transform.target_view_id],transform)
                    established.add(transform.target_view_id); pending.remove(transform); progress=True
            for value in list(views.values()):
                if value.view_id not in established and value.state_id in active:
                    known=next((views[k] for k in established if views[k].state_id==value.state_id and graph_signature(views[k])==graph_signature(value)),None)
                    if known:
                        compare_views(known,value); established.add(value.view_id); progress=True
    try:
        initial=states.get(doc.initial_state_id)
        if initial is None or initial.primary_view_id not in views:
            raise FormalError("Missing initial state/primary view")
        primary=views[initial.primary_view_id]
        scan(primary,"initial"); validate_view(primary,initial); active.add(initial.state_id); resolve()
        for event in doc.events:
            current_event=event
            resolve()
            if event.reactant_view_id not in established:
                raise FormalError("Event reactant view is not an established state representation")
            source=views[event.reactant_view_id]; origin=states[source.state_id]
            validate_view(source,origin)
            if _backend is not None and hasattr(_backend,'prepare_source'):
                prepared,link=_backend.prepare_source(source,event)
                if link is not None:
                    validate_resonance(source,prepared,link)
                    if prepared.view_id in views:compare_views(prepared,views[prepared.view_id])
                    views[prepared.view_id]=prepared;established.add(prepared.view_id)
                    transforms[link.transform_id]=link;event.reactant_view_id=prepared.view_id;source=prepared
            if event.product_state_id in active:
                raise FormalError("ReactionEvent must produce a new state_id")
            if event.event_type in {"E1_ELIMINATION","E1CB_ELIMINATION"}:
                required="E1_IONIZATION" if event.event_type=="E1_ELIMINATION" else "E1CB_DEPROTONATION"
                parent=completed.get(event.depends_on_event)
                if parent is None or parent.event_type!=required or parent.product_state_id!=source.state_id or parent.bindings!=event.bindings:
                    raise FormalError("Elimination predecessor and bindings do not match")
            elif event.depends_on_event and event.depends_on_event not in completed:
                raise FormalError("Unknown event dependency")
            target,report=execute_event(source,event)
            multiplicity,spin=product_spin(origin,event,target)
            report["spin"]=spin
            if target.view_id in views:
                compare_views(target,views[target.view_id])
            expected=states.get(target.state_id)
            if expected:
                if expected.primary_view_id!=target.view_id:
                    raise FormalError("Target state's primary view does not match event")
                if expected.multiplicity is not None and multiplicity is not None and expected.multiplicity!=multiplicity:
                    raise Unsupported("requires_spin_transition: asserted target multiplicity contradicts conserved spin")
                validate_view(target,expected)
            elif _backend is not None:
                expected=_backend.make_state(origin,event,target,multiplicity,report)
                validate_view(target,expected)
            else:
                expected=ChemicalState(state_id=target.state_id,primary_view_id=target.view_id,
                    total_charge=sum(a.formal_charge for a in target.atoms),multiplicity=multiplicity,
                    electronic_mode="RADICAL" if any(a.radical_electrons for a in target.atoms) else "CLOSED_SHELL",
                    conditions={**origin.conditions,**event.conditions},
                    annotations={**origin.annotations,"template_constraints":report.get("relative_constraints",[])})
                validate_view(target,expected)
            states[target.state_id]=expected
            views.setdefault(target.view_id,target)
            active.add(target.state_id); established.add(target.view_id); completed[event.event_id]=event
            # Commit completed event before the separately reported analysis.
            result.events.append(report)
            if event.region_atom_ids:
                center=set(event.region_atom_ids)
            else:
                center=set()
                for action in event.actions or []:
                    if action.op=="FISHHOOK":
                        for endpoint in (action.source,action.destination):
                            center.update([endpoint.atom_id] if hasattr(endpoint,"atom_id") else [endpoint.a,endpoint.b])
                    else:
                        center.update(getattr(action,key) for key in ("a","b","c") if hasattr(action,key))
            analysis=scan(target,f"after.{event.event_id}",center)
            report["resonance_analysis"]=analysis.status
        current_event=None; resolve()
        if pending or set(views)-established or set(states)-active:
            raise FormalError("Unreachable or unjustified state/view/resonance assertion")
        for analysis in doc.analyses:
            if (analysis.state_id and analysis.state_id not in active) or (analysis.view_id and analysis.view_id not in established):
                raise FormalError("Analysis references unknown state/view")
    except (FormalError,ValueError,KeyError) as exc:
        result.status="unsupported" if isinstance(exc,Unsupported) else "failed"
        result.errors=[str(exc)]
        if current_event:
            result.failed_event_report={"event_id":current_event.event_id,"formal_execution":result.status,
                "template_scope":"unsupported" if isinstance(exc,Unsupported) else "mismatch",
                "orbital_symmetry":"unknown","geometry_status":"unverified","stereochemistry":"unknown",
                "mechanism_evidence":current_event.mechanism_evidence,"kinetics_selectivity":"unknown",
                "requires_spin_transition":"requires_spin_transition" in str(exc),**getattr(exc,"report",{})}
    doc.states=list(states.values()); doc.views=list(views.values())
    doc.resonance_transforms=list(transforms.values()); doc.analyses=list(analyses.values())
    return result
