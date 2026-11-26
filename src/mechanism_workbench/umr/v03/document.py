from ..chemistry import FormalError, Unsupported
from ..document import RunResult
from ..v02.ast import ChemicalState
from .ast import AnalysisRecord
from .execution import execute_event
from .ledger import compare_views, signature, validate_view
from .resonance import resonance_scan, validate_resonance


def run_document(document,*,max_views=32,max_expansions=128):
    result=RunResult("passed",document.model_copy(deep=True))
    doc=result.document
    states={s.state_id:s for s in doc.states};views={v.view_id:v for v in doc.views}
    analyses={a.analysis_id:a for a in doc.analyses};links={r.transform_id:r for r in doc.resonance_transforms}
    active=set();established=set();completed={};pending=list(doc.resonance_transforms)

    def scan(view,label):
        found,transforms,analysis=resonance_scan(view,state=states[view.state_id],max_views=max_views,
            max_expansions=max_expansions,analysis_id=f"auto.scan.{label}")
        analyses[analysis.analysis_id]=analysis
        for v in found:
            if v.view_id in views: compare_views(v,views[v.view_id])
            else: views[v.view_id]=v
            established.add(v.view_id)
        for r in transforms: links[r.transform_id]=r
        try:
            ledger=validate_view(view,states[view.state_id])
            a=AnalysisRecord(analysis_id=f"auto.metal.{label}",operation="METAL_ELECTRON_COUNT",state_id=view.state_id,
                view_id=view.view_id,status="complete",output={"centers":ledger["metal_counts"],
                    "chemical_feasibility":"unknown","electrons_added_to_global_ledger":0})
            analyses[a.analysis_id]=a
        except Unsupported: pass

    def resolve():
        progress=True
        while progress:
            progress=False
            for r in pending[:]:
                if r.source_view_id not in established: continue
                if r.target_view_id not in views: raise FormalError("Resonance target view is missing")
                validate_resonance(views[r.source_view_id],views[r.target_view_id],r)
                established.add(r.target_view_id);pending.remove(r);progress=True
            for v in list(views.values()):
                if v.view_id in established or v.state_id not in active: continue
                for known_id in list(established):
                    known=views[known_id]
                    if known.state_id==v.state_id and signature(known)==signature(v):
                        compare_views(known,v);established.add(v.view_id);progress=True;break
    try:
        initial=states.get(doc.initial_state_id)
        if initial is None or initial.primary_view_id not in views: raise FormalError("Missing initial state/primary view")
        scan(views[initial.primary_view_id],"initial")
        validate_view(views[initial.primary_view_id],initial);active.add(initial.state_id);resolve()
        for e in doc.events:
            resolve()
            if e.reactant_view_id not in established: raise FormalError("Event reactant view has not been established")
            if e.product_state_id in active: raise FormalError("ReactionEvent must create a new state")
            source=views[e.reactant_view_id]
            validate_view(source,states[source.state_id])
            if e.event_type in {"E1_ELIMINATION","E1CB_ELIMINATION"}:
                required="E1_IONIZATION" if e.event_type=="E1_ELIMINATION" else "E1CB_DEPROTONATION"
                parent=completed.get(e.depends_on_event)
                if not parent or parent.event_type!=required or parent.product_state_id!=source.state_id or parent.bindings!=e.bindings:
                    raise FormalError("Elimination requires its matching distinct precursor event")
            elif e.depends_on_event and e.depends_on_event not in completed: raise FormalError("Unknown/forward event dependency")
            target,report=execute_event(source,e)
            if target.view_id in views: compare_views(target,views[target.view_id])
            origin=states[source.state_id]
            multiplicity=origin.multiplicity
            if not source.metal_centers:
                from ..v02.execution import product_spin
                multiplicity,spin=product_spin(origin,e,target)
                report["spin"]=spin
            elif e.product_multiplicity is not None and e.product_multiplicity!=origin.multiplicity:
                raise Unsupported("Metal adapter does not implement spin-changing reaction events")
            product_state=states.get(target.state_id)
            if product_state:
                if product_state.primary_view_id!=target.view_id: raise FormalError("Product primary view reference mismatch")
            else:
                product_state=ChemicalState(state_id=target.state_id,primary_view_id=target.view_id,total_charge=origin.total_charge,
                    multiplicity=multiplicity,electronic_mode="RADICAL" if any(a.radical_electrons for a in target.atoms) else "CLOSED_SHELL",
                    conditions={**origin.conditions,**e.conditions})
            validate_view(target,product_state)
            views[target.view_id]=target;states[target.state_id]=product_state
            established.add(target.view_id);active.add(target.state_id);completed[e.event_id]=e;result.events.append(report)
            scan(target,f"after.{e.event_id}")
        resolve()
        if pending: raise FormalError("Unreachable resonance transform")
        if any(v.view_id not in established for v in doc.views): raise FormalError("Unestablished asserted view")
        if any(s.state_id not in active for s in doc.states): raise FormalError("Unreachable asserted state")
    except Unsupported as exc: result.status="unsupported";result.errors.append(str(exc))
    except FormalError as exc: result.status="failed";result.errors.append(str(exc))
    doc.states=list(states.values());doc.views=list(views.values());doc.analyses=list(analyses.values());doc.resonance_transforms=list(links.values())
    return result
