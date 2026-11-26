from __future__ import annotations

from dataclasses import dataclass, field

from .ast import ChemicalState, LewisView, MechanismDocument
from .chemistry import (FormalError, Unsupported, compare_views, graph_signature, validate_view)
from .engine import execute_event
from .resonance import resonance_scan, validate_resonance


@dataclass
class RunResult:
    status: str
    document: MechanismDocument
    events: list[dict] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    def to_dict(self):
        return {"status": self.status, "document": self.document.model_dump(mode="json"),
                "event_reports": self.events, "errors": self.errors,
                "real_reaction_steps": len(self.events),
                "representation_derivations": len(self.document.resonance_transforms),
                "chemical_feasibility": "unknown", "human_approved": False}


def run_document(document: MechanismDocument, *, max_views: int | None = None,
                 max_expansions: int | None = None) -> RunResult:
    # New research defaults are opt-in by schema; historical versions keep 32/128.
    if max_views is None:max_views=16 if document.schema_version in {'umr/0.6','umr/0.7','umr/0.8','umr/0.9','umr/1.0'} else 32
    if max_expansions is None:max_expansions=64 if document.schema_version in {'umr/0.6','umr/0.7','umr/0.8','umr/0.9','umr/1.0'} else 128
    if document.schema_version == 'umr/1.0':
        from .v10.document import run_document as run10
        return run10(document,max_views=max_views,max_expansions=max_expansions)
    if document.schema_version == 'umr/0.9':
        from .v09.document import run_document as run09
        return run09(document,max_views=max_views,max_expansions=max_expansions)
    if document.schema_version == 'umr/0.8':
        from .v08.document import run_document as run08
        return run08(document,max_views=max_views,max_expansions=max_expansions)
    if document.schema_version == 'umr/0.7':
        from .v07.document import run_document as run07
        return run07(document,max_views=max_views,max_expansions=max_expansions)
    if document.schema_version == "umr/0.6":
        from .v06.document import run_document as run06
        return run06(document,max_views=max_views,max_expansions=max_expansions)
    if document.schema_version == "umr/0.5":
        from .v05.document import run_document as run05
        return run05(document,max_views=max_views,max_expansions=max_expansions)
    if document.schema_version == "umr/0.4":
        from .v04.document import run_document as run04
        return run04(document,max_views=max_views,max_expansions=max_expansions)
    if document.schema_version == "umr/0.3":
        from .v03.document import run_document as run03
        return run03(document,max_views=max_views,max_expansions=max_expansions)
    if document.schema_version == "umr/0.2":
        from .v02.document import run_document as run02
        return run02(document,max_views=max_views,max_expansions=max_expansions)
    result = RunResult("passed", document.model_copy(deep=True))
    doc = result.document
    states = {s.state_id: s for s in doc.states}
    views = {v.view_id: v for v in doc.views}
    transforms = {r.transform_id: r for r in doc.resonance_transforms}
    analyses = {a.analysis_id: a for a in doc.analyses}
    established = set()
    active_states = set()
    completed = {}
    pending_transforms = list(doc.resonance_transforms)

    def scan(view, label, seeds=None):
        found, links, analysis = resonance_scan(view, seeds=seeds, max_views=max_views,
                                               max_expansions=max_expansions, analysis_id=f"auto.scan.{label}")
        for value in found:
            if value.view_id in views:
                compare_views(value, views[value.view_id])
            else:
                views[value.view_id] = value
            established.add(value.view_id)
        for link in links:
            transforms[link.transform_id] = link
        analyses[analysis.analysis_id] = analysis

    def resolve_views():
        # Explicit resonance derivations can be written before or after their target declaration.
        progress = True
        while progress:
            progress = False
            for transform in pending_transforms[:]:
                if transform.source_view_id in established:
                    if transform.target_view_id not in views:
                        raise FormalError("Missing resonance target view")
                    validate_resonance(views[transform.source_view_id], views[transform.target_view_id], transform)
                    established.add(transform.target_view_id)
                    pending_transforms.remove(transform)
                    progress = True
            # An explicitly asserted view can also match a discovered contributor.
            for value in list(views.values()):
                if value.view_id not in established and value.state_id in active_states:
                    if any(views[key].state_id == value.state_id and graph_signature(views[key]) == graph_signature(value)
                           for key in established):
                        known = next(views[key] for key in established
                                     if views[key].state_id == value.state_id and graph_signature(views[key]) == graph_signature(value))
                        compare_views(known, value)
                        established.add(value.view_id)
                        progress = True

    try:
        initial = states.get(doc.initial_state_id)
        if initial is None or initial.primary_view_id not in views:
            raise FormalError("Initial state/primary view is missing")
        primary = views[initial.primary_view_id]
        # Mandatory even when an initial graph is outside supported chemistry.
        scan(primary, "initial")
        validate_view(primary, initial)
        active_states.add(initial.state_id)
        resolve_views()
        for event in doc.events:
            resolve_views()
            if event.reactant_view_id not in established:
                raise FormalError(f"Unestablished reactant view {event.reactant_view_id}")
            source = views[event.reactant_view_id]
            validate_view(source, states[source.state_id])
            if event.product_state_id in active_states:
                raise FormalError("A ReactionEvent must produce a distinct, new state_id")
            if event.event_type in {"E1_ELIMINATION", "E1CB_ELIMINATION"}:
                required = "E1_IONIZATION" if event.event_type == "E1_ELIMINATION" else "E1CB_DEPROTONATION"
                parent = completed.get(event.depends_on_event)
                if parent is None or parent.event_type != required or parent.product_state_id != source.state_id or parent.bindings != event.bindings:
                    raise FormalError(f"{event.event_type} must depend on the matching {required} event and intermediate")
            elif event.depends_on_event and event.depends_on_event not in completed:
                raise FormalError("Unknown or forward event dependency")
            target, report = execute_event(source, event)
            if event.product_view_id in views:
                compare_views(target, views[event.product_view_id])
            else:
                views[event.product_view_id] = target
            expected_state = states.get(event.product_state_id)
            if expected_state:
                if expected_state.primary_view_id != target.view_id:
                    raise FormalError("Asserted product state's primary view differs from event output")
                validate_view(target, expected_state)
            else:
                origin = states[source.state_id]
                states[target.state_id] = ChemicalState(
                    state_id=target.state_id, primary_view_id=target.view_id,
                    total_charge=sum(a.formal_charge for a in target.atoms),
                    conditions={**origin.conditions, **event.conditions},
                    annotations=origin.annotations.copy(),
                )
            active_states.add(target.state_id)
            established.add(target.view_id)
            completed[event.event_id] = event
            result.events.append(report)
            center = {value for action in event.actions for key, value in action.model_dump().items() if key in {"a", "b", "c"}}
            scan(target, f"after.{event.event_id}", center)
        resolve_views()
        if pending_transforms:
            raise FormalError("Unreachable resonance transform source")
        if set(views) - established or set(states) - active_states:
            raise FormalError("An asserted state/view was not produced or justified by resonance")
        for analysis in doc.analyses:
            if analysis.state_id and analysis.state_id not in active_states:
                raise FormalError("Analysis references unknown state")
            if analysis.view_id and analysis.view_id not in established:
                raise FormalError("Analysis references unknown view")
    except Unsupported as exc:
        result.status, result.errors = "unsupported", [str(exc)]
    except (FormalError, ValueError, KeyError) as exc:
        result.status, result.errors = "failed", [str(exc)]
    doc.states, doc.views = list(states.values()), list(views.values())
    doc.resonance_transforms, doc.analyses = list(transforms.values()), list(analyses.values())
    return result
