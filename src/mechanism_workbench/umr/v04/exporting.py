from pathlib import Path

from ..exporting import export_bundle as shared_export, write_json


def export_bundle(result,directory:Path,*,draw=True):
    manifest=shared_export(result,directory,draw=False,_dispatch=False)
    if draw:
        from ...drawing.models import RenderConfig
        from ...drawing.render import draw_panel, render_document
        config=RenderConfig(max_relayouts=0,routing_seconds=1.0,max_resonance_panels=2)
        for view in result.document.views:
            try:
                state=next(s for s in result.document.states if s.state_id==view.state_id)
                draw_panel(view,None,[],view.view_id,
                    f"{view.view_id} · M={state.multiplicity if state.multiplicity is not None else 'unknown'} · {view.electronic_state.character}",
                    {},directory/"structures"/view.view_id,config,state_only=True)
            except Exception as exc:
                manifest["derivative_errors"].append({"view_id":view.view_id,"status":"unsupported","reason":str(exc)})
        try:
            rendering=render_document(result.document,directory/"drawing",config)
            manifest["drawing_status"]=rendering["status"]
            for panel in rendering["panels"]:
                if panel["status"] in {"render_failed","render_unsupported"}:
                    manifest["derivative_errors"].append({"drawing":panel["event_id"],"reason":panel.get("reason")})
        except Exception as exc:
            manifest["derivative_errors"].append({"drawing":"document","reason":str(exc)})
    manifest.update(diagram_scope="deterministic structure/electron-flow panels; state character and explicit occupation retained",
        electronic_state_authority="mechanism.json and mapped-SMILES sidecars",spin_transitions_executed=0)
    write_json(directory/"manifest.json",manifest)
    return manifest
