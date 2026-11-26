import json
from pathlib import Path

from .ast import MechanismDocument
from .document import RunResult
from .dsl import dump_dsl
from .pka import PKAQuery
from .rendering import draw_overview, draw_view
from .smiles import RDKIT_VERSION, export_mapped_smiles


def write_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def mechsmiles_for(event, source) -> str:
    from .chemistry import Unsupported
    if any(action.op=="FISHHOOK" for action in event.actions or []):
        raise Unsupported("Original MechSMILES cannot losslessly represent single-electron FISHHOOK actions")
    arrows = []
    for action in event.actions:
        if action.op == "PAIR":
            arrows.append(f"({action.a},{action.b})")
        elif action.op == "SPLIT":
            arrows.append(f"(({action.a},{action.b}),{action.b})")
        else:
            arrows.append(f"(({action.a},{action.b}),{action.c})")
    return export_mapped_smiles(source) + "|" + ";".join(arrows)


def export_bundle(result: RunResult, directory: Path, *, draw=True, _dispatch=True) -> dict:
    if _dispatch and result.document.schema_version=="umr/1.0":
        from .v10.exporting import export_bundle as export10
        return export10(result,directory,draw=draw)
    if _dispatch and result.document.schema_version=='umr/0.9':
        from .v09.exporting import export_bundle as export09
        return export09(result,directory,draw=draw)
    if _dispatch and result.document.schema_version=='umr/0.8':
        from .v08.exporting import export_bundle as export08
        return export08(result,directory,draw=draw)
    if _dispatch and result.document.schema_version=='umr/0.7':
        from .v07.exporting import export_bundle as export07
        return export07(result,directory,draw=draw)
    if _dispatch and result.document.schema_version in {"umr/0.5","umr/0.6"}:
        from .v04.exporting import export_bundle as export05
        return export05(result,directory,draw=draw)
    if _dispatch and result.document.schema_version=="umr/0.4":
        from .v04.exporting import export_bundle as export04
        return export04(result,directory,draw=draw)
    if result.document.schema_version=="umr/0.3":
        from .v03.exporting import export_bundle as export03
        return export03(result,directory,draw=draw)
    directory.mkdir(parents=True, exist_ok=False)
    doc = result.document
    smiles_export=export_mapped_smiles
    if doc.schema_version=='umr/0.5':
        from .v05.smiles import export_mapped_smiles as smiles_export
    if doc.schema_version=='umr/0.6':
        from .v06.smiles import export_mapped_smiles as smiles_export
    write_json(directory / "result.json", result.to_dict())
    write_json(directory / "mechanism.json", doc.model_dump(mode="json"))
    (directory / "mechanism.jsonl").write_text(json.dumps(doc.model_dump(mode="json"), ensure_ascii=False)+"\n", encoding="utf-8")
    (directory / "mechanism.umr").write_text(dump_dsl(doc), encoding="utf-8")
    views = {v.view_id: v for v in doc.views}
    derivatives, errors, mechsmiles_unavailable = [], [], []
    write_json(directory / "actions.json", [{"event_id":e.event_id,"event_type":e.event_type,
        "actions":[a.model_dump(mode="json") for a in e.actions or []]} for e in doc.events])
    for view in doc.views:
        try:
            smiles = smiles_export(view)
            folder = directory / "smiles"
            folder.mkdir(exist_ok=True)
            (folder / f"{view.view_id}.smi").write_text(smiles+"\n", encoding="utf-8")
            # Full AST accompanies all SMILES, including unknown labels, stereo and n.
            metadata=view.model_dump(mode="json")
            if doc.schema_version=="umr/0.5":
                from .v05.smiles import sidecar
                metadata=sidecar(next(s for s in doc.states if s.state_id==view.state_id),view,smiles)
            if doc.schema_version=="umr/0.6":
                from .v06.smiles import sidecar
                metadata=sidecar(next(s for s in doc.states if s.state_id==view.state_id),view,smiles)
            if doc.schema_version=="umr/0.4":
                from .v04.smiles import sidecar
                metadata=sidecar(next(s for s in doc.states if s.state_id==view.state_id),view,smiles)
            write_json(folder / f"{view.view_id}.metadata.json", metadata)
            if draw:
                region={i for e in doc.events if view.view_id in {e.reactant_view_id,e.product_view_id}
                        for i in getattr(e,"region_atom_ids",[])}
                draw_view(view, directory / "structures" / view.view_id,region_atom_ids=region)
            derivatives.append(view.view_id)
        except Exception as exc:
            errors.append({"view_id": view.view_id, "status": "unsupported", "reason": str(exc)})
    for event in doc.events:
        # Never advertise an unexecuted/unsupported event as validated MechSMILES.
        if event.event_id not in {r["event_id"] for r in result.events}:
            continue
        folder = directory / "mechsmiles"
        folder.mkdir(exist_ok=True)
        from .chemistry import Unsupported
        try:
            if doc.schema_version in {'umr/0.5','umr/0.6'} and (event.element_template or any(a.element in {'S','Se','P','Cl','Br','I','B','Na','Mg','Al','Cr','Pd','Ru','Hg'} for a in views[event.reactant_view_id].atoms)):
                raise Unsupported('Element/coordination conventions require the UMR AST and sidecar; no standalone MechSMILES claim')
            text=mechsmiles_for(event,views[event.reactant_view_id])
            (folder / f"{event.event_id}.txt").write_text(text+"\n",encoding="utf-8")
        except Unsupported as exc:
            unavailable={"event_id":event.event_id,"status":"unsupported","reason":str(exc),"umr_retained":True}
            mechsmiles_unavailable.append(unavailable)
            write_json(folder/f"{event.event_id}.unsupported.json",unavailable)
        write_json(folder / f"{event.event_id}.metadata.json", {
            "event": event.model_dump(mode="json"), "schema_version": doc.schema_version,
            "states": [s.model_dump(mode="json") for s in doc.states],
            "authority": "../mechanism.json", "lossless_without_sidecar": False,
            "note": "Conditions, event classification, evidence and grouping require this companion metadata.",
        })
    if draw:
        draw_overview(doc, directory)
    manifest = {"schema_version": "umr-export/"+doc.schema_version.split("/")[1], "status": result.status,
                "rdkit_version": RDKIT_VERSION, "structures_exported": len(derivatives),
                "derivative_errors": errors,"mechsmiles_unavailable":mechsmiles_unavailable,
                "diagram_scope":"structures, reaction/resonance arrows and action table; no curved electron-arrow layout",
                "chemical_feasibility": "unknown", "human_approved": False}
    write_json(directory / "manifest.json", manifest)
    return manifest


def write_schemas(directory: Path):
    directory.mkdir(parents=True, exist_ok=True)
    write_json(directory / "umr-v0.1.schema.json", MechanismDocument.model_json_schema())
    write_json(directory / "pka-query-v0.1.schema.json", PKAQuery.model_json_schema())
    from .v02.ast import MechanismDocument as Document02
    write_json(directory / "umr-v0.2.schema.json",Document02.model_json_schema())
    from .v03.ast import MechanismDocument as Document03
    write_json(directory / "umr-v0.3.schema.json",Document03.model_json_schema())
    from .v04.ast import MechanismDocument as Document04
    write_json(directory / "umr-v0.4.schema.json",Document04.model_json_schema())
    from .v05.ast import MechanismDocument as Document05
    write_json(directory / "umr-v0.5.schema.json",Document05.model_json_schema())
    from .v06.ast import MechanismDocument as Document06
    write_json(directory / "umr-v0.6.schema.json",Document06.model_json_schema())
    from .v07.ast import MechanismDocument as Document07, EvidenceException
    write_json(directory / 'umr-v0.7.schema.json',Document07.model_json_schema())
    write_json(directory / 'umr-evidence-exception-v0.1.schema.json',EvidenceException.model_json_schema())
    from .v08.ast import MechanismDocument as Document08
    write_json(directory / 'umr-v0.8.schema.json',Document08.model_json_schema())
    from .v09.ast import MechanismDocument as Document09
    write_json(directory / 'umr-v0.9.schema.json',Document09.model_json_schema())
    from .v10.ast import MechanismDocument as Document10
    write_json(directory / 'umr-v1.0.schema.json',Document10.model_json_schema())
