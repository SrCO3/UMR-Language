import json

from ..dsl import dump_dsl
from ..exporting import write_json
from ..rendering import draw_overview
from ..smiles import RDKIT_VERSION
from .adapters import draw_view, export_mapped_smiles


def export_bundle(result,directory,*,draw=True):
    directory.mkdir(parents=True,exist_ok=False)
    doc=result.document
    write_json(directory/"result.json",result.to_dict())
    write_json(directory/"mechanism.json",doc.model_dump(mode="json"))
    (directory/"mechanism.jsonl").write_text(json.dumps(doc.model_dump(mode="json"),ensure_ascii=False)+"\n",encoding="utf-8")
    (directory/"mechanism.umr").write_text(dump_dsl(doc),encoding="utf-8")
    write_json(directory/"actions.json",[{"event":e.model_dump(mode="json"),
        "execution":next((r for r in result.events if r["event_id"]==e.event_id),None)} for e in doc.events])
    errors=[];derivatives=[];unavailable=[]
    for view in doc.views:
        # Write authority even if a derivative cannot be constructed.
        write_json(directory/"smiles"/f"{view.view_id}.metadata.json",{"view":view.model_dump(mode="json"),
            "authority":"../mechanism.json","lossless_without_sidecar":False,
            "eta_encoding":"Retained in sidecar; no fabricated ordinary or duplicate dative bonds."})
        try:
            smiles=export_mapped_smiles(view)
            (directory/"smiles"/f"{view.view_id}.smi").write_text(smiles+"\n",encoding="utf-8")
            if draw: draw_view(view,directory/"structures"/view.view_id)
            derivatives.append(view.view_id)
        except Exception as exc:
            errors.append({"view_id":view.view_id,"status":"unsupported","reason":str(exc)})
    for event in doc.events:
        if not any(r["event_id"]==event.event_id for r in result.events): continue
        report={"event_id":event.event_id,"status":"unsupported",
                "reason":"Metal coordination and ionic-count semantics require UMR; no lossless original MechSMILES encoding is asserted.",
                "umr_retained":True}
        unavailable.append(report)
        write_json(directory/"mechsmiles"/f"{event.event_id}.unsupported.json",report)
        write_json(directory/"mechsmiles"/f"{event.event_id}.metadata.json",{
            "schema_version":doc.schema_version,"event":event.model_dump(mode="json"),
            "authority":"../mechanism.json","execution":next(r for r in result.events if r["event_id"]==event.event_id)})
    if draw: draw_overview(doc,directory)
    manifest={"schema_version":"umr-export/0.3","status":result.status,"rdkit_version":RDKIT_VERSION,
        "structures_exported":len(derivatives),"derivative_errors":errors,"mechsmiles_unavailable":unavailable,
        "diagram_scope":"RDKit structures; directed dative bonds; dashed eta2 centroid links; separate reaction/resonance arrows. No curved electron arrow layout.",
        "chemical_feasibility":"unknown","human_approved":False,"real_model_called":False}
    write_json(directory/"manifest.json",manifest)
    return manifest
