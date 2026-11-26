"""Execute a UMR document, display event outcomes and optionally draw its mechanism."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"src"))
from mechanism_workbench.umr.dsl import parse_dsl
from mechanism_workbench.umr.v09.ast import MechanismDocument as Document09
from mechanism_workbench.umr.v09.document import run_document as run09
from mechanism_workbench.umr.v10.ast import MechanismDocument as Document10
from mechanism_workbench.umr.v10.document import run_document as run10

VERSIONS = {"umr/0.9": (Document09,run09), "umr/1.0": (Document10,run10)}

def execute(path):
    source = path.read_text(encoding="utf-8-sig")
    if path.suffix == ".json":
        raw = json.loads(source)
        if raw.get("schema_version") not in VERSIONS:
            raise ValueError("Expected a UMR 0.9 or 1.0 JSON document")
        doc = VERSIONS[raw["schema_version"]][0].model_validate(raw)
    else:
        doc = parse_dsl(source)
    engine = VERSIONS[doc.schema_version][1]
    result = engine(doc, propose_quench=False, max_views=1, max_expansions=1)
    if not result.execution_complete and (doc.resonance_transforms or len(doc.views)>1):
        result = engine(doc, propose_quench=False, max_views=16, max_expansions=64)
    return result

def summary(result):
    return dict(status=result.status,events=len(result.events),
                strict=bool(result.to_dict().get("full_electron_path_verified")))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="A native .umr or UMR JSON document")
    parser.add_argument("--draw", type=Path, metavar="DIRECTORY", help="Render the executed path")
    parser.add_argument("--receipt", type=Path, metavar="JSON", help="Save the full execution receipt")
    args = parser.parse_args()
    result = execute(args.input)
    print(json.dumps({**summary(result),"errors":result.errors},ensure_ascii=False))
    for event in result.events:
        print(event["event_id"],event["formal_status"],event.get("execution_mode",""))
    if args.receipt:
        args.receipt.write_text(json.dumps(result.to_dict(),indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    if args.draw:
        from mechanism_workbench.drawing.models import RenderConfig
        from mechanism_workbench.umr.v07.exporting import render_result
        render_result(result,args.draw,RenderConfig(show_atom_ids=True))
    sys.exit(0 if result.execution_complete else 1)
