"""Browse the two textbook catalogues and replay their linked UMR programs."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CATALOG = ROOT/"named-reactions/catalog.json"

def verify(entry):
    from umr import execute, summary
    for program in entry["programs"]:
        path = ROOT/program["path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != program["sha256"]:
            raise ValueError(f"{entry['id']} has a different input digest")
        result = summary(execute(path))
        expected = {key:program[key] for key in result}
        if result != expected:
            raise ValueError(f"{entry['id']} {path.name}: expected {expected}, obtained {result}")
        yield f"{path.name}  {result['status']}  events={result['events']}  strict={result['strict']}"

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--find",metavar="NAME",help="Look up an ID or reaction name in the printed catalogues")
    mode.add_argument("--show",metavar="ID",help="Show the catalogue entry and linked input files")
    mode.add_argument("--replay",metavar="ID",help="Check the files associated with this entry")
    mode.add_argument("--check-all",action="store_true",help="Compare all linked programs with their recorded outcomes")
    args=parser.parse_args()
    entries=json.loads(CATALOG.read_text("utf-8"))
    if args.find:
        word=args.find.casefold()
        matches=[r for r in entries if any(word in r[k].casefold() for k in ("id","name","name_zh"))]
        for r in matches:print(f"{r['id']}  p.{r['page']}  {r['name']}  {r['status']}")
        return 0 if matches else 1
    if args.check_all:
        count=sum(1 for r in entries for _ in verify(r))
        print(f"Matched the recorded outcomes for {count} programs across {len(entries)} catalogue records.")
        return 0
    selected=next((r for r in entries if r["id"].casefold()==(args.show or args.replay).casefold()),None)
    if not selected:parser.error("Unknown catalogue ID")
    if args.show:
        print(json.dumps(selected,indent=2,ensure_ascii=False))
    elif selected["programs"]:
        print("\n".join(verify(selected)))
    else:
        print(f"{selected['id']}: {selected['status']}. See its printed page and related family entries.")
        return 1
    return 0

if __name__=="__main__":sys.exit(main())
