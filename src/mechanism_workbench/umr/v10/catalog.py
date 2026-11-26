"""Page-linked native records shipped in the source tree and wheel."""
import json
from importlib import resources

def _folder():
    return resources.files("mechanism_workbench").joinpath("resources","named_examples")

def catalog():
    return json.loads(_folder().joinpath("index.json").read_text(encoding="utf-8"))

def example_text(name):
    if name not in catalog()["instances"]:
        raise KeyError(name)
    return _folder().joinpath(name+".umr").read_text(encoding="utf-8")
