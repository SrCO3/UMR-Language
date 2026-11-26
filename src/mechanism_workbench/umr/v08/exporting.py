"""v0.8 export retains source assertions and never draws net ledger recipes as arrows."""
from pathlib import Path
from ..v07.exporting import export_bundle as export07
from ..exporting import write_json


def export_bundle(result,directory,*,draw=True):
    manifest=export07(result,directory,draw=draw)
    manifest.update(schema_version='umr-export/0.8',full_electron_path_verified=result.status=='passed',
                    net_transformation_steps=sum(e.get('execution_mode')=='net_transformation' for e in result.events),
                    source_stereo_is_execution_proof=False)
    write_json(Path(directory)/'manifest.json',manifest)
    index=Path(directory)/'index.html'
    index.write_text(index.read_text(encoding='utf-8').replace('<title>UMR 0.7</title>','<title>UMR 0.8</title>'),encoding='utf-8')
    return manifest
