"""Export executed views and separately labelled observed organic products."""
from html import escape
from pathlib import Path

from ..exporting import write_json
from ..v09.exporting import export_bundle as export_views
from ..v09.local_rules import use_rules


def export_bundle(result, directory, *, draw=True):
    directory=Path(directory)
    with use_rules(result.input_document):
        manifest=export_views(result,directory,draw=draw)
    report=result.to_dict()
    manifest.update(schema_version="umr-export/1.0", rule_version=report["rule_version"],
                    observed_net_steps=report["observed_net_steps"],
                    full_electron_path_verified=report["full_electron_path_verified"])
    rows=report["observations"]
    write_json(directory/"observations.json",rows)
    html=(directory/"index.html").read_text("utf-8").replace("UMR 0.9","UMR 1.0")
    for number,row in enumerate(rows,1):
        folder=directory/"observations"
        folder.mkdir(exist_ok=True)
        stem=f"organic-net-{number}"
        (folder/(stem+".smi")).write_text(row["mapped_product"]+"\n",encoding="utf-8")
        html+=f'<section><h2>Observed organic product · {escape(row["id"])}</h2>'
        html+=f'<p>Recorded net transformation. The atom inventory of the unresolved remainder is in <a href="observations.json">observations.json</a>.</p>'
        if draw:
            from rdkit import Chem
            from rdkit.Chem.Draw import rdMolDraw2D
            mol=Chem.MolFromSmiles(row["mapped_product"])
            canvas=rdMolDraw2D.MolDraw2DSVG(760,340)
            canvas.DrawMolecule(mol);canvas.FinishDrawing()
            (folder/(stem+".svg")).write_text(canvas.GetDrawingText(),encoding="utf-8")
            html+=f'<img src="observations/{stem}.svg" alt="Mapped observed organic product">'
        html+="</section>"
    (directory/"index.html").write_text(html,encoding="utf-8")
    write_json(directory/"manifest.json",manifest)
    return manifest
