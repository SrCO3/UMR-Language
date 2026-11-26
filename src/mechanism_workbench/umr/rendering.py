from pathlib import Path
from xml.etree import ElementTree as ET

from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from .ast import LewisView
from .smiles import to_rdkit


def draw_view(view: LewisView, output: Path, show_ids: bool = True, region_atom_ids=()) -> tuple[Path, Path]:
    output.parent.mkdir(parents=True, exist_ok=True)
    mol = to_rdkit(view)
    for atom in mol.GetAtoms():
        number = atom.GetAtomMapNum()
        atom.SetAtomMapNum(0)
        if show_ids:
            atom.SetProp("atomNote", str(number))
        original = next(a for a in view.atoms if a.atom_id == number)
        if original.element == "*" and original.label:
            atom.SetProp("atomLabel", original.label)
    rdDepictor.Compute2DCoords(mol, canonOrient=True, clearConfs=True)
    svg = rdMolDraw2D.MolDraw2DSVG(420, 250)
    highlights=[a.GetIdx() for a in mol.GetAtoms() if view.atoms[a.GetIdx()].atom_id in region_atom_ids]
    svg.DrawMolecule(mol, legend=f"{view.state_id} / {view.view_id}",highlightAtoms=highlights)
    svg.FinishDrawing()
    svg_path, png_path = Path(str(output)+".svg"), Path(str(output)+".png")
    svg_path.write_text(svg.GetDrawingText(), encoding="utf-8")
    png = rdMolDraw2D.MolDraw2DCairo(420, 250)
    png.DrawMolecule(mol, legend=f"{view.state_id} / {view.view_id}",highlightAtoms=highlights)
    png.FinishDrawing()
    png_path.write_bytes(png.GetDrawingText())
    return svg_path, png_path


def draw_overview(document, directory: Path) -> Path:
    """Solid single reaction arrows; blue double-headed resonance arrows."""
    namespace = "http://www.w3.org/2000/svg"
    ET.register_namespace("", namespace)
    def element(tag, **attrs): return ET.Element(f"{{{namespace}}}{tag}", {k.replace("_", "-"): str(v) for k, v in attrs.items()})
    rows = [(e.reactant_view_id, e.product_view_id, "reaction", e.event_type) for e in document.events]
    # One explicit/discovered derivation per source/target pair, up to 8 in the overview.
    rows += [(r.source_view_id, r.target_view_id, "resonance", r.rule) for r in document.resonance_transforms[:8]]
    if not rows:
        rows = [(document.states[0].primary_view_id, None, "state", "initial state")]
    root = element("svg", width=1050, height=60+len(rows)*290, viewBox=f"0 0 1050 {60+len(rows)*290}")
    root.append(element("rect", width="100%", height="100%", fill="white"))
    caption = element("text", x=20, y=25, font_size=14, fill="#333")
    caption.text = "Reaction: solid arrow; resonance: blue two-headed arrow. Finite view overview; all structures exported separately."
    root.append(caption)
    defs = element("defs")
    for name, color in (("reaction-tip", "#222"), ("resonance-tip", "#2364a0")):
        marker = element("marker", id=name, markerWidth=8, markerHeight=8, refX=6, refY=3, orient="auto-start-reverse")
        marker.append(element("path", d="M0,0 L6,3 L0,6 Z", fill=color))
        defs.append(marker)
    root.append(defs)
    for index, (left, right, kind, title) in enumerate(rows):
        y = 45 + index*290
        for view_id, x in ((left, 0), (right, 620)):
            if view_id is None:
                continue
            path = directory / "structures" / f"{view_id}.svg"
            if not path.exists():
                continue
            child = ET.fromstring(path.read_text(encoding="utf-8"))
            child.set("x", str(x)); child.set("y", str(y))
            root.append(child)
        if right:
            line = element("line", x1=435, x2=605, y1=y+120, y2=y+120,
                           stroke="#2364a0" if kind == "resonance" else "#222", stroke_width=2)
            line.set("marker-end", f"url(#{'resonance-tip' if kind == 'resonance' else 'reaction-tip'})")
            if kind == "resonance":
                line.set("marker-start", "url(#resonance-tip)")
            root.append(line)
            text = element("text", x=520, y=y+95, text_anchor="middle", font_size=12)
            text.text = title
            root.append(text)
    output = directory / "overview.svg"
    ET.ElementTree(root).write(output, encoding="utf-8", xml_declaration=True)
    return output
