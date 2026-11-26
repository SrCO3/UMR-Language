"""Graphical keys for the abbreviations actually used in a figure."""
from functools import lru_cache
import math
import re
from xml.etree import ElementTree as ET

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D

from .geometry import samples, segment_distance, unit
from .svg import element, label, text_box


def chelate_links(root,coords,folds,obstacles):
    """An outer curve denotes the folded ligand backbone between donor atoms."""
    curves=[]
    for fold in folds:
        paths=fold.get("layout_backbone_paths",[])
        if not paths:
            continue
        m=coords[fold["attachment_atom_ids"][0]]
        for n,path in enumerate(paths):
            a,b=coords[path[0]],coords[path[-1]]
            middle=((a[0]+b[0])/2,(a[1]+b[1])/2)
            u=unit((b[0]-a[0],b[1]-a[1]))
            perp=(-u[1],u[0])
            if (middle[0]-m[0])*perp[0]+(middle[1]-m[1])*perp[1]<0:
                perp=(-perp[0],-perp[1])
            external=[o for o in obstacles if not set(o["atom_ids"]) & {path[0],path[-1]}]
            def make(bulge):
                control=(middle[0]+perp[0]*bulge,middle[1]+perp[1]*bulge)
                ends=[]
                for x,y in (a,b):
                    w=unit((control[0]-x,control[1]-y))
                    ends.append((x+w[0]*11,y+w[1]*11))
                curve=(ends[0],control,control,ends[1])
                points=samples(curve,20)
                x,y=points[len(points)//2]
                tx,ty=x+perp[0]*16,y+perp[1]*16+4
                hits=sum(o["box"][0]-2<px<o["box"][2]+2 and o["box"][1]-2<py<o["box"][3]+2
                         for px,py in points[2:-2] for o in external if o["kind"]=="text")
                hits+=sum(segment_distance(p,start,end)<2 for p in points[3:-3:2]
                          for o in external for start,end in o.get("segments",[]))
                clearance=min(math.dist(p,m) for p in points)
                return (hits,max(0,23-clearance),abs(bulge)),curve,points,(tx,ty)
            _,curve,points,(tx,ty)=min((make(v) for v in (26,39,56,-34,-52)),key=lambda r:r[0])
            p,c,d,q=curve
            cmd=f"M {p[0]:.2f} {p[1]:.2f} C {c[0]:.2f} {c[1]:.2f} {d[0]:.2f} {d[1]:.2f} {q[0]:.2f} {q[1]:.2f}"
            root.append(element("path",d=cmd,fill="none",stroke="#476255",stroke_width=1.8,
                                data_kind="chelate-backbone",data_ligand=fold["label"],
                                data_donor_ids=f"{path[0]} {path[-1]}",data_backbone_ids=" ".join(map(str,path))))
            box=[min(p[0] for p in points)-1,min(p[1] for p in points)-1,
                 max(p[0] for p in points)+1,max(p[1] for p in points)+1]
            key=f"chelate.{path[0]}.{path[-1]}"
            obstacles.append(dict(id=key,kind="bond",segments=list(zip(points,points[1:])),width=2,
                                  atom_ids=[path[0],path[-1]],box=box))
            if n==0:
                label(root,fold["label"],tx,ty,size=13,text_anchor="middle",fill="#365948",data_kind="chelate-label")
                width=text_box(fold["label"],0,0,13)[2]
                caption=text_box(fold["label"],tx-width/2,ty,13)
                obstacles.append(dict(id=key+".name",kind="text",atom_ids=[],box=caption))
            curves.append(dict(label=fold["label"],backbone_atom_ids=path,box=box))
    return curves


KEY_W,KEY_H=260,146


@lru_cache(maxsize=256)
def legend_molecule(smiles,ports,ligand):
    """Add external stubs on a display copy so implicit H respects the cut."""
    base=Chem.MolFromSmiles(smiles)
    work=Chem.RWMol(base)
    external=work.AddAtom(Chem.Atom(0))
    for i in ports:
        work.AddBond(i,external,Chem.BondType.DATIVE if ligand else Chem.BondType.SINGLE)
    mol=work.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(mol)
    Chem.Kekulize(mol,clearAromaticFlags=True)
    rdDepictor.Compute2DCoords(mol,useRingTemplates=True)
    if ligand:
        conf=mol.GetConformer()
        tip=conf.GetAtomPosition(external)
        donors=[conf.GetAtomPosition(i) for i in ports]
        u=unit((tip.x-sum(p.x for p in donors)/len(donors),tip.y-sum(p.y for p in donors)/len(donors)))
        conf.SetAtomPosition(external,Point3D(tip.x+1.15*u[0],tip.y+1.15*u[1],0))
    drawer=rdMolDraw2D.MolDraw2DSVG(KEY_W,KEY_H,-1,-1,True)
    opts=drawer.drawOptions()
    opts.padding=.13
    opts.fixedFontSize=13
    opts.clearBackground=False
    opts.includeRadicals=False
    opts.prepareMolsBeforeDrawing=False
    opts.atomLabels[external]="M" if ligand else " "
    if ligand:
        for i in ports:
            if mol.GetAtomWithIdx(i).GetAtomicNum()==6:
                opts.atomLabels[i]="C"
    Chem.WedgeMolBonds(mol,mol.GetConformer())
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    root=ET.fromstring(drawer.GetDrawingText())
    end=drawer.GetDrawCoords(external)
    for i in ports:
        bond=mol.GetBondBetweenAtoms(i,external).GetIdx()
        for node in root.iter():
            if f"bond-{bond}" in node.get("class","").split():
                node.set("data-key-port",str(i))
                if ligand:
                    node.set("style",re.sub(r"(stroke|fill):#[A-Fa-f0-9]{6}",r"\1:#26739c",node.get("style","")))
                    node.set("data-kind","donor-to-metal")
                    if "Z" in node.get("d","").upper():
                        numbers=[float(v) for v in re.findall(r"-?[\d.]+",node.get("d"))]
                        coords=list(zip(numbers[::2],numbers[1::2]))
                        tx,ty=min(coords,key=lambda p:math.dist(p,(end.x,end.y)))
                        start=drawer.GetDrawCoords(i)
                        ux,uy=unit((tx-start.x,ty-start.y))
                        a=(tx-5.7*ux-2.4*uy,ty-5.7*uy+2.4*ux)
                        b=(tx-5.7*ux+2.4*uy,ty-5.7*uy-2.4*ux)
                        node.set("d",f"M {a[0]} {a[1]} L {tx} {ty} L {b[0]} {b[1]} Z")
                        node.set("style","fill:#26739c;stroke:none")
                        node.set("data-role","donor-arrowhead")
        if ligand:
            continue
        start=drawer.GetDrawCoords(i)
        ux,uy=unit((end.x-start.x,end.y-start.y))
        px,py=-uy,ux
        cx,cy=end.x-ux*2,end.y-uy*2
        points=[]
        for k in range(33):
            t=k/32
            along=(t-.5)*17
            wobble=2.2*math.sin(t*4*math.pi)
            points.append((cx+px*along+ux*wobble,cy+py*along+uy*wobble))
        root.append(element("polyline",points=" ".join(f"{x:.2f},{y:.2f}" for x,y in points),
                            fill="none",stroke="#314b43",stroke_width=1.5,data_kind="wavy-cut",data_key_port=i))
    return ET.tostring(root,encoding="unicode")


def append_keys(root,folds):
    unique={}
    for f in folds:
        if "legend_ports" not in f or "ligand" not in f:
            # Older render receipts do not identify the donor/cut atom in
            # the pattern. Re-rendering their UMR supplies that correspondence.
            continue
        key=(f["label"],f["smiles"],tuple(f["legend_ports"]),f["ligand"])
        unique.setdefault(key,f)
    if not unique:
        return []
    _,_,width,top=map(float,root.get("viewBox").split())
    columns=min(len(unique),max(1,int((width-40)//(KEY_W+16))))
    row_height=KEY_H+30
    height=top+50+row_height*math.ceil(len(unique)/columns)
    group=element("g",id="abbreviation-key",data_count=len(unique))
    root.append(group)
    group.append(element("line",x1=24,y1=top+3,x2=width-24,y2=top+3,stroke="#c4d2cb",stroke_width=1))
    label(group,"缩写结构",24,top+26,size=16,fill="#25493c")
    results=[]
    for n,((name,smiles,ports,ligand),fold) in enumerate(unique.items()):
        x=(width-(columns*(KEY_W+16)-16))/2+(n%columns)*(KEY_W+16)
        y=top+39+(n//columns)*row_height
        entry=element("g",transform=f"translate({x} {y})",data_abbreviation=name,data_smiles=smiles,
                      data_port_indices=" ".join(map(str,ports)))
        group.append(entry)
        label(entry,name,KEY_W/2,15,size=15,text_anchor="middle",fill="#203e35")
        body=element("g",transform="translate(0 18)")
        for node in ET.fromstring(legend_molecule(smiles,ports,ligand)):
            body.append(node)
        entry.append(body)
        results.append(dict(label=name,smiles=smiles,ports=list(ports),ligand=ligand,source=fold["source"]))
    scale=min(float(root.get("width"))/width,math.sqrt(15_000_000/(width*height)))
    root.set("viewBox",f"0 0 {width:g} {height:g}")
    root.set("width",str(round(width*scale)))
    root.set("height",str(round(height*scale)))
    return results
