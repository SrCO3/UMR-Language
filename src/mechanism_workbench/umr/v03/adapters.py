"""Derived RDKit/SMILES adapter. A metal AST sidecar is always mandatory."""
from pathlib import Path
from xml.etree import ElementTree as ET

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D

from ..ast import Atom, Bond, LewisView as View01
from ..chemistry import VALENCE, FormalError, Unsupported
from ..smiles import to_rdkit as base_rdkit
from ..v02.ast import ChemicalState
from .ast import Ligand, LewisView, MetalBinding, MetalCenter
from .ledger import GROUP, METALS


def to_rdkit(view):
    plain=View01(**{k:v for k,v in view.model_dump().items() if k in View01.model_fields})
    rw=Chem.RWMol(base_rdkit(plain))
    ids={a.GetAtomMapNum():a.GetIdx() for a in rw.GetAtoms()}
    for binding in view.metal_bindings:
        if binding.mode=="dative":
            rw.AddBond(ids[binding.donor_atom_ids[0]],ids[binding.metal_atom_id],Chem.BondType.DATIVE)
    mol=rw.GetMol();Chem.SanitizeMol(mol)
    return mol


def export_mapped_smiles(view):
    return Chem.MolToSmiles(to_rdkit(view),canonical=True,isomericSmiles=True)


def import_mapped_smiles(text,state_id="s0",view_id="v0"):
    params=Chem.SmilesParserParams();params.removeHs=False
    mol=Chem.MolFromSmiles(text,params)
    if mol is None: raise FormalError("Invalid mapped SMILES")
    maps=[a.GetAtomMapNum() for a in mol.GetAtoms()]
    if min(maps,default=0)<=0 or len(maps)!=len(set(maps)): raise FormalError("All explicit atoms require stable unique positive maps")
    Chem.Kekulize(mol,clearAromaticFlags=True)
    bonds=[];datives=[]
    for b in mol.GetBonds():
        a,c=maps[b.GetBeginAtomIdx()],maps[b.GetEndAtomIdx()]
        if b.GetBondType()==Chem.BondType.DATIVE:
            datives.append((a,c))
        elif b.GetBondType() in {Chem.BondType.SINGLE,Chem.BondType.DOUBLE,Chem.BondType.TRIPLE}:
            stereo={"status":"unknown"}
            if b.GetStereo()!=Chem.BondStereo.STEREONONE:
                stereo={"status":"specified","rdkit_stereo":str(b.GetStereo()),"reference_atoms":[maps[i] for i in b.GetStereoAtoms()]}
            bonds.append(Bond(a=a,b=c,order=int(b.GetBondTypeAsDouble()),stereo=stereo))
        else: raise Unsupported("Unsupported mapped bond type; retain original input")
    atoms=[];mids=set()
    for a in mol.GetAtoms():
        i=a.GetAtomMapNum();element=a.GetSymbol();h=a.GetTotalNumHs(includeNeighbors=False)
        v=h+sum(b.order for b in bonds if i in {b.a,b.b})
        group=GROUP.get(element,VALENCE.get(element))
        n=None if group is None else group-a.GetFormalCharge()-v
        stereo={"status":"unknown"}
        if str(a.GetChiralTag())!="CHI_UNSPECIFIED":
            stereo={"status":"specified","rdkit_tag":str(a.GetChiralTag()),"neighbor_order":[x.GetAtomMapNum() for x in a.GetNeighbors()]}
        atoms.append(Atom(atom_id=i,element=element,isotope=a.GetIsotope(),formal_charge=a.GetFormalCharge() if element!="*" else None,
            nonbonding_electrons=n,hydrogens=h,radical_electrons=a.GetNumRadicalElectrons() if element not in METALS else 0,
            label="R" if element=="*" else None,stereo=stereo))
        if element in METALS-{"Li","Na","K"}: mids.add(i)
    if any(m not in mids or d in mids for d,m in datives): raise Unsupported("Dative direction must be organic donor -> metal")
    # Preserve ligand fragments without inferring their ionic class or oxidation state.
    organic=set(maps)-mids;components=[]
    while organic:
        pending=[organic.pop()];seen=set(pending)
        while pending:
            i=pending.pop()
            for b in bonds:
                if i in (b.a,b.b):
                    j=b.b if b.a==i else b.a
                    if j in organic: organic.remove(j);seen.add(j);pending.append(j)
        components.append(seen)
    ligands=[];bindings=[]
    for members in components:
        contacts=[(d,m,"dative") for d,m in datives if d in members]
        contacts += [(b.a if b.a in members else b.b,b.b if b.a in members else b.a,"sigma_covalent")
                     for b in bonds if bool({b.a,b.b}&mids) and bool({b.a,b.b}&members)]
        if not contacts: continue
        lid=f"ligand_{min(members)}";donors=sorted({d for d,_,_ in contacts})
        ligands.append(Ligand(ligand_id=lid,atom_ids=sorted(members),donor_atom_ids=donors,
            identity_note="Imported connectivity only: ligand ionic class and charge remain unknown; classify before execution."))
        for index,(d,m,mode) in enumerate(contacts):
            bindings.append(MetalBinding(binding_id=f"bind_{lid}_{index}",metal_atom_id=m,ligand_id=lid,donor_atom_ids=[d],mode=mode))
    view=LewisView(state_id=state_id,view_id=view_id,atoms=atoms,bonds=bonds,
        metal_centers=[MetalCenter(atom_id=i) for i in sorted(mids)],ligands=ligands,metal_bindings=bindings,
        annotations={"raw_mapped_smiles":text,"metal_assignment_status":"incomplete_unknown" if mids else "not_applicable",
                     "eta_information":"Not inferable from bare SMILES; authoritative sidecar required."})
    state=ChemicalState(state_id=state_id,primary_view_id=view_id,total_charge=sum(a.GetFormalCharge() for a in mol.GetAtoms()),
        multiplicity=None if any(a.radical_electrons for a in atoms) else 1,
        electronic_mode="RADICAL" if any(a.radical_electrons for a in atoms) else "CLOSED_SHELL")
    return state,view


def draw_view(view,output:Path,show_ids=True,region_atom_ids=()):
    output.parent.mkdir(parents=True,exist_ok=True)
    mol=to_rdkit(view);ids={a.GetAtomMapNum():a.GetIdx() for a in mol.GetAtoms()}
    for a in mol.GetAtoms():
        number=a.GetAtomMapNum();a.SetAtomMapNum(0)
        if show_ids: a.SetProp("atomNote",str(number))
    rdDepictor.Compute2DCoords(mol,canonOrient=True,clearConfs=True)
    eta=[b for b in view.metal_bindings if b.mode=="eta"]
    legend=f"{view.state_id} / {view.view_id}"
    if eta: legend+="; dashed = eta2 pi binding (one 2e unit)"
    paths=[]
    for cls,suffix in ((rdMolDraw2D.MolDraw2DSVG,".svg"),(rdMolDraw2D.MolDraw2DCairo,".png")):
        drawing=cls(420,250);drawing.DrawMolecule(mol,legend=legend)
        drawing.SetColour((0.1,0.35,0.6))
        for binding in eta:
            start=drawing.GetDrawCoords(ids[binding.metal_atom_id])
            p,q=[drawing.GetDrawCoords(ids[i]) for i in binding.donor_atom_ids]
            end=Point2D((p.x+q.x)/2,(p.y+q.y)/2)
            for t in range(0,12,2):
                u,v=t/12,(t+1)/12
                drawing.DrawLine(Point2D(start.x+(end.x-start.x)*u,start.y+(end.y-start.y)*u),
                    Point2D(start.x+(end.x-start.x)*v,start.y+(end.y-start.y)*v),rawCoords=True)
        drawing.FinishDrawing();path=Path(str(output)+suffix)
        data=drawing.GetDrawingText()
        if suffix==".svg":
            root=ET.fromstring(data)
            description=ET.Element("{http://www.w3.org/2000/svg}desc")
            description.text=legend+"; coordination: "+str([b.model_dump() for b in view.metal_bindings])
            root.insert(0,description)
            data=ET.tostring(root,encoding="unicode")
        path.write_text(data,encoding="utf-8") if suffix==".svg" else path.write_bytes(data)
        paths.append(path)
    return tuple(paths)
