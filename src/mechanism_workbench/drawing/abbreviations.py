"""Mapped display labels with one attachment or two recorded chelate donors."""
from functools import lru_cache

from rdkit import Chem

from .models import GroupAbbreviation

GROUPS = (
    ("dppe","P(c1ccccc1)(c1ccccc1)CCP(c1ccccc1)c1ccccc1",True),
    ("dppp","P(c1ccccc1)(c1ccccc1)CCCP(c1ccccc1)c1ccccc1",True),
    ("dppb","P(c1ccccc1)(c1ccccc1)CCCCP(c1ccccc1)c1ccccc1",True),
    ("SIMes","[C]1N(c2c(C)cc(C)cc2C)CCN1c1c(C)cc(C)cc1C",True),
    ("PPh3","P(c1ccccc1)(c1ccccc1)c1ccccc1",True),
    ("PCy3","P(C1CCCCC1)(C1CCCCC1)C1CCCCC1",True),
    ("P(tBu)3","P(C(C)(C)C)(C(C)(C)C)C(C)(C)C",True),
    ("TBDPS","[Si](c1ccccc1)(c1ccccc1)C(C)(C)C",False),
    ("TIPS","[Si](C(C)C)(C(C)C)C(C)C",False),
    ("TBS","[Si](C)(C)C(C)(C)C",False),
    ("TMS","[Si](C)(C)C",False),
    ("Boc","C(=O)OC(C)(C)C",False),
    ("Bn","Cc1ccccc1",False),("Ph","c1ccccc1",False),("tBu","C(C)(C)C",False),
)


@lru_cache(maxsize=256)
def pattern(smiles, anchor):
    q = Chem.MolFromSmiles(smiles)
    if q is None or len(Chem.GetMolFrags(q)) != 1 or not 0 <= anchor < q.GetNumAtoms() or q.GetNumAtoms() < 2:
        raise ValueError("Drawing abbreviation requires a connected group and an existing anchor")
    # An unsatisfied attachment valence in an isolated pattern (e.g. silyl,
    # NHC) must not demand a radical on the actual bound donor.
    q.GetAtomWithIdx(anchor).SetNumRadicalElectrons(0)
    return q


def find_labels(mol, view, visible, touched, definitions):
    work = Chem.Mol(mol)
    Chem.GetSymmSSSR(work)
    try:
        Chem.SetAromaticity(work)
    except (RuntimeError,ValueError):
        work = Chem.Mol(mol)
    atoms = {a.atom_id:a for a in view.atoms}
    protected = {a.atom_id for a in view.atoms if a.element in {"H","*"} or a.isotope or a.radical_electrons
                 or a.formal_charge not in (0,None) or a.label}
    for obj in (*view.atoms,*view.bonds):
        if obj.stereo.get("status") == "specified":
            protected.update([obj.atom_id] if hasattr(obj,"atom_id") else [obj.a,obj.b])
            protected.update(obj.stereo.get("reference_atoms",[]))
            protected.update(obj.stereo.get("neighbor_order",[]))
    protected.update(c.atom_id for c in getattr(view,"reactive_centers",[]))
    protected.update(i for b in getattr(view,"metal_bindings",[]) for i in (b.metal_atom_id,*b.donor_atom_ids))
    external = {a.atom_id:set() for a in view.atoms}
    for b in view.bonds:
        external[b.a].add(b.b)
        external[b.b].add(b.a)
    for b in getattr(view,"metal_bindings",[]):
        for d in b.donor_atom_ids:
            external[d].add(b.metal_atom_id)
            external[b.metal_atom_id].add(d)
    custom = [GroupAbbreviation.model_validate(d) for d in definitions]
    catalog = custom + [GroupAbbreviation(label=l,smiles=s,ligand=g,source="UMR depiction dictionary")
                        for l,s,g in GROUPS]
    candidates = []
    chelates={}
    for ligand in getattr(view,"ligands",[]):
        bound=[b for b in view.metal_bindings if b.ligand_id==ligand.ligand_id]
        donors=set(ligand.donor_atom_ids)
        if (len(donors)>=2 and bound and all(b.mode=="dative" for b in bound) and
                len({b.metal_atom_id for b in bound})==1 and donors=={i for b in bound for i in b.donor_atom_ids}):
            chelates[frozenset(ligand.atom_ids)]=donors
    for order,rule in enumerate(catalog):
        query = pattern(rule.smiles,rule.anchor_index)
        for match in work.GetSubstructMatches(query,uniquify=False,maxMatches=2048):
            ordered = [work.GetAtomWithIdx(i).GetAtomMapNum() for i in match]
            anchor = ordered[rule.anchor_index]
            members = frozenset(ordered)
            if not members <= visible or (rule.atom_ids and members != set(rule.atom_ids)):
                continue
            ports=chelates.get(members,{anchor}) if rule.ligand else {anchor}
            if anchor not in ports or (len(ports)>1 and ports<=touched):
                continue
            inside = members-ports
            if inside & (protected|touched) or (not rule.ligand and anchor in protected|touched):
                continue
            if any(atoms[i].isotope or atoms[i].radical_electrons or atoms[i].stereo.get("status")=="specified"
                   or atoms[i].label for i in ports):
                continue
            boundary = {(i,j) for i in members for j in external[i]-members}
            if any(i not in ports for i,j in boundary) or len(boundary)>len(ports) or (not boundary and not rule.ligand):
                continue
            # Exact induced topology, including unsaturation/charge on internal atoms.
            if sum(len(external[i]&members) for i in members) != 2*query.GetNumBonds():
                continue
            if any(atoms[i].formal_charge != query.GetAtomWithIdx(k).GetFormalCharge() for k,i in enumerate(ordered)):
                continue
            if rule.ligand and boundary:
                if any(not any(b.mode=="dative" and b.metal_atom_id==j and i in b.donor_atom_ids
                               for b in getattr(view,"metal_bindings",[])) for i,j in boundary):
                    continue
            if boundary:
                i,j=next(iter(boundary))
                b=next((b for b in view.bonds if {b.a,b.b}=={i,j}),None)
                if b and b.order != 1:
                    continue
            row=dict(label=rule.label,anchor_atom_id=anchor,anchor_atom_ids=sorted(ports),atom_ids=sorted(members),source=rule.source,
                     smiles=rule.smiles,attachment_atom_ids=sorted({j for _,j in boundary}))
            row["legend_ports"]=sorted(ordered.index(i) for i in ports)
            row["ligand"]=rule.ligand
            if len(ports)>1:
                # Connect the ports along the shortest actual ligand paths. A
                # spanning tree prevents extra lines between nonadjacent arms.
                paths=[]
                for i in sorted(ports):
                    for j in sorted(k for k in ports if k>i):
                        path=Chem.GetShortestPath(query,ordered.index(i),ordered.index(j))
                        paths.append([ordered[k] for k in path])
                sets={i:{i} for i in ports}
                links=[]
                for path in sorted(paths,key=lambda p:(len(p),p)):
                    first,last=path[0],path[-1]
                    if last in sets[first]:
                        continue
                    combined=sets[first]|sets[last]
                    for i in combined:
                        sets[i]=combined
                    links.append(path)
                row["layout_backbone_paths"]=links
                row["layout_backbone_atom_ids"]=list(dict.fromkeys(i for path in links for i in path))
            candidates.append((order>=len(custom),-len(members),order,anchor,tuple(sorted(members)),row))
    taken=set()
    labels=[]
    for *_,row in sorted(candidates,key=lambda item:item[:-1]):
        if taken.isdisjoint(row["atom_ids"]):
            labels.append(row)
            taken.update(row["atom_ids"])
    return labels
