"""Fixed, map-independent 2D projections for isolated carbon cage systems."""
from functools import lru_cache
from math import dist
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Geometry import Point2D, Point3D

# Reference order is the atom order of each SMILES. Coordinates are drawing
# conventions in bond-length units. They carry no stereochemical assignment.
REFERENCES = (
    ("adamantane","C1C2CC3CC1CC(C2)C3",
     [(-1.95,1.05),(0,1.85),(1.95,1.05),(1.78,-.82),(.05,-2.05),(-1.78,-.82),(-.82,-.95),(0,-.12),(-.62,.84),(1.06,.10)]),
    ("bornane","CC1(C)C2CCC1(C)CC2",
     [(-.70,2.30),(0,1.10),(.83,2.18),(-1.60,.02),(-1.03,-1.72),(.73,-1.75),(1.58,-.02),(2.92,.48),(.70,-.55),(-.79,-.66)]),
    ("pinane","CC1CCC2CC1C2(C)C",
     [(.84,1.75),(1.17,.93),(1.58,-.52),(1.11,-1.59),(-.02,-1.08),(.46,.21),(0,.95),(-1.13,-.17),(-1.65,-.59),(-2.37,.10)]),
    ("bicyclo[2.2.1]heptane","C1CC2CCC1C2",
     [(1.39,-.39),(.98,1.09),(-.29,1.06),(-1.43,.41),(-1.01,-1.07),(.32,-1.08),(.05,-.03)]),
    ("bicyclo[3.1.1]heptane","C1CC2CC(C1)C2",
     [(-.14,-1.2),(1.16,-.97),(1.1,.36),(.01,.33),(-.98,.65),(-1.37,-.62),(.21,1.45)]),
)


@lru_cache(maxsize=1)
def patterns():
    rows = []
    for name,smiles,points in REFERENCES:
        base = Chem.MolFromSmiles(smiles)
        # Match the carbon skeleton in single/double bond derivatives.
        options = Chem.AdjustQueryParameters.NoAdjustments()
        options.makeBondsGeneric = True
        rows.append((name,Chem.AdjustQueryProperties(base,options),points))
    return rows


def depict(mol, enabled):
    rdDepictor.Compute2DCoords(mol,canonOrient=True,useRingTemplates=enabled)
    if not enabled:
        return []
    Chem.GetSymmSSSR(mol)
    if mol.GetRingInfo().NumRings()<2:
        return []
    mappings=[]
    frags=Chem.GetMolFrags(mol,asMols=True,sanitizeFrags=False,fragsMolAtomMapping=mappings)
    dest=mol.GetConformer()
    applied=[]
    for frag,indices in zip(frags,mappings):
        Chem.GetSymmSSSR(frag)
        if frag.GetRingInfo().NumRings()<2:
            continue
        for name,query,points in patterns():
            matches=[]
            for match in frag.GetSubstructMatches(query,uniquify=False,maxMatches=512):
                members=set(match)
                ring_edges=[e for e in frag.GetBonds() if e.IsInRing() and
                            (e.GetBeginAtomIdx() in members or e.GetEndAtomIdx() in members)]
                if any(e.GetBeginAtomIdx() not in members or e.GetEndAtomIdx() not in members
                       or e.GetBondTypeAsDouble() not in (1.,2.) for e in ring_edges):
                    continue
                if sum(1 for e in frag.GetBonds() if e.GetBeginAtomIdx() in members and e.GetEndAtomIdx() in members)!=query.GetNumBonds():
                    continue
                matches.append(match)
            if not matches:
                continue
            center=tuple(sum(p[k] for p in points)/len(points) for k in (0,1))
            def access_key(match):
                # Put external substituents at accessible perimeter positions.
                # Symmetry is resolved consistently with stable UMR atom IDs.
                members=set(match)
                exposed=sum(dist(points[k],center) for k,i in enumerate(match)
                            if any(a.GetIdx() not in members for a in frag.GetAtomWithIdx(i).GetNeighbors()))
                return -round(exposed,6),tuple(frag.GetAtomWithIdx(i).GetAtomMapNum() for i in match)
            match=min(matches,key=access_key)
            fixed={i:Point2D(*p) for i,p in zip(match,points)}
            rdDepictor.Compute2DCoords(frag,canonOrient=False,coordMap=fixed,forceRDKit=True)
            pos=frag.GetConformer()
            error=max(dist((pos.GetAtomPosition(i).x,pos.GetAtomPosition(i).y),p) for i,p in zip(match,points))
            if error>1e-4:
                raise ValueError("RDKit did not preserve locked cage coordinates")
            # Restore the fragment's previous centroid (important for eta assemblies).
            dx=sum(dest.GetAtomPosition(i).x-pos.GetAtomPosition(j).x for j,i in enumerate(indices))/len(indices)
            dy=sum(dest.GetAtomPosition(i).y-pos.GetAtomPosition(j).y for j,i in enumerate(indices))/len(indices)
            for j,i in enumerate(indices):
                p=pos.GetAtomPosition(j)
                dest.SetAtomPosition(i,Point3D(p.x+dx,p.y+dy,0))
            applied.append(dict(template=name,atom_ids=[frag.GetAtomWithIdx(i).GetAtomMapNum() for i in match],
                                max_locked_error=error))
            break
    return applied
