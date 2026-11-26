"""Projected coordination sites. All changes are on the private RDKit conformer."""
from itertools import permutations
import math

from rdkit import Chem
from rdkit.Geometry import Point3D

from .svg import element

# site name, drawing x/y, depth relative to metal, ideal unit-vector direction.
# Opposite/adjacent tests use the 3D direction, not the angle on the page.
GEOMETRIES = {
    "square_planar": ("平面正方形", [
        ("p1", .82,.82,0,(1,0,0)), ("p2",-.82,.82,0,(0,1,0)),
        ("p3",-.82,-.82,0,(-1,0,0)), ("p4",.82,-.82,0,(0,-1,0))]),
    "tetrahedral": ("正四面体", [
        ("t1",-.95,.60,0,(1,1,1)), ("t2",.95,.60,0,(-1,-1,1)),
        ("t3",-.47,-.83,1,(-1,1,-1)), ("t4",.54,-.73,-1,(1,-1,-1))]),
    "trigonal_bipyramidal": ("三角双锥", [
        ("ax1",0,1.14,0,(0,0,1)), ("ax2",0,-1.14,0,(0,0,-1)),
        ("eq1",1.12,0,0,(1,0,0)), ("eq2",-.8,-.48,1,(-.5,.8660254,0)),
        ("eq3",-.78,.49,-1,(-.5,-.8660254,0))]),
    "square_pyramidal": ("四方锥", [
        ("apical",0,1.24,0,(0,0,1)), ("basal1",.98,.36,-1,(1,0,0)),
        ("basal2",-.89,.48,-1,(0,1,0)), ("basal3",-.98,-.36,1,(-1,0,0)),
        ("basal4",.89,-.48,1,(0,-1,0))]),
    "octahedral": ("八面体", [
        ("ax1",0,1.21,0,(0,0,1)), ("ax2",0,-1.21,0,(0,0,-1)),
        ("eq1",.99,.34,-1,(1,0,0)), ("eq2",-.86,.50,-1,(0,1,0)),
        ("eq3",-.99,-.34,1,(-1,0,0)), ("eq4",.86,-.50,1,(0,-1,0))]),
}


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def unit(v):
    length=math.sqrt(dot(v,v)) or 1.
    return tuple(x/length for x in v)


def assignments(sites,contacts,center,directions):
    index={b.binding_id:k for k,b in enumerate(contacts)}
    names=[s[0] for s in sites]
    fixed={k:names.index(b.site_id) for k,b in enumerate(contacts) if b.site_id in names}
    relations=[(index[r.first_binding_id],index[r.second_binding_id],r.relation)
               for r in center.relative_positions if r.relation!="unknown" and
               r.first_binding_id in index and r.second_binding_id in index]
    same_ligand=[(i,j) for j,b in enumerate(contacts) for i in range(j)
                 if contacts[i].ligand_id==b.ligand_id]
    angles=[[dot(unit(a[4]),unit(b[4])) for b in sites] for a in sites]
    choices=[]
    for order in permutations(range(len(sites))):
        if any(order[i]!=s for i,s in fixed.items()):
            continue
        if any(abs(angles[order[i]][order[j]]-(-1 if rel=="trans" else 0))>.02 for i,j,rel in relations):
            continue
        # Adjacent chelate donors are a depiction preference. A stated trans
        # relation takes precedence, and no assignment is written to the state.
        opposite=sum(angles[order[i]][order[j]] < -.99 for i,j in same_ligand)
        motion=sum(1-dot(directions[i],unit(sites[s][1:3])) for i,s in enumerate(order))
        choices.append((opposite,motion,order))
    return min(choices)[2] if choices else None


def ligand_components(mol,metal_idx):
    groups={}
    remaining=set(range(mol.GetNumAtoms()))-{metal_idx}
    while remaining:
        members=set()
        queue=[next(iter(remaining))]
        while queue:
            i=queue.pop()
            if i not in remaining:
                continue
            remaining.remove(i)
            members.add(i)
            queue.extend(n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors())
        for i in members:
            groups[i]=frozenset(members)
    return groups


def project_centers(mol,view):
    conf=mol.GetConformer()
    indices={a.GetAtomMapNum():a.GetIdx() for a in mol.GetAtoms()}
    centers={c.atom_id for c in getattr(view,"metal_centers",[])}
    reports=[]
    for center in getattr(view,"metal_centers",[]):
        metal=center.atom_id
        bindings=sorted((b for b in view.metal_bindings if b.metal_atom_id==metal),key=lambda b:b.binding_id)
        contacts=[]
        ambiguous_refs=set()
        for b in bindings:
            if b.mode=="dative" and len(b.donor_atom_ids)>1:
                ambiguous_refs.add(b.binding_id)
                contacts.extend(b.model_copy(update={"donor_atom_ids":[i]}) for i in sorted(b.donor_atom_ids))
            else:
                contacts.append(b)
        geometry=center.geometry
        hints=view.annotations.get("drawing_coordination_geometry",{})
        hint=hints.get(str(metal)) if isinstance(hints,dict) else None
        schematic=geometry=="unknown" and ((isinstance(hint,str) and hint in GEOMETRIES) or
                                            (hint is None and len(contacts)==6))
        if schematic:
            geometry=hint or "octahedral"
        record=dict(metal_atom_id=metal,declared=center.geometry,geometry=geometry,sites=[],
                    status="schematic" if schematic else "declared",
                    source="drawing_hint" if schematic and hint else "six_site_convention" if schematic else "metal_center")
        reports.append(record)
        if geometry not in GEOMETRIES:
            record["status"]="unspecified" if geometry=="unknown" else "general_layout"
            continue
        sites=GEOMETRIES[geometry][1]
        related={bid for r in center.relative_positions if r.relation!="unknown" for bid in (r.first_binding_id,r.second_binding_id)}
        if ambiguous_refs & (related | {b.binding_id for b in bindings if b.site_id}):
            record["status"]="grouped_donor_site_unresolved"
            continue
        if len(sites)!=len(contacts) or any(b.mode!="eta" and len(b.donor_atom_ids)!=1 for b in contacts):
            record["status"]="site_count_or_binding_unresolved"
            continue
        mi=indices[metal]
        described={i for b in contacts for i in b.donor_atom_ids}
        if any(n.GetAtomMapNum() not in described for n in mol.GetAtomWithIdx(mi).GetNeighbors()):
            record["status"]="additional_contacts_unresolved"
            continue
        if mol.GetAtomWithIdx(mi).GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED:
            record["status"]="native_stereochemical_layout"
            continue
        groups=ligand_components(mol,mi)
        metal_indices={indices[i] for i in centers if i in indices}
        arms=[]
        for b in contacts:
            donors=[indices[i] for i in b.donor_atom_ids]
            group=groups.get(donors[0],frozenset())
            if not group or not set(donors)<=group or group & metal_indices:
                break
            arms.append((group,donors))
        if len(arms)!=len(sites):
            record["status"]="shared_or_multinuclear_ligand"
            continue
        positions={i:(conf.GetAtomPosition(i).x,conf.GetAtomPosition(i).y) for i in range(mol.GetNumAtoms())}
        origin=positions[mi]
        def mean(points):
            return tuple(sum(p[k] for p in points)/len(points) for k in (0,1))
        anchors=[mean([positions[i] for i in ids]) for _,ids in arms]
        directions=[unit((p[0]-origin[0],p[1]-origin[1])) for p in anchors]
        order=assignments(sites,contacts,center,directions)
        if order is None:
            record["status"]="incompatible_declared_sites"
            continue
        bundles={}
        for k,(group,_) in enumerate(arms):
            bundles.setdefault(group,[]).append(k)
        xy=[sites[s][1:3] for s in order]
        radii=[math.dist(anchors[a],anchors[b])/math.dist(xy[a],xy[b])
               for ks in bundles.values() for j,b in enumerate(ks) for a in ks[:j]]
        radius=sum(radii)/len(radii) if radii else 2.35
        targets=[(origin[0]+radius*x,origin[1]+radius*y) for x,y in xy]
        moves={}
        for group,ks in bundles.items():
            old=mean([anchors[k] for k in ks])
            new=mean([targets[k] for k in ks])
            if len(ks)==1:
                k=ks[0]
                center_of_group=mean([positions[i] for i in group])
                heading=(center_of_group[0]-old[0],center_of_group[1]-old[1])
                if math.hypot(*heading)<.1:
                    heading=directions[k]
                angle=math.atan2(xy[k][1],xy[k][0])-math.atan2(heading[1],heading[0])
            else:
                pairs=[((anchors[k][0]-old[0],anchors[k][1]-old[1]),
                        (targets[k][0]-new[0],targets[k][1]-new[1])) for k in ks]
                angle=math.atan2(sum(a[0]*b[1]-a[1]*b[0] for a,b in pairs),sum(dot(a,b) for a,b in pairs))
            co,si=math.cos(angle),math.sin(angle)
            def move(point):
                x,y=point[0]-old[0],point[1]-old[1]
                return (new[0]+x*co-y*si,new[1]+x*si+y*co)
            if any(math.dist(move(anchors[k]),targets[k])>.25 for k in ks):
                break
            moves.update((i,move(positions[i])) for i in group)
        else:
            for i,(x,y) in moves.items():
                conf.SetAtomPosition(i,Point3D(x,y,0))
            orders={frozenset((b.a,b.b)):b.order for b in view.bonds}
            record["sites"]=[dict(binding_id=b.binding_id,donors=list(b.donor_atom_ids),mode=b.mode,
                                  order=orders.get(frozenset((metal,b.donor_atom_ids[0])),1),
                                  position=sites[s][0],depth=sites[s][3]) for b,s in zip(contacts,order)]
            record["site_assignment"]="specified" if all(b.site_id in {s[0] for s in sites} for b in contacts) else "schematic"
            continue
        record["status"]="rigid_chelate_fit_unresolved"
    return reports


def depth_marks(root,coords,records):
    """Draw depth on the middle of real contacts; dative arrowheads remain visible."""
    for r in records:
        metal=r["metal_atom_id"]
        for s in r["sites"]:
            if not s["depth"] or s["mode"]=="eta" or s["order"]!=1:
                continue
            m=coords[metal]
            donor=s["donors"][0]
            p=coords[donor]
            u=unit((p[0]-m[0],p[1]-m[1]))
            perp=(-u[1],u[0])
            def point(t,w=0):
                return (m[0]+(p[0]-m[0])*t+perp[0]*w,m[1]+(p[1]-m[1])*t+perp[1]*w)
            color="#26739c" if s["mode"]=="dative" else "#252f30"
            g=element("g",data_binding_id=s["binding_id"],data_depth="front" if s["depth"]>0 else "back")
            root.append(g)
            if s["depth"]>0:
                points=" ".join(f"{x:.2f},{y:.2f}" for x,y in (point(.32),point(.69,3.1),point(.69,-3.1)))
                g.append(element("polygon",points=points,fill=color))
            else:
                a,b=point(.32),point(.69)
                g.append(element("line",x1=a[0],y1=a[1],x2=b[0],y2=b[1],stroke="white",stroke_width=3))
                for k in range(6):
                    t=.32+k*.074
                    a,b=point(t,.6+k*.49),point(t,-.6-k*.49)
                    g.append(element("line",x1=a[0],y1=a[1],x2=b[0],y2=b[1],stroke=color,stroke_width=1.35))
