"""Rigid, reaction-site-aware placement of disconnected molecular assemblies."""
import math
from rdkit.Geometry import Point3D
from .geometry import segment_distance,intersects


def centroid(points):
    return tuple(sum(p[k] for p in points)/len(points) for k in (0,1))


def limits(points):
    return (min(x for x,y in points),min(y for x,y in points),
            max(x for x,y in points),max(y for x,y in points))


def place_partners(mol, components, atom_ids, pairs, expansion):
    if len(components)<2 or not pairs:
        return []
    conf=mol.GetConformer()
    groups=[{atom_ids[i]:i for i in c} for c in components]
    heavy=lambda g:sum(mol.GetAtomWithIdx(i).GetAtomicNum()!=1 for i in g.values())
    positions={a:(conf.GetAtomPosition(i).x,conf.GetAtomPosition(i).y) for g in groups for a,i in g.items()}
    connections=sorted({tuple(sorted((a,b))) for a,b in pairs if any(a in g and b not in g for g in groups)
                        and a in positions and b in positions})
    if not connections:
        return []
    groups.sort(key=lambda g:(not any(a in g or b in g for a,b in connections),-heavy(g),min(g)))
    # Orient two comparable partners with their reacting faces towards one
    # another using proper rotations only (no stereochemical reflection).
    if len(groups)==2 and min(map(heavy,groups))>=5 and 2*min(map(heavy,groups))>=max(map(heavy,groups)):
        for k,g in enumerate(groups):
            mid=centroid([positions[a] for a in g])
            site=centroid([positions[a] for pair in connections for a in pair if a in g])
            angle=k*math.pi-math.atan2(site[1]-mid[1],site[0]-mid[0])
            co,si=math.cos(angle),math.sin(angle)
            for a,i in g.items():
                x,y=positions[a][0]-mid[0],positions[a][1]-mid[1]
                positions[a]=(mid[0]+co*x-si*y,mid[1]+si*x+co*y)
                conf.SetAtomPosition(i,Point3D(*positions[a],0))
    placed=set(groups[0])
    boxes=[limits([positions[a] for a in placed])]
    records=[]
    for n in range(1,len(groups)):
        # Follow the contact graph first, including networks with three partners.
        k=next((k for k in range(n,len(groups)) if any(
            (a in placed and b in groups[k]) or (b in placed and a in groups[k]) for a,b in connections)),n)
        groups[n],groups[k]=groups[k],groups[n]
        group=groups[n]
        edges=[(a,b) if a in placed else (b,a) for a,b in connections
               if (a in placed and b in group) or (b in placed and a in group)]
        points=[positions[a] for a in group]
        donor=centroid([positions[b] for a,b in edges]) if edges else centroid(points)
        site=centroid([positions[a] for a,b in edges]) if edges else (max(b[2] for b in boxes),0)
        large=heavy(group)>=5 and heavy(group)*2>=heavy(groups[0])
        dxs,dys=zip(*[(x-donor[0],y-donor[1]) for x,y in points])
        if large or not edges:
            gap=3.9*expansion
            point=(max(b[2] for b in boxes)+gap-min(dxs),site[1])
            mode="horizontal_partners" if edges else "inventory_row"
        else:
            # Occupancy includes the whole assembly plus room for charges and
            # labels. Enumerate its four outer edges; no maximum-radius cutoff.
            gap=2.7*expansion
            candidates=[]
            xs=[site[0]];ys=[site[1]]
            for x0,y0,x1,y1 in boxes:
                xs.extend((x0-gap-max(dxs),x1+gap-min(dxs)))
                ys.extend((y0-gap-max(dys),y1+gap-min(dys)))
            for x in xs:
                for y in ys:
                    if all(x+max(dxs)<=u-gap+1e-6 or x+min(dxs)>=w+gap-1e-6 or
                           y+max(dys)<=v-gap+1e-6 or y+min(dys)>=z+gap-1e-6 for u,v,w,z in boxes):
                        candidates.append((sum(math.dist(positions[a],(positions[b][0]-donor[0]+x,positions[b][1]-donor[1]+y)) for a,b in edges),x,y))
            trial={(x,y) for _,x,y in candidates}
            for radius in (2.6,3.15,3.85,4.7):
                for j in range(16):
                    angle=j*math.pi/8
                    trial.add((site[0]+radius*expansion*math.cos(angle),site[1]+radius*expansion*math.sin(angle)))
            bonds=[(atom_ids[e.GetBeginAtomIdx()],atom_ids[e.GetEndAtomIdx()]) for e in mol.GetBonds()]
            bonds=[(a,b) for a,b in bonds if a in placed and b in placed]
            def merit(point):
                wx,wy=point[0]-donor[0],point[1]-donor[1]
                moved={i:(positions[i][0]+wx,positions[i][1]+wy) for i in group}
                clearance=min(math.dist(p,positions[i]) for p in moved.values() for i in placed)
                # Use a conservative envelope for molecular reagents. Single
                # atoms may occupy an open sector beyond an oblique ring edge.
                overlap=len(group)>1 and any(not (point[0]+max(dxs)<u-1.4 or point[0]+min(dxs)>w+1.4 or
                    point[1]+max(dys)<v-1.4 or point[1]+min(dys)>z+1.4) for u,v,w,z in boxes)
                crossings=0
                for a,b in edges:
                    for i,j in bonds:
                        if min(segment_distance(q,positions[i],positions[j]) for q in moved.values())<1.45:
                            overlap=True
                        if a!=i and a!=j:
                            crossings+=intersects(positions[a],moved[b],positions[i],positions[j])
                    crossings+=sum(segment_distance(positions[i],positions[a],moved[b])<.62
                                   for i in placed if i!=a)
                length=sum(math.dist(positions[a],moved[b]) for a,b in edges)
                return (int(overlap or clearance<2.05),crossings,round(length,5),-round(clearance,4),point)
            x,y=min(trial,key=merit)
            point=(x,y)
            mode="small_partner_at_site"
        shift=(point[0]-donor[0],point[1]-donor[1])
        for a,i in group.items():
            p=(positions[a][0]+shift[0],positions[a][1]+shift[1])
            conf.SetAtomPosition(i,Point3D(*p,0))
            positions[a]=p
        boxes.append(limits([positions[a] for a in group]))
        placed.update(group)
        records.append(dict(method=mode,atom_ids=sorted(group),contact_pairs=edges,
                            site_distances=[math.dist(positions[a],positions[b]) for a,b in edges]))
    components[:]=[tuple(g.values()) for g in groups]
    return records
