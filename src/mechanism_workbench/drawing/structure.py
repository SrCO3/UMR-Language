"""Stage 1: existing UMR -> RDKit SVG and explicit view-scoped coordinates."""
import math,re
from xml.etree import ElementTree as ET

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D

from .geometry import add,distance,mul,normal,sub,unit
from .svg import element,label,text_box


def native_signature(mol):
    # WedgeMolBonds can reverse the storage direction of an ordinary bond when
    # placing its narrow end at a stereocenter. Covalent edges are undirected;
    # dative bonds retain their donor->acceptor direction. Chiral tags stay checked.
    def ends(b):
        pair=(b.GetBeginAtomIdx(),b.GetEndAtomIdx())
        return pair if str(b.GetBondType())=='DATIVE' else tuple(sorted(pair))
    return [(a.GetAtomMapNum(),a.GetSymbol(),a.GetIsotope(),a.GetFormalCharge(),a.GetNumExplicitHs(),a.GetNumRadicalElectrons(),str(a.GetChiralTag())) for a in mol.GetAtoms()], [(*ends(b),str(b.GetBondType()),str(b.GetStereo())) for b in mol.GetBonds()]


def spread_triphenylphosphines(mol,view):
    """Rigidly move complete monodentate PPh3 ligands away from crowded metals.

    This is depiction geometry, not a 3D coordination/stereochemistry claim.
    No atom, edge, aromaticity, chiral tag or authoritative coordinate is edited.
    Bridging/chelating or substituted phosphines remain outside this small fix.
    """
    ids={a.GetAtomMapNum():a.GetIdx() for a in mol.GetAtoms()};conf=mol.GetConformer()
    ligands={l.ligand_id:l for l in getattr(view,'ligands',[])};moved=[]
    for b in getattr(view,'metal_bindings',[]):
        if b.mode!='dative' or len(b.donor_atom_ids)!=1:continue
        donor=b.donor_atom_ids[0];lig=ligands.get(b.ligand_id)
        if not lig or len(lig.atom_ids)!=19 or mol.GetAtomWithIdx(ids[donor]).GetSymbol()!='P':continue
        members=set(lig.atom_ids)
        if any(mol.GetAtomWithIdx(ids[i]).GetSymbol()!='C' for i in members-{donor}):continue
        # Never separate a ligand that has another external contact.
        contacts=[edge for edge in mol.GetBonds() if (edge.GetBeginAtom().GetAtomMapNum() in members)!=(edge.GetEndAtom().GetAtomMapNum() in members)]
        if len(contacts)!=1 or str(contacts[0].GetBondType())!='DATIVE':continue
        p=conf.GetAtomPosition(ids[donor]);m=conf.GetAtomPosition(ids[b.metal_atom_id])
        dx,dy=p.x-m.x,p.y-m.y;length=math.hypot(dx,dy)
        if length<.1:continue
        factor=max(0,6.5/length-1)
        for i in members:
            old=conf.GetAtomPosition(ids[i]);conf.SetAtomPosition(ids[i],Point3D(old.x+dx*factor,old.y+dy*factor,old.z))
        moved.append(dict(ligand_id=lig.ligand_id,atom_ids=sorted(members),method='rigid_translation_extended_dative_bond'))
    return moved


def spread_terminal_metal_ligands(mol,view):
    """Display-only spacing for W/Mo alkoxides and axial Co-salen ligands.

    Rigid translation preserves each ligand and stereochemistry. Extended metal
    contacts are schematic: these coordinates do not assert coordination geometry.
    Chelates, bridging ligands and multiply attached groups are never separated.
    """
    from .geometry import segment_distance
    ids={a.GetAtomMapNum():a.GetIdx() for a in mol.GetAtoms()};conf=mol.GetConformer()
    ligands={l.ligand_id:l for l in getattr(view,'ligands',[])};moved=[]
    bindings=getattr(view,'metal_bindings',[])
    for b in bindings:
        if len(b.donor_atom_ids)!=1 or b.mode not in {'dative','sigma_covalent'}:continue
        donor=b.donor_atom_ids[0];lig=ligands.get(b.ligand_id);symbol=mol.GetAtomWithIdx(ids[b.metal_atom_id]).GetSymbol()
        if not lig:continue
        cobalt=(symbol=='Co' and any(x.metal_atom_id==b.metal_atom_id and ligands[x.ligand_id].ionic_class=='L2X2' for x in bindings))
        tungsten=(symbol in {'W','Mo'} and mol.GetAtomWithIdx(ids[donor]).GetSymbol()=='O' and lig.ionic_class=='X')
        if not (cobalt or tungsten):continue
        members=set(lig.atom_ids)
        contacts=[e for e in mol.GetBonds() if (e.GetBeginAtom().GetAtomMapNum() in members)!=(e.GetEndAtom().GetAtomMapNum() in members)]
        if len(contacts)!=1 or contacts[0].GetBondTypeAsDouble()!=1:continue
        component=next(c for c in Chem.GetMolFrags(mol) if ids[donor] in c)
        other=[i for i in component if mol.GetAtomWithIdx(i).GetAtomMapNum() not in members and i!=ids[b.metal_atom_id]]
        old={i:conf.GetAtomPosition(ids[i]) for i in members};m=conf.GetAtomPosition(ids[b.metal_atom_id]);p=old[donor]
        candidates=[]
        for radius in (4.,5.5,7.,9.):
            for k in range(32):
                angle=k*math.pi/16;dx=m.x+radius*math.cos(angle)-p.x;dy=m.y+radius*math.sin(angle)-p.y
                penalty=sum(max(0,1.8-math.hypot(q.x+dx-conf.GetAtomPosition(j).x,q.y+dy-conf.GetAtomPosition(j).y))**2 for q in old.values() for j in other)
                # A long display contact must not pass through another atom label.
                penalty+=sum(max(0,.6-segment_distance((conf.GetAtomPosition(j).x,conf.GetAtomPosition(j).y),(m.x,m.y),(p.x+dx,p.y+dy)))**2*4 for j in other)
                candidates.append((penalty*1000+math.hypot(dx,dy)+radius*.2+k*.001,dx,dy))
        _,dx,dy=min(candidates)
        for i,p in old.items():conf.SetAtomPosition(ids[i],Point3D(p.x+dx,p.y+dy,p.z))
        moved.append(dict(ligand_id=lig.ligand_id,atom_ids=sorted(members),method='rigid_translation_schematic_metal_contact'))
    return moved


def spread_phosphine_chelates(mol,view):
    """Depict a C/H/P diphosphine with all of its actual donor contacts.

    RDKit's metal-containing ring layout can pile BINAP phenyl rings on top of
    each other. Generate the disconnected ligand's coordinates, rigidly place
    them next to the metal, then draw the unchanged original dative graph.
    Bridged/multimetal and other ligand classes are outside this display fix.
    """
    from .geometry import segment_distance,intersects
    ids={a.GetAtomMapNum():a.GetIdx() for a in mol.GetAtoms()};conf=mol.GetConformer();moved=[]
    for lig in getattr(view,'ligands',[]):
        contacts=[b for b in getattr(view,'metal_bindings',[]) if b.ligand_id==lig.ligand_id]
        donors=sorted(set(lig.donor_atom_ids))
        bound_donors=sorted({i for b in contacts for i in b.donor_atom_ids})
        metals={b.metal_atom_id for b in contacts};members=set(lig.atom_ids)
        if members==set(donors):continue
        if len(contacts)>2 or len(donors)!=2 or len(metals)>1 or any(b.mode!='dative' for b in contacts):continue
        if any(mol.GetAtomWithIdx(ids[i]).GetSymbol()!='P' for i in donors):continue
        if any(mol.GetAtomWithIdx(ids[i]).GetSymbol() not in {'C','H','P'} for i in members):continue
        boundary=[b for b in mol.GetBonds() if (b.GetBeginAtom().GetAtomMapNum() in members)!=(b.GetEndAtom().GetAtomMapNum() in members)]
        if len(boundary)!=len(contacts) or any(str(b.GetBondType())!='DATIVE' for b in boundary):continue
        submol=Chem.RWMol(mol)
        for i in reversed(range(mol.GetNumAtoms())):
            if mol.GetAtomWithIdx(i).GetAtomMapNum() not in members:submol.RemoveAtom(i)
        submol.UpdatePropertyCache(strict=False)
        # CoordGen separates the two fused-ring halves of BINAP; the general
        # metal-ring layout and default free-ligand depiction can overlap them.
        # Its coordinates are used only on this private, disconnected copy.
        from rdkit.Chem import rdCoordGen
        submol.RemoveAllConformers();rdCoordGen.AddCoords(submol)
        # Place terminal phenyl substituents separately: a free BINAP depiction
        # can itself overlap these rings. Ring geometry is moved rigidly and
        # the original P-C edges remain; only display bond lengths may grow.
        sc=submol.GetConformer();idx={a.GetAtomMapNum():a.GetIdx() for a in submol.GetAtoms()}
        for donor in donors:
            pi=idx[donor]
            for neighbor in submol.GetAtomWithIdx(pi).GetNeighbors():
                anchor=neighbor.GetIdx();group={anchor};todo=[anchor]
                while todo:
                    i=todo.pop()
                    for a in submol.GetAtomWithIdx(i).GetNeighbors():
                        j=a.GetIdx()
                        if j!=pi and j not in group:group.add(j);todo.append(j)
                if len(group)!=6 or any(submol.GetAtomWithIdx(i).GetSymbol()!='C' for i in group):continue
                if sum(b.GetBeginAtomIdx() in group and b.GetEndAtomIdx() in group for b in submol.GetBonds())!=6:continue
                old={i:(sc.GetAtomPosition(i).x,sc.GetAtomPosition(i).y) for i in group}
                ax,ay=old[anchor];px,py=sc.GetAtomPosition(pi).x,sc.GetAtomPosition(pi).y
                heading=math.atan2(sum(y for x,y in old.values())/6-ay,sum(x for x,y in old.values())/6-ax)
                outside={i:(sc.GetAtomPosition(i).x,sc.GetAtomPosition(i).y) for i in range(submol.GetNumAtoms()) if i not in group}
                edges=[(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in submol.GetBonds()]
                trials=[]
                for k in range(32):
                    angle=k*math.pi/16;rot=angle-heading;co,si=math.cos(rot),math.sin(rot)
                    for radius in (1.5,2.5,3.8,5.):
                        tx,ty=px+radius*math.cos(angle),py+radius*math.sin(angle)
                        pos={i:((x-ax)*co-(y-ay)*si+tx,(x-ax)*si+(y-ay)*co+ty) for i,(x,y) in old.items()}
                        allpos={**outside,**pos}
                        penalty=sum(max(0,1.15-math.dist(p,q))**2 for p in pos.values() for j,q in outside.items() if j!=pi)
                        changed=[(a,b) for a,b in edges if a in group or b in group]
                        fixed=[(a,b) for a,b in edges if a not in group and b not in group]
                        penalty+=8*sum(intersects(allpos[a],allpos[b],allpos[c],allpos[d]) for a,b in changed for c,d in fixed if len({a,b,c,d})==4)
                        trials.append((penalty*1000+radius+k*.0001,pos))
                _,best=min(trials,key=lambda x:x[0])
                for i,(x,y) in best.items():sc.SetAtomPosition(i,Point3D(x,y,0))
        local={a.GetAtomMapNum():submol.GetConformer().GetAtomPosition(a.GetIdx()) for a in submol.GetAtoms()}
        local_edges=[(b.GetBeginAtom().GetAtomMapNum(),b.GetEndAtom().GetAtomMapNum()) for b in submol.GetBonds()]
        if not contacts:
            # The free incoming ligand must also be readable. Move only its
            # private depiction; do not invent either P->metal connection.
            dx=sum(conf.GetAtomPosition(ids[i]).x-local[i].x for i in members)/len(members)
            dy=sum(conf.GetAtomPosition(ids[i]).y-local[i].y for i in members)/len(members)
            for i,p in local.items():conf.SetAtomPosition(ids[i],Point3D(p.x+dx,p.y+dy,0))
            moved.append(dict(ligand_id=lig.ligand_id,atom_ids=sorted(members),donors=donors,bound_donors=[],metal=None,
                              method='independent_free_diphosphine_layout_no_new_contacts'))
            continue
        center=[sum(getattr(local[i],k) for i in donors)/2 for k in ('x','y')]
        metal=next(iter(metals));m=conf.GetAtomPosition(ids[metal])
        component=next(c for c in Chem.GetMolFrags(mol) if ids[metal] in c)
        other=[conf.GetAtomPosition(j) for j in component if mol.GetAtomWithIdx(j).GetAtomMapNum() not in members]
        trials=[]
        # Every candidate is a proper rotation + translation, never a reflection.
        for angle in range(16):
            theta=angle*math.pi/8;co,si=math.cos(theta),math.sin(theta)
            xy={i:((p.x-center[0])*co-(p.y-center[1])*si,(p.x-center[0])*si+(p.y-center[1])*co) for i,p in local.items()}
            for direction in range(16):
                phi=direction*math.pi/8
                for radius in (7.,10.,13.):
                    dx,dy=m.x+radius*math.cos(phi),m.y+radius*math.sin(phi)
                    placed={i:(x+dx,y+dy) for i,(x,y) in xy.items()}
                    penalty=sum(max(0,1.6-math.hypot(x-q.x,y-q.y))**2 for x,y in placed.values() for q in other)
                    penalty+=4*sum(max(0,.55-segment_distance(p,(m.x,m.y),placed[d]))**2 for d in bound_donors for i,p in placed.items() if i!=d)
                    penalty+=8*sum(intersects((m.x,m.y),placed[d],placed[a],placed[b]) for d in bound_donors for a,b in local_edges if d not in (a,b))
                    trials.append((penalty*1000+radius+angle*.0001+direction*.00001,placed))
        _,best=min(trials,key=lambda x:x[0])
        for i,(x,y) in best.items():conf.SetAtomPosition(ids[i],Point3D(x,y,0))
        moved.append(dict(ligand_id=lig.ligand_id,atom_ids=sorted(members),donors=donors,bound_donors=bound_donors,metal=metal,
                          method='independent_whole_ligand_layout_preserving_original_dative_contacts'))
    return moved


def keep_mapped_atoms(mol, atom_ids):
    if len(atom_ids)==mol.GetNumAtoms():
        return mol
    copy=Chem.RWMol(mol)
    for i in reversed(range(copy.GetNumAtoms())):
        if copy.GetAtomWithIdx(i).GetAtomMapNum() not in atom_ids:
            copy.RemoveAtom(i)
    result=copy.GetMol()
    result.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(result)
    return result


def structure(view,config,*,touched=(),expansion=1.0,align=None,inactive=(),required_hydrogens=(),
              visible_ids=None,contact_pairs=(),definitions=()):
    from ..benchmark.umr_bridge import digest
    from ..umr.smiles import to_rdkit
    if getattr(view,'metal_centers',[]):from ..umr.v03.adapters import to_rdkit
    if getattr(view,'element_rule_version',None):from ..umr.v05.smiles import to_rdkit
    if hasattr(view,'material_boundaries'):from ..umr.v09.smiles import to_rdkit
    authority=view
    frozen=digest(view.model_dump(mode='json'));mol=to_rdkit(view)
    from .abbreviations import find_labels
    from .materials import subset
    from .placement import place_partners
    from .scaffolds import depict
    labels=[]
    chelate_scaffolds=None
    if visible_ids is not None:
        keep=set(visible_ids)
        if config.abbreviate_groups:
            labels=find_labels(mol,view,keep,set(touched)|set(required_hydrogens),
                [*config.abbreviations,*definitions,*view.annotations.get('drawing_abbreviations',[])])
        for group in labels:
            keep.difference_update(set(group['atom_ids'])-set(group['anchor_atom_ids']))
        view=subset(view,keep)
        # Use the real donor-to-donor path while determining coordinates.
        # Pruning it earlier can put opposite ends of the same small chelate
        # on opposite sides of the metal. No virtual P–P bond is introduced.
        backbone={i for group in labels for i in group.get('layout_backbone_atom_ids',())}
        if backbone-keep:
            mol=keep_mapped_atoms(mol,keep|backbone)
            chelate_scaffolds=depict(mol,config.use_scaffold_templates)
        mol=keep_mapped_atoms(mol,keep)
    labels_by_anchor={row['anchor_atom_id']:row['label'] for row in labels}
    for row in (r for r in labels if len(r['anchor_atom_ids'])>1):
        for i in row['anchor_atom_ids']:
            labels_by_anchor[i]=next(a.element for a in view.atoms if a.atom_id==i)
    if not mol.GetNumAtoms():
        empty=element('g')
        label(empty,'本步产物已记入物料区',260,140,size=16,text_anchor='middle',fill='#62716c')
        return dict(root=empty,coords={},charges=[],hydrogen_display={'hidden_atom_ids':[],'explicit_visible':[],
            'display_hydrogens':[],'attached_hydrogen_labels':[]},width=520,height=280,obstacles=[],components=[],
            atom_indices={},alignment={},mol_coords={},input_hash=frozen,chemical_signature_unchanged=True,
            display_view=view,labels=[],scaffolds=[],partner_layout=[],phosphine_layout=[],terminal_metal_layout=[],
            coordination=[],chelate_links=[])
    from .hydrogens import display_policy
    hydrogen_display=display_policy(view,set(required_hydrogens)|set(touched))
    hidden=set(hydrogen_display['hidden_atom_ids'])
    # Draw the authority's actual Kekule Lewis view, not aromatic dashed bonds.
    original_orders={tuple(sorted((b.a,b.b))):b.order for b in view.bonds}
    for atom in mol.GetAtoms():atom.SetIsAromatic(False)
    for bond in mol.GetBonds():
        key=tuple(sorted((bond.GetBeginAtom().GetAtomMapNum(),bond.GetEndAtom().GetAtomMapNum())))
        if key in original_orders:
            bond.SetIsAromatic(False);bond.SetBondType({1:Chem.BondType.SINGLE,2:Chem.BondType.DOUBLE,3:Chem.BondType.TRIPLE}[original_orders[key]])
    signature=native_signature(mol)
    scaffolds=depict(mol,config.use_scaffold_templates) if chelate_scaffolds is None else chelate_scaffolds
    chelate_layout=spread_phosphine_chelates(mol,view)
    phosphine_layout=spread_triphenylphosphines(mol,view)
    terminal_metal_layout=spread_terminal_metal_ligands(mol,view)
    from .coordination_layout import project_centers,depth_marks
    coordination=project_centers(mol,view)
    # Rigid rotation only: determinant +1, no reflected stereochemical drawings.
    alignment={'method':'stable_atom_id_rigid_rotation_no_reflection','common_atoms':0,'chelate_layout':chelate_layout}
    if align:
        import numpy as np
        ids={a.GetAtomMapNum():a.GetIdx() for a in mol.GetAtoms()};common=sorted(set(ids)&set(align))
        if len(common)>=2:
            conf=mol.GetConformer();a=np.array([[conf.GetAtomPosition(ids[i]).x,conf.GetAtomPosition(ids[i]).y] for i in common]);b=np.array([align[i] for i in common])
            u,_,vt=np.linalg.svd((a-a.mean(axis=0)).T@(b-b.mean(axis=0)));r=u@vt
            if np.linalg.det(r)<0:u[:,-1]*=-1;r=u@vt
            for atom in mol.GetAtoms():
                p=conf.GetAtomPosition(atom.GetIdx());xy=np.array([p.x,p.y])@r;conf.SetAtomPosition(atom.GetIdx(),Point3D(float(xy[0]),float(xy[1]),0))
            alignment['common_atoms']=len(common)
    conf=mol.GetConformer();components=list(Chem.GetMolFrags(mol));atom_ids={a.GetIdx():a.GetAtomMapNum() for a in mol.GetAtoms()}
    # Eta binding joins a coordination assembly without inventing RDKit bonds.
    for binding in getattr(view,'metal_bindings',[]):
        if binding.mode=='eta':
            linked={binding.metal_atom_id,*binding.donor_atom_ids}
            joined=[c for c in components if linked & {atom_ids[i] for i in c}]
            components=[c for c in components if c not in joined]+[tuple(i for c in joined for i in c)]
    components.sort(key=lambda c:(not bool({atom_ids[i] for i in c}&set(touched)),bool({atom_ids[i] for i in c}&set(inactive)),min(atom_ids[i] for i in c)))
    partner_layout=place_partners(mol,components,atom_ids,contact_pairs,expansion) if config.site_layout else []
    boxes=[];cursor_x=cursor_y=row_height=0.0;gap=3.3*expansion
    widths=[]
    for comp in components:
        xs=[conf.GetAtomPosition(i).x for i in comp];widths.append(max(xs)-min(xs)+3)
    lane=max(20*expansion,max(widths,default=4))
    for comp in components:
        pts=[conf.GetAtomPosition(i) for i in comp];xmin,xmax=min(p.x for p in pts),max(p.x for p in pts);ymin,ymax=min(p.y for p in pts),max(p.y for p in pts)
        w,h=max(2.5,xmax-xmin),max(2.5,ymax-ymin)
        if not partner_layout:
            if cursor_x and cursor_x+w>lane:cursor_y+=row_height+gap;cursor_x=row_height=0
            for i in comp:
                p=conf.GetAtomPosition(i);conf.SetAtomPosition(i,Point3D(p.x-xmin+cursor_x,p.y-ymin+cursor_y,0))
        boxes.append({'rdkit_indices':list(comp),'atom_ids':[atom_ids[i] for i in comp],'inactive':bool(set(inactive)&{atom_ids[i] for i in comp})})
        cursor_x+=w+gap;row_height=max(row_height,h)
    xs=[conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())];ys=[conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    scale=config.bond_pixels/1.5
    width=max(520,math.ceil((max(xs)-min(xs)+4)*scale));height=max(280,math.ceil((max(ys)-min(ys)+4)*scale))
    if width*height>config.max_panel_pixels:raise ValueError('Structure panel exceeds pixel budget; record retained, not redrawn with altered chemistry')
    drawer=rdMolDraw2D.MolDraw2DSVG(width,height,-1,-1,True)
    options=drawer.drawOptions();options.fixedFontSize=config.font_pixels;options.minFontSize=config.font_pixels;options.maxFontSize=config.font_pixels
    # Keep fluorine legible on the white print background.
    options.updateAtomPalette({9:(0.0,107/255,101/255)})
    options.fixedBondLength=config.bond_pixels;options.padding=.13;options.clearBackground=False;options.includeRadicals=False
    options.prepareMolsBeforeDrawing=False
    for a in mol.GetAtoms():
        original=next(x for x in view.atoms if x.atom_id==a.GetAtomMapNum())
        if original.element=='*' and original.label:a.SetProp('atomLabel',original.label)
        if original.element=='C' and (original.formal_charge or original.isotope):
            isotope=f'<sup>{original.isotope}</sup>' if original.isotope else ''
            # An isolated [CH3-] has no skeletal bonds from which readers can
            # infer its hydrogens. Keep them in its formula without inventing
            # numbered H atoms or changing the connected-carbon display policy.
            attached=('H'+(f'<sub>{original.hydrogens}</sub>' if original.hydrogens>1 else '')) if a.GetDegree()==0 and original.hydrogens else ''
            options.atomLabels[a.GetIdx()]=isotope+'C'+attached
        if any(c.atom_id==original.atom_id for c in getattr(view,'reactive_centers',[])):
            # Skeletal endpoints normally imply ordinary carbon valence. Explicit
            # CH/CH2 at a low-valent center prevents a carbene from looking like
            # an ordinary saturated end group; mapped transferring H stays separate.
            isotope=f'<sup>{original.isotope}</sup>' if original.isotope else ''
            attached='H'+(f'<sub>{original.hydrogens}</sub>' if original.hydrogens>1 else '') if original.hydrogens else ''
            options.atomLabels[a.GetIdx()]=isotope+original.element+attached
        if original.atom_id in labels_by_anchor:
            options.atomLabels[a.GetIdx()]=labels_by_anchor[original.atom_id]
        a.SetAtomMapNum(0)
    Chem.WedgeMolBonds(mol,mol.GetConformer())
    drawer.DrawMolecule(mol);drawer.FinishDrawing();root=ET.fromstring(drawer.GetDrawingText())
    dative={b.GetIdx() for b in mol.GetBonds() if str(b.GetBondType())=='DATIVE'}
    for node in root.iter():
        if any(f'bond-{i} ' in node.get('class','')+' ' for i in dative):
            node.set('style',re.sub(r'(stroke|fill):#[0-9A-Fa-f]{6}',r'\1:#26739c',node.get('style','')))
            node.set('data-connection-kind','dative')
    coords={atom_ids[i]:[drawer.GetDrawCoords(i).x,drawer.GetDrawCoords(i).y] for i in range(mol.GetNumAtoms())}
    atoms_by_id={a.atom_id:a for a in view.atoms}
    charged={i:atoms_by_id[aid] for i,aid in atom_ids.items() if atoms_by_id[aid].formal_charge}
    for parent in root.iter():
        for node in list(parent):
            indices=[int(i) for i in re.findall(r'atom-(\d+)',node.get('class',''))]
            text=''.join(node.itertext()).strip()
            if indices and any(atom_ids[i] in hidden for i in indices):
                parent.remove(node);continue
            if node.tag.endswith('text') and len(indices)==1 and indices[0] in charged:
                atom=charged[indices[0]]
                if text in {'+','-','−'} or (atom.element!='C' and abs(atom.formal_charge)>1 and text==str(abs(atom.formal_charge)) and float(node.get('x','0'))>coords[atom.atom_id][0]):
                    parent.remove(node)
    # Recheck identity/isotope/H/charge/r/stereo after depiction (bond wedge
    # direction is layout metadata, not a replacement for AST stereochemistry).
    for a in mol.GetAtoms():a.SetAtomMapNum(atom_ids[a.GetIdx()])
    if native_signature(mol)!=signature:raise ValueError('RDKit depiction preparation changed chemical identity or stereo')
    if frozen!=digest(authority.model_dump(mode='json')):raise ValueError('Drawing mutated authoritative UMR')
    obstacles=[]
    for j,node in enumerate(root.iter()):
        tag=node.tag.rsplit('}',1)[-1];klass=node.attrib.get('class','')
        indices=[int(i) for i in re.findall(r'atom-(\d+)',klass)];ids=[atom_ids[i] for i in indices]
        if tag=='text':
            x=float(node.get('x','0'));y=float(node.get('y','0'));style=node.get('style','')
            match=re.search(r'font-size:([\d.]+)px',style);size=float(match.group(1)) if match else config.font_pixels
            text=''.join(node.itertext());box=text_box(text,x,y,size)
            obstacles.append({'id':f'rdkit.text.{j}','kind':'text','box':box,'atom_ids':ids})
        elif tag=='path':
            values=[float(v) for v in re.findall(r'-?\d+(?:\.\d+)?(?:e[+-]?\d+)?',node.get('d',''),re.I)]
            points=list(zip(values[::2],values[1::2]))
            if len(points)>=2:
                match=re.search(r'stroke-width:([\d.]+)',node.get('style',''));stroke=float(match.group(1)) if match else 2
                obstacles.append({'id':f'rdkit.path.{j}','kind':'bond','box':[min(p[0] for p in points)-stroke,min(p[1] for p in points)-stroke,max(p[0] for p in points)+stroke,max(p[1] for p in points)+stroke],
                    'segments':list(zip(points,points[1:])),'width':stroke,'atom_ids':ids,
                    'filled':bool(re.search(r'fill:(?!none)',node.get('style',''))) and 'Z' in node.get('d','').upper()})
    from .figure_keys import chelate_links
    depth_marks(root,coords,coordination)
    links=chelate_links(root,coords,labels,obstacles)
    from .hydrogens import draw_attached_hydrogens
    obstacles+=draw_attached_hydrogens(view,coords,root,obstacles,config,hydrogen_display)
    for box in boxes:
        pts=[coords[i] for i in box['atom_ids'] if i not in hidden]
        pts += [h['local_xy'] for h in hydrogen_display['display_hydrogens'] if h['parent_atom_id'] in box['atom_ids']]
        box['box']=[min(p[0] for p in pts)-25,min(p[1] for p in pts)-30,max(p[0] for p in pts)+25,max(p[1] for p in pts)+30]
    charges=[]
    from .geometry import segment_distance
    for atom in view.atoms:
        if not atom.formal_charge:continue
        p=coords[atom.atom_id];neighbors=[coords[b.b if b.a==atom.atom_id else b.a] for b in view.bonds if atom.atom_id in {b.a,b.b}]
        component=next(b for b in boxes if atom.atom_id in b['atom_ids'])
        center=mul(tuple(sum(coords[i][k] for i in component['atom_ids']) for k in (0,1)),1/len(component['atom_ids']))
        candidates=[]
        for radius in (22,30,40,52):
            for k in range(16):
                angle=-math.pi/2+k*math.pi/8;point=add(p,(radius*math.cos(angle),radius*math.sin(angle)));box=[point[0]-11,point[1]-11,point[0]+11,point[1]+11]
                overlap=sum(1 for o in obstacles if not (box[2]<o['box'][0] or box[0]>o['box'][2] or box[3]<o['box'][1] or box[1]>o['box'][3]))
                toward=sum(max(0,1-(distance(point,q)/max(distance(p,q),1))) for q in neighbors)
                outward=distance(point,center)-distance(p,center)
                score=overlap*10000+toward*100+radius*.8-outward*.2+k*.01
                candidates.append((score,point,box))
        _,point,box=min(candidates,key=lambda v:v[0]);charges.append({'atom_id':atom.atom_id,'formal_charge':atom.formal_charge,'position':point,'box':box})
        obstacles.append({'id':f'charge.{atom.atom_id}','kind':'charge','atom_ids':[atom.atom_id],'box':box,'position':point})
    return {'root':root,'coords':coords,'charges':charges,'hydrogen_display':hydrogen_display,'width':width,'height':height,'obstacles':obstacles,'components':boxes,
        'atom_indices':{atom_ids[i]:i for i in atom_ids},'alignment':alignment,'mol_coords':{atom_ids[i]:[conf.GetAtomPosition(i).x,conf.GetAtomPosition(i).y] for i in atom_ids},
        'input_hash':frozen,'chemical_signature_unchanged':True,'phosphine_layout':phosphine_layout,'terminal_metal_layout':terminal_metal_layout,
        'display_view':view,'labels':labels,'scaffolds':scaffolds,'partner_layout':partner_layout,
        'coordination':coordination,'chelate_links':links}


def place_layout(layout,parent,offset,view,side,event_id,ids_layer,show_ids=False):
    view=layout['display_view']
    group=element('g',transform=f'translate({offset[0]} {offset[1]})',data_view_id=view.view_id,data_side=side)
    for child in list(layout['root']):group.append(child)
    parent.append(group);coords={i:add(p,offset) for i,p in layout['coords'].items()};obstacles=[]
    for o in layout['obstacles']:
        d=dict(o);d['id']=f'{side}.{o["id"]}';x0,y0,x1,y1=o['box'];d['box']=[x0+offset[0],y0+offset[1],x1+offset[0],y1+offset[1]]
        if 'segments' in o:d['segments']=[(add(a,offset),add(b,offset)) for a,b in o['segments']]
        if 'position' in o:d['position']=add(o['position'],offset)
        obstacles.append(d)
    for charge in layout.get('charges',[]):
        p=add(charge['position'],offset);q=charge['formal_charge']
        parent.append(element('circle',cx=p[0],cy=p[1],r=8.5,fill='none',stroke='#244b5c',stroke_width=1.35,**{'class':'formal-charge','data-atom-id':charge['atom_id'],'data-side':side}))
        label(parent,('+' if q>0 else '−') if abs(q)==1 else str(abs(q))+('+' if q>0 else '−'),p[0],p[1]+4,size=13,fill='#244b5c',text_anchor='middle',**{'class':'formal-charge','data-atom-id':charge['atom_id'],'data-side':side})
    neighbors={i:[] for i in coords}
    for b in view.bonds:neighbors[b.a].append(coords[b.b]);neighbors[b.b].append(coords[b.a])
    for i,p in coords.items():
        if i in layout['hydrogen_display']['hidden_atom_ids']:continue
        # Estimate free angular sector, then reserve its ID label rectangle.
        angles=[math.atan2(q[1]-p[1],q[0]-p[0]) for q in neighbors[i]]
        choices=[k*math.pi/4 for k in range(8)]
        angle=max(choices,key=lambda a:min((abs(math.atan2(math.sin(a-b),math.cos(a-b))) for b in angles),default=math.pi))
        x=p[0]+math.cos(angle)*22;y=p[1]+math.sin(angle)*22
        node=label(ids_layer,i,x,y,size=11,fill='#666',text_anchor='middle',data_atom_id=i,data_view_id=view.view_id,data_side=side)
        box=text_box(str(i),x-len(str(i))*3,y,11)
        if show_ids:obstacles.append({'id':f'{side}.id.{i}','kind':'atom_id','atom_ids':[i],'box':box})
    for binding in getattr(view,'metal_bindings',[]):
        if binding.mode=='eta':
            start=coords[binding.metal_atom_id];pts=[coords[i] for i in binding.donor_atom_ids]
            end=tuple(sum(p[k] for p in pts)/len(pts) for k in (0,1))
            start=add(start,mul(unit(sub(end,start)),13))
            parent.append(element('line',x1=start[0],y1=start[1],x2=end[0],y2=end[1],stroke='#26739c',stroke_width=2,stroke_dasharray='5 4',data_binding_id=binding.binding_id,data_side=side))
            note='eta2 (one 2e unit)' if binding.hapticity==2 else 'eta'+str(binding.hapticity)+' (haptic region)';label(parent,note,end[0]+10,end[1]-10,size=11,fill='#26739c',data_side=side)
            obstacles.extend([{'id':side+'.eta.'+binding.binding_id,'kind':'bond','atom_ids':[binding.metal_atom_id,*binding.donor_atom_ids],
                'box':[min(start[0],end[0]),min(start[1],end[1]),max(start[0],end[0]),max(start[1],end[1])],'segments':[(start,end)],'width':2},
                {'id':side+'.eta.label.'+binding.binding_id,'kind':'text','atom_ids':[],'box':text_box(note,end[0]+10,end[1]-10,11)}])
    return coords,obstacles,[{'atom_id':i,'view_id':view.view_id,'state_id':view.state_id,'event_id':event_id,'side':side,
        'rdkit_atom_index':layout['atom_indices'][i],'local_xy':layout['coords'][i],'canvas_xy':list(coords[i]),'visible':i not in layout['hydrogen_display']['hidden_atom_ids'],
        'transform':{'scale':1.0,'translate':list(offset)},'index_scope':f'{event_id}:{view.view_id}:{side}'} for i in coords]
