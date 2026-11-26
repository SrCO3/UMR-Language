"""Local geometric hypotheses for an electron-executed addition or elimination."""
from rdkit import Chem
from ..chemistry import FormalError
from .smiles import to_rdkit

KEY = 'local_geometry_hypotheses_v1'
BASIS = 'enumerated_geometry_not_selectivity_prediction'
SUBSTITUTION_KEY = 'substitution_geometry_centers_v1'


def neighbours(view):
    adjacent = {a.atom_id:{} for a in view.atoms}
    for b in view.bonds:
        adjacent[b.a][b.b] = adjacent[b.b][b.a] = b.order
    return adjacent


def choices(source, target, substitution_ids=()):
    """Return both parities only for a newly tetrahedral carbon or a new stereogenic pi bond."""
    before, after = neighbours(source), neighbours(target)
    old = {a.atom_id:a for a in source.atoms}
    new = {a.atom_id:a for a in target.atoms}
    mol = to_rdkit(target)
    ids = [a.GetAtomMapNum() for a in mol.GetAtoms()]
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    Chem.SanitizeMol(mol)
    source_mol=to_rdkit(source)
    source_ids=[a.GetAtomMapNum() for a in source_mol.GetAtoms()]
    for atom in source_mol.GetAtoms():atom.SetAtomMapNum(0)
    Chem.SanitizeMol(source_mol)
    old_potential={source_ids[s.centeredOn] for s in Chem.FindPotentialStereo(source_mol)
                   if s.type==Chem.StereoType.Atom_Tetrahedral}
    changed={i for i in before if before[i]!=after[i]}
    alternatives = []
    for s in Chem.FindPotentialStereo(mol, cleanIt=False, flagPossible=True):
        if s.specified != Chem.StereoSpecified.Unspecified:
            continue
        if s.type == Chem.StereoType.Atom_Tetrahedral:
            i = ids[s.centeredOn]
            a, b = old[i], new[i]
            lost=set(before[i])-set(after[i])
            acquired=set(after[i])-set(before[i])
            copper_coupling=(len(lost)==len(acquired)==1 and old[next(iter(lost))].element=='Cu' and
                             new[next(iter(acquired))].element in {'C','Si'})
            trigonal=(len(before[i])+a.hydrogens==3 and (2 in before[i].values() or a.formal_charge in (-1,1)))
            substitution=(i in substitution_ids and len(lost)==len(acquired)==1 and
                          len(before[i])+a.hydrogens==4 and all(v==1 for v in before[i].values()) and
                          old[next(iter(lost))].element in {'O','N','S','F','Cl','Br','I'} and
                          new[next(iter(acquired))].element in {'C','N','O','S'})
            desym=(i not in old_potential and before[i]==after[i] and a.hydrogens==b.hydrogens and
                   any(j in changed for j in before[i]))
            if (a.element!='C' or b.formal_charge!=0 or b.radical_electrons!=0 or
                    (a.stereo.get('status')=='specified' and not substitution) or b.stereo.get('status')=='specified' or
                    not (trigonal or copper_coupling or substitution or desym) or len(after[i])+b.hydrogens!=4 or
                    not all(v==1 for v in after[i].values())):
                continue
            refs = sorted(after[i])
            kind=('saturated_substitution' if substitution else 'tetrahedral_desymmetrization' if desym
                  else 'copper_carbon_coupling' if copper_coupling else 'tetrahedral_addition')
            alternatives.append([{'kind':kind,'atom_ids':[i],'reference_atoms':refs,
                                  'orientation':p} for p in ('CW','CCW')])
        elif s.type == Chem.StereoType.Bond_Double:
            bond = mol.GetBondWithIdx(s.centeredOn)
            a,b = ids[bond.GetBeginAtomIdx()],ids[bond.GetEndAtomIdx()]
            if new[a].element!='C' or new[b].element!='C' or after[a].get(b)!=2:
                continue
            kind='alkene_formation'
            if before[a].get(b)==2:
                # Replacing one of two identical geminal halogens creates
                # E/Z information even though the C=C order is unchanged.
                # Enumerate both orientations; never select with a reference.
                desym=False
                for center in (a,b):
                    hs=[i for i in before[center] if old[i].element in {'Cl','Br','I'}]
                    if len(hs)==2 and (old[hs[0]].element,old[hs[0]].isotope)==(old[hs[1]].element,old[hs[1]].isotope):
                        missing=set(before[center])-set(after[center]);arriving=set(after[center])-set(before[center])
                        if len(missing)==len(arriving)==1 and missing<=set(hs) and all(new[i].element in {'Pd','Ni','Cu'} for i in arriving):desym=True
                if not desym:continue
                kind='vinyl_halide_desymmetrization'
            elif before[a].get(b)!=1:continue
            ends = [sorted(set(after[a])-{b}),sorted(set(after[b])-{a})]
            if not all(ends):
                continue
            alternatives.append([{'kind':kind,'atom_ids':[a,b],
                                  'reference_atoms':[ends[0][0],ends[1][0]],'orientation':p}
                                 for p in ('cis','trans')])
    return alternatives


def mark(view, branch):
    ids, refs, orientation = branch['atom_ids'],branch['reference_atoms'],branch['orientation']
    if branch['kind'] in {'tetrahedral_addition','tetrahedral_desymmetrization','saturated_substitution','copper_carbon_coupling','retained_initial_tetrahedron'}:
        atom = next(a for a in view.atoms if a.atom_id==ids[0])
        atom.stereo={'status':'specified','rdkit_tag':'CHI_TETRAHEDRAL_'+orientation,
                     'neighbor_order':list(refs),'basis':BASIS}
    else:
        bond=next(b for b in view.bonds if {b.a,b.b}==set(ids))
        bond.stereo={'status':'specified','rdkit_stereo':'STEREO'+orientation.upper(),
                     'reference_atoms':refs if bond.a==ids[0] else list(reversed(refs)),'basis':BASIS}


def apply(source,target,event,report):
    branches=event.conditions[KEY]
    da=(event.event_type=='PERICYCLIC' and event.template_id=='DA_4_2' and event.direction=='forward' and
        len(event.bindings.get('diene',[]))==4 and len(event.bindings.get('dienophile',[]))==2 and
        report.get('atomic_commit') and len(report.get('compiled_actions',[]))==3)
    polar=event.actions and event.event_type in {'POLAR','E2','CARBONYL_ADDITION','E1_ELIMINATION','E1CB_ELIMINATION','PROTON_TRANSFER'}
    if not isinstance(branches,list) or not branches or len(branches)>4 or not (da or polar):
        raise FormalError('Geometry hypothesis needs an eligible electron-executed polar event')
    allowed=choices(source,target,event.conditions.get(SUBSTITUTION_KEY,()))
    used=set()
    for branch in branches:
        def equivalent(option):
            if not isinstance(branch,dict) or set(branch)!=set(option):
                return False
            if branch['orientation']!=option['orientation'] or branch['kind']!=option['kind']:
                return False
            ids,refs=branch['atom_ids'],branch['reference_atoms']
            if not isinstance(ids,list) or not isinstance(refs,list) or not all(type(i) is int for i in [*ids,*refs]):
                return False
            if branch['kind'] in {'tetrahedral_addition','tetrahedral_desymmetrization','saturated_substitution','copper_carbon_coupling'}:
                return ids==option['atom_ids'] and sorted(refs)==sorted(option['reference_atoms'])
            if set(ids)!=set(option['atom_ids']) or len(ids)!=2 or len(refs)!=2:
                return False
            adj=neighbours(target)
            return refs[0] in adj[ids[0]] and refs[1] in adj[ids[1]] and not set(refs).intersection(ids)
        if not any(equivalent(option) for options in allowed for option in options):
            raise FormalError('Unproved or nonlocal stereogenic site or invalid local geometric reference')
        key=tuple(sorted(branch['atom_ids']))
        if key in used:
            raise FormalError('Duplicate local geometric constraint')
        used.add(key)
        mark(target,branch)
    report['local_geometry']={'version':KEY,'hypotheses':branches,'basis':BASIS,
                              'chemical_selectivity_verified':False,'origin':'proposed_completion'}
    if da:report['local_geometry']['pericyclic_geometry']='enumerated hypotheses; suprafacial correlations and endo/exo preference not established'

