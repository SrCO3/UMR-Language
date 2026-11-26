"""Preserve explicit Lewis bonds before RDKit cleanup; retain spin in sidecars."""
import hashlib
from rdkit import Chem, rdBase
from ..ast import Bond
from ..chemistry import FormalError, Unsupported, bond_orders
from ..smiles import to_rdkit as shared_rdkit
from .ast import Atom, ChemicalState, LewisView, MetalCenter, Ligand, MetalBinding, LewisAdduct, ElementDeclaration
from .elements import MAIN, METAL_GROUP, VALENCE_ELECTRONS, oxo_oxidation

FLAGS=Chem.SanitizeFlags.SANITIZE_ALL
for _flag in ('SANITIZE_CLEANUP','SANITIZE_CLEANUP_ORGANOMETALLICS','SANITIZE_PROPERTIES',
              'SANITIZE_FINDRADICALS','SANITIZE_SETAROMATICITY','SANITIZE_ADJUSTHS','SANITIZE_CLEANUPCHIRALITY'):
    FLAGS ^= getattr(Chem.SanitizeFlags,_flag)


def import_mapped_smiles(smiles,state_id='s0',view_id='v0',labels=None,metadata=None,*,_stereo_perception=False,
                         _preserve_source_hydrogens=False,_valence=VALENCE_ELECTRONS,_metals=METAL_GROUP,_main=MAIN):
    VALENCE_ELECTRONS, METAL_GROUP, MAIN = _valence, _metals, _main
    params=Chem.SmilesParserParams();params.removeHs=False;params.sanitize=False
    mol=Chem.MolFromSmiles(smiles,params)
    if mol is None:raise FormalError('Invalid mapped SMILES')
    # Kekulize(clearAromaticFlags=True) can move [nH]'s explicit H into
    # RDKit's implicit-H bookkeeping. UMR must retain the input bracket count,
    # not reconstruct it from RDKit's post-normalization valence defaults.
    source_hydrogens={a.GetIdx():a.GetNumExplicitHs() for a in mol.GetAtoms()}
    mol.UpdatePropertyCache(strict=False)
    try:Chem.SanitizeMol(mol,sanitizeOps=FLAGS);Chem.Kekulize(mol,clearAromaticFlags=True)
    except Exception as e:raise Unsupported('Integer Lewis graph could not be preserved: '+type(e).__name__) from e
    if _preserve_source_hydrogens:
        for a in mol.GetAtoms():
            a.SetNumExplicitHs(source_hydrogens[a.GetIdx()])
            a.SetNoImplicit(True)
        mol.UpdatePropertyCache(strict=False)
    if _stereo_perception:Chem.AssignStereochemistry(mol,cleanIt=False,force=True)
    ids=[a.GetAtomMapNum() for a in mol.GetAtoms()]
    if min(ids,default=0)<=0 or len(set(ids))!=len(ids):raise FormalError('Every explicit atom needs a unique positive atom_id')
    bonds=[];datives=[];atoms=[]
    for b in mol.GetBonds():
        a,c=ids[b.GetBeginAtomIdx()],ids[b.GetEndAtomIdx()]
        if b.GetBondType()==Chem.BondType.DATIVE:datives.append((a,c));continue
        if b.GetBondType() not in {Chem.BondType.SINGLE,Chem.BondType.DOUBLE,Chem.BondType.TRIPLE}:raise Unsupported('Unsupported bond type')
        stereo={'status':'unknown'}
        if b.GetStereo()!=Chem.BondStereo.STEREONONE:
            stereo={'status':'specified','rdkit_stereo':str(b.GetStereo()),'reference_atoms':[ids[i] for i in b.GetStereoAtoms()]}
        bonds.append(Bond(a=a,b=c,order=int(b.GetBondTypeAsDouble()),stereo=stereo))
    for a in mol.GetAtoms():
        i=a.GetAtomMapNum();el=a.GetSymbol();h=a.GetNumExplicitHs();q=a.GetFormalCharge()
        # Mapped bracket hydrogens are explicit counts. Never fill a metal or a
        # low-coordinate atom using RDKit's default implicit-H valence.
        b=h+sum(x.order for x in bonds if i in {x.a,x.b});v=VALENCE_ELECTRONS.get(el)
        n=v-q-b if v is not None else None
        if n is not None and n<0:raise FormalError(f'Negative nonbonding count at {i}')
        r=(n%2 if n is not None else None)
        if (el in {'Cr','Ru'} and n!=0) or (el in {'Pd','Hg','Zn'} and n not in {0,8,10,12}):r=None
        if q==0 and ((el=='C' and (b,n) in {(2,2),(1,3)}) or (el=='N' and (b,n)==(1,4))):r=None
        stereo={'status':'unknown'}
        if a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED:
            stereo={'status':'specified','rdkit_tag':str(a.GetChiralTag()),'neighbor_order':[x.GetAtomMapNum() for x in a.GetNeighbors()]}
        atoms.append(Atom(atom_id=i,element=el,isotope=a.GetIsotope(),hydrogens=h,formal_charge=q if el!='*' else None,
                          nonbonding_electrons=n,radical_electrons=r,stereo=stereo,label=(labels or {}).get(i,'R' if el=='*' else None)))
    amap={a.atom_id:a for a in atoms};mids={i for i,a in amap.items() if a.element in METAL_GROUP}
    todo=set(ids)-mids;components=[]
    while todo:
        stack=[todo.pop()];seen=set(stack)
        while stack:
            i=stack.pop()
            for b in bonds:
                if i not in {b.a,b.b}:continue
                j=b.b if b.a==i else b.a
                if j in todo:todo.remove(j);seen.add(j);stack.append(j)
        components.append(seen)
    ligands=[];bindings=[];adducts=[]
    for d,a in datives:
        if a not in mids:
            if amap[a].element!='B':raise Unsupported('Dative acceptor outside metal/B scope')
            adducts.append(LewisAdduct(binding_id=f'adduct.{d}.{a}',donor_atom_id=d,acceptor_atom_id=a))
    for members in components:
        contacts=[(d,m,'dative',0) for d,m in datives if d in members and m in mids]
        contacts += [(b.a if b.a in members else b.b,b.b if b.a in members else b.a,'sigma_covalent',b.order)
                     for b in bonds if {b.a,b.b}&mids and {b.a,b.b}&members]
        if not contacts:continue
        iq=sum(amap[i].formal_charge for i in members)-sum(o for _,_,_,o in contacts)
        cls='unknown'
        if len(contacts)==1:
            d,m,mode,o=contacts[0]
            if mode=='dative' and iq==0:cls='L'
            elif mode=='sigma_covalent' and iq==-1 and o==1:cls='X'
            elif mode=='sigma_covalent' and iq==-2 and len(members)==1 and amap[d].element=='O':cls='X2'
        lid=f'ligand.{min(members)}'
        ligands.append(Ligand(ligand_id=lid,atom_ids=sorted(members),donor_atom_ids=sorted({d for d,_,_,_ in contacts}),
                              ionic_class=cls,ionic_charge=iq if cls!='unknown' else None,
                              assignment_basis='inferred finite ionic partition from mapped graph; not measured charge'))
        for k,(d,m,mode,o) in enumerate(contacts):
            bindings.append(MetalBinding(binding_id=f'bind.{lid}.{k}',metal_atom_id=m,ligand_id=lid,donor_atom_ids=[d],mode=mode))
    view=LewisView(state_id=state_id,view_id=view_id,atoms=atoms,bonds=bonds,ligands=ligands,metal_bindings=bindings,
                   metal_centers=[MetalCenter(atom_id=i) for i in sorted(mids)],lewis_adducts=adducts,
                   annotations={'raw_mapped_smiles':smiles,'rdkit_version':rdBase.rdkitVersion,
                                'rdkit_cleanup_applied':False,'implicit_hydrogens_added':False,
                                'spin_assignment':'formal finite local occupation or unknown, not a ground-state prediction'})
    view.element_declarations=[ElementDeclaration(atom_id=a.atom_id,oxidation_state=oxo_oxidation(view,a.atom_id))
                               for a in atoms if a.element in MAIN]
    state=ChemicalState(state_id=state_id,primary_view_id=view_id,total_charge=sum(a.GetFormalCharge() for a in mol.GetAtoms()),
                        multiplicity=None,electronic_mode='UNKNOWN')
    if metadata is not None:
        if metadata.get('schema_version')!='umr-smiles/0.5' or metadata.get('smiles_sha256')!=hashlib.sha256(smiles.encode()).hexdigest():
            raise FormalError('SMILES sidecar version/hash mismatch')
        s=ChemicalState.model_validate(metadata['state']);v=LewisView.model_validate(metadata['view'])
        def carrier(x):
            return ({a.atom_id:(a.element,a.isotope,a.hydrogens,a.formal_charge) for a in x.atoms},bond_orders(x),
                    sorted((b.donor_atom_ids[0],b.metal_atom_id) for b in x.metal_bindings if b.mode=='dative'),
                    sorted((b.donor_atom_id,b.acceptor_atom_id) for b in x.lewis_adducts))
        if carrier(v)!=carrier(view):raise FormalError('Sidecar contradicts mapped graph, charge, H or dative direction')
        def stereo_carrier(x):
            clone=x.model_copy(deep=True)
            for atom in clone.atoms:atom.radical_electrons=0
            return Chem.MolToSmiles(to_rdkit(clone),canonical=True,isomericSmiles=True,kekuleSmiles=True)
        if stereo_carrier(v)!=stereo_carrier(view):raise FormalError('Sidecar contradicts mapped stereochemistry')
        if v.state_id!=s.state_id or v.view_id!=metadata.get('represented_view_id') or v.electronic_state!=s.electronic_state:raise FormalError('Sidecar state/view identity mismatch')
        return s,v
    return state,view


def to_rdkit(view):
    carrier=view.model_copy(deep=True);unknown=[]
    for a in carrier.atoms:
        if a.radical_electrons is None:unknown.append(a.atom_id);a.radical_electrons=0
    rw=Chem.RWMol(shared_rdkit(carrier,_sanitize_ops=FLAGS));ids={a.GetAtomMapNum():a.GetIdx() for a in rw.GetAtoms()}
    pairs=[(b.donor_atom_ids[0],b.metal_atom_id) for b in view.metal_bindings if b.mode=='dative']
    pairs += [(b.donor_atom_id,b.acceptor_atom_id) for b in view.lewis_adducts]
    for d,a in pairs:rw.AddBond(ids[d],ids[a],Chem.BondType.DATIVE)
    mol=rw.GetMol();mol.UpdatePropertyCache(strict=False)
    # This is carrier metadata. RDKit draws atomNote as text at the atom,
    # which would suggest a localized spin assignment for an unknown state.
    for i in unknown:mol.GetAtomWithIdx(ids[i]).SetBoolProp('_UMRSpinOccupationUnknown',True)
    return mol


def export_mapped_smiles(view):return Chem.MolToSmiles(to_rdkit(view),canonical=False,isomericSmiles=True,kekuleSmiles=True)


def sidecar(state,view,smiles):
    return {'schema_version':'umr-smiles/0.5','smiles_sha256':hashlib.sha256(smiles.encode()).hexdigest(),
            'state':state.model_dump(mode='json'),'view':view.model_dump(mode='json'),'represented_view_id':view.view_id,
            'lossless_without_sidecar':False,'bare_smiles_electronic_state':'ambiguous',
            'representation_policy':'No RDKit cleanup or implicit-H completion of authoritative Lewis graph'}
