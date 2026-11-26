"""Mapped carrier plus mandatory metadata for coordination/electronic meaning."""
import hashlib
from rdkit import Chem
from ..v05.smiles import import_mapped_smiles as shared, to_rdkit as rdkit05
from ..chemistry import FormalError,bond_orders
from .ast import ChemicalState
from .elements import current_elements

def to_rdkit(view):
    # Old RDKit adapter accepted only the first atom of a donor list.
    c=view.model_copy(deep=True);expanded=[]
    for b in c.metal_bindings:
        if b.mode=='dative':
            for i in b.donor_atom_ids:
                x=b.model_copy(deep=True);x.donor_atom_ids=[i];expanded.append(x)
        else:expanded.append(b)
    c.metal_bindings=expanded
    return rdkit05(c)

def export_mapped_smiles(view):
    return Chem.MolToSmiles(to_rdkit(view),canonical=False,isomericSmiles=True,kekuleSmiles=True)

def sidecar(state,view,smiles):
    return {'schema_version':'umr-smiles/0.9','smiles_sha256':hashlib.sha256(smiles.encode()).hexdigest(),
      'state':state.model_dump(mode='json'),'view':view.model_dump(mode='json'),'represented_view_id':view.view_id,
      'lossless_without_sidecar':False,'eta_and_spin_in_sidecar':True,'bare_smiles_is_not_execution_authority':True}

def import_mapped_smiles(smiles,*,metadata=None,**kwargs):
    s,old=shared(smiles,_stereo_perception=True,_preserve_source_hydrogens=True,
       _valence=current_elements().valences,_metals=current_elements().metals,_main=current_elements().main,**kwargs)
    rules=current_elements()
    raw=old.model_dump(mode='json');raw['element_rule_version']=rules.version;v=rules.view_type.model_validate(raw)
    for a in v.atoms:
        if a.element in current_elements().metals:
            # For the newly registered Ca2+ valence shell N=0 fixes r=0.
            # Do not turn unknown occupations of legacy metals into zero.
            a.radical_electrons=0 if a.element=='Ca' and a.nonbonding_electrons==0 else None
    # Chelation is a topology observation, not a claim of complex stability.
    for l in v.ligands:
        bs=[b for b in v.metal_bindings if b.ligand_id==l.ligand_id]
        if len(bs)==2 and len({b.metal_atom_id for b in bs})==1 and all(b.mode=='dative' for b in bs) and sum(a.formal_charge for a in v.atoms if a.atom_id in l.atom_ids)==0:
            l.ionic_class='L2';l.ionic_charge=0
    if metadata is None:return s,v
    if metadata.get('schema_version')!='umr-smiles/0.9' or metadata.get('smiles_sha256')!=hashlib.sha256(smiles.encode()).hexdigest():raise FormalError('SMILES sidecar version/hash mismatch')
    target=rules.view_type.model_validate(metadata['view']);state=ChemicalState.model_validate(metadata['state'])
    def carrier(x):
        return ({a.atom_id:(a.element,a.isotope,a.hydrogens,a.formal_charge) for a in x.atoms},bond_orders(x),
          sorted((i,b.metal_atom_id) for b in x.metal_bindings if b.mode=='dative' for i in b.donor_atom_ids),
          sorted((b.donor_atom_id,b.acceptor_atom_id) for b in x.lewis_adducts))
    if carrier(v)!=carrier(target):raise FormalError('Sidecar contradicts mapped graph, H, charge or donor direction')
    def stereo(x):
        c=x.model_copy(deep=True)
        for a in c.atoms:a.radical_electrons=0
        return Chem.MolToSmiles(to_rdkit(c),canonical=True,isomericSmiles=True,kekuleSmiles=True)
    if stereo(v)!=stereo(target):raise FormalError('Sidecar contradicts stereochemistry')
    if state.state_id!=target.state_id or state.electronic_state!=target.electronic_state or metadata.get('represented_view_id')!=target.view_id:raise FormalError('Sidecar state/view mismatch')
    return state,target
