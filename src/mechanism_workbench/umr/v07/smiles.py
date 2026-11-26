import hashlib
from ..v05.smiles import import_mapped_smiles as import05,to_rdkit,export_mapped_smiles
from ..chemistry import FormalError,bond_orders
from .ast import ChemicalState,LewisView
from .elements import VALENCE,METAL_GROUP,MAIN


def import_mapped_smiles(smiles,*,metadata=None,**kwargs):
    state,old=import05(smiles,_stereo_perception=True,_preserve_source_hydrogens=True,
                      _valence=VALENCE,_metals=METAL_GROUP,_main=MAIN,**kwargs)
    raw=old.model_dump(mode='json');raw['element_rule_version']='umr-elements/0.7';view=LewisView.model_validate(raw)
    for a in view.atoms:
        if a.element=='Fe':a.radical_electrons=None
        if a.element=='Cu':a.radical_electrons=None
    view.annotations['input_hydrogen_policy']='pre-Kekulization bracket H retained'
    if metadata:
        if metadata.get('schema_version')!='umr-smiles/0.7' or metadata.get('smiles_sha256')!=hashlib.sha256(smiles.encode()).hexdigest():
            raise FormalError('SMILES sidecar version/hash mismatch')
        s=ChemicalState.model_validate(metadata['state']);v=LewisView.model_validate(metadata['view'])
        def carrier(x):return ({a.atom_id:(a.element,a.isotope,a.hydrogens,a.formal_charge) for a in x.atoms},bond_orders(x),
            sorted((tuple(b.donor_atom_ids),b.metal_atom_id,b.mode) for b in x.metal_bindings if b.mode=='dative'),
            sorted((b.donor_atom_id,b.acceptor_atom_id) for b in x.lewis_adducts))
        if carrier(v)!=carrier(view):raise FormalError('Sidecar contradicts mapped graph, charge, H or dative direction')
        from rdkit import Chem
        def stereo(x):
            clone=x.model_copy(deep=True)
            for a in clone.atoms:a.radical_electrons=0
            return Chem.MolToSmiles(to_rdkit(clone),canonical=True,isomericSmiles=True,kekuleSmiles=True)
        if stereo(v)!=stereo(view):raise FormalError('Sidecar contradicts mapped stereochemistry')
        if s.state_id!=v.state_id or s.electronic_state!=v.electronic_state or s.primary_view_id!=metadata.get('primary_view_id',s.primary_view_id) or v.view_id!=metadata['represented_view_id']:
            raise FormalError('Sidecar references disagree')
        return s,v
    return state,view


def sidecar(state,view,smiles):
    return {'schema_version':'umr-smiles/0.7','smiles_sha256':hashlib.sha256(smiles.encode()).hexdigest(),
            'state':state.model_dump(mode='json'),'view':view.model_dump(mode='json'),
            'represented_view_id':view.view_id,'primary_view_id':state.primary_view_id,'lossless_without_sidecar':False}
