import hashlib
from ..v05.smiles import import_mapped_smiles as shared, export_mapped_smiles, to_rdkit
from ..chemistry import FormalError, bond_orders
from .ast import ChemicalState,LewisView
from .elements import VALENCE,METAL_GROUP,MAIN,NEW_METALS


def import_mapped_smiles(smiles,*,metadata=None,**kwargs):
    state,old=shared(smiles,_stereo_perception=True,_preserve_source_hydrogens=True,
                     _valence=VALENCE,_metals=METAL_GROUP,_main=MAIN,**kwargs)
    raw=old.model_dump(mode='json');raw['element_rule_version']='umr-elements/0.8'
    view=LewisView.model_validate(raw)
    for a in view.atoms:
        if a.element in {'Ce','Sm','Ti','Mn','Os','Rh','Co','Pt','Fe','Cu'}:a.radical_electrons=None
    if metadata:
        if metadata.get('schema_version')!='umr-smiles/0.8' or metadata.get('smiles_sha256')!=hashlib.sha256(smiles.encode()).hexdigest():
            raise FormalError('SMILES sidecar version/hash mismatch')
        s=ChemicalState.model_validate(metadata['state']);v=LewisView.model_validate(metadata['view'])
        carrier=lambda x:({a.atom_id:(a.element,a.isotope,a.hydrogens,a.formal_charge) for a in x.atoms},bond_orders(x),
                         sorted((tuple(b.donor_atom_ids),b.metal_atom_id,b.mode) for b in x.metal_bindings if b.mode=='dative'),
                         sorted((b.donor_atom_id,b.acceptor_atom_id) for b in x.lewis_adducts))
        if carrier(v)!=carrier(view):raise FormalError('Sidecar contradicts graph/material/charge')
        from rdkit import Chem
        def stereo(x):
            c=x.model_copy(deep=True)
            for a in c.atoms:a.radical_electrons=0
            return Chem.MolToSmiles(to_rdkit(c),canonical=True,isomericSmiles=True,kekuleSmiles=True)
        if stereo(v)!=stereo(view) or s.state_id!=v.state_id or s.electronic_state!=v.electronic_state or metadata.get('represented_view_id')!=v.view_id:
            raise FormalError('Sidecar contradicts stereo or state/view identities')
        return s,v
    return state,view


def sidecar(state,view,smiles):
    return {'schema_version':'umr-smiles/0.8','smiles_sha256':hashlib.sha256(smiles.encode()).hexdigest(),
            'state':state.model_dump(mode='json'),'view':view.model_dump(mode='json'),
            'represented_view_id':view.view_id,'lossless_without_sidecar':False}
