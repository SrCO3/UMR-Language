from ..v05.smiles import import_mapped_smiles as import05, export_mapped_smiles, to_rdkit, sidecar as sidecar05
from ..chemistry import FormalError

def import_mapped_smiles(smiles,**kwargs):
    metadata=kwargs.get('metadata')
    if metadata:
        if metadata.get('schema_version')!='umr-smiles/0.6':raise FormalError('Expected UMR 0.6 SMILES sidecar')
        kwargs['metadata']={**metadata,'schema_version':'umr-smiles/0.5'}
    state,view=import05(smiles,_stereo_perception=True,_preserve_source_hydrogens=True,**kwargs)
    if not metadata:
        view.annotations['input_hydrogen_policy']='source bracket counts preserved before Kekulization; no implicit-H completion'
        from .stereo import analyze
        view.annotations['stereochemistry_analysis']=analyze(view)
    return state,view

def sidecar(state,view,smiles):
    return {**sidecar05(state,view,smiles),'schema_version':'umr-smiles/0.6'}
