"""Stable-map local parity is authoritative; CIP is recomputed, never flipped by name."""
from rdkit import Chem, rdBase
from rdkit.Chem import rdCIPLabeler
from ..chemistry import FormalError, Unsupported, atoms_by_id
from ..v05.smiles import to_rdkit

RULE_VERSION='umr-local-stereo/0.6'


def analyze(view):
    try:
        mol=to_rdkit(view)
        potential=Chem.FindPotentialStereo(mol,cleanIt=False,flagPossible=True)
        rdCIPLabeler.AssignCIPLabels(mol,maxRecursiveIterations=1000000)
    except (ValueError,RuntimeError,Unsupported) as exc:
        return {'status':'unsupported','reason':str(exc),'rdkit_version':rdBase.rdkitVersion,
                'local_stereo_retained':True,'absolute_configuration_inferred':False}
    atom_potential={s.centeredOn for s in potential if s.type==Chem.StereoType.Atom_Tetrahedral}
    atoms=[]
    for a in mol.GetAtoms():
        original=atoms_by_id(view)[a.GetAtomMapNum()]
        specified=original.stereo.get('status')=='specified'
        atoms.append({'atom_id':a.GetAtomMapNum(),
            'status':'specified' if specified and a.GetIdx() in atom_potential else 'unknown' if a.GetIdx() in atom_potential else 'not_stereocenter',
            'local_relation':original.stereo,
            'cip':a.GetProp('_CIPCode') if a.HasProp('_CIPCode') and specified else None})
    bonds=[{'atom_ids':[b.GetBeginAtom().GetAtomMapNum(),b.GetEndAtom().GetAtomMapNum()],
            'rdkit_stereo':str(b.GetStereo()),'cip':b.GetProp('_CIPCode') if b.HasProp('_CIPCode') else None}
           for b in mol.GetBonds() if b.GetBondType()==Chem.BondType.DOUBLE]
    return {'status':'complete','atoms':atoms,'double_bonds':bonds,'sample_composition':view.annotations.get('sample_stereochemistry','unknown'),
            'rdkit_version':rdBase.rdkitVersion,'rule_version':RULE_VERSION,
            'method':'FindPotentialStereo + rdCIPLabeler.AssignCIPLabels; max iterations 1000000',
            'absolute_configuration_inferred':False,'unsupported_scope':['axial chirality','metal coordination stereochemistry']}


def apply_event(source,target,event):
    before=atoms_by_id(source);after=atoms_by_id(target);effects=[]
    if event.event_type=='SN2':
        c=event.bindings['electrophile'];x=event.bindings['leaving_group'];nu=event.bindings['nucleophile']
        old=before[c].stereo
        if old.get('status')=='specified' and old.get('rdkit_tag') in {'CHI_TETRAHEDRAL_CW','CHI_TETRAHEDRAL_CCW'}:
            refs=old.get('neighbor_order',[])
            if x not in refs:raise FormalError('SN2 leaving atom absent from local stereo reference')
            after[c].stereo={**old,'neighbor_order':[nu if i==x else i for i in refs],
                'rdkit_tag':'CHI_TETRAHEDRAL_CCW' if old['rdkit_tag']=='CHI_TETRAHEDRAL_CW' else 'CHI_TETRAHEDRAL_CW',
                'basis':'template-relative backside substitution; CIP independently recalculated'}
            after[c].stereo.pop('cip',None)
            effects.append({'atom_id':c,'effect':'inversion','basis':'local neighbor replacement and odd parity; not R/S toggle'})
        else:
            effects.append({'atom_id':c,'effect':'unknown','basis':'SN2 relative inversion constraint; reactant absolute configuration unspecified'})
    actual=analyze(target)
    prior=analyze(source)
    if actual['status']=='complete' and prior['status']=='complete':
        previous={r['atom_id']:r for r in prior['atoms']}
        for row in actual['atoms']:
            i=row['atom_id'];old=previous[i]
            if any(e['atom_id']==i for e in effects):continue
            effect=None
            if old['status']!='not_stereocenter' and row['status']=='not_stereocenter':effect='destroyed'
            elif old['status']=='not_stereocenter' and row['status']=='unknown':effect='created_unspecified'
            elif old['status']=='specified' and row['status']=='specified' and before[i].stereo==after[i].stereo:effect='retention'
            elif old['status']!='not_stereocenter' or row['status']!='not_stereocenter':effect='unknown'
            if effect:effects.append({'atom_id':i,'effect':effect,'basis':'mapped local relationship and potential-center reassessment; no sample configuration invented'})
    # Declared effects are constraints to check, not permissions to invent geometry.
    for assertion in event.stereo_effects:
        if assertion.atom_id not in before:raise FormalError('Stereo effect references absent atom')
        match=next((e for e in effects if e['atom_id']==assertion.atom_id),None)
        if assertion.effect!='unknown' and (match is None or match['effect']!=assertion.effect):
            raise Unsupported('Declared stereo effect has no implemented event-relative proof')
    return {'effects':effects,'product_analysis':actual,'rule_version':RULE_VERSION}
