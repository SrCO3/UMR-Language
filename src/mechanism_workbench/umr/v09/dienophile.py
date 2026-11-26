"""Two explicitly activated heterodienophile leaves of the existing ordered [4+2] compiler."""
from collections import Counter
import json
from rdkit import Chem
from ..chemistry import FormalError,Unsupported,atoms_by_id
from ..v08.pericyclic import compile_profile
from .smiles import to_rdkit

VERSION='0.9-dienophile/1'
PROFILES={
    'AZODICARBOXYLATE_DA':('42',('CCCCNN',),'neutral azodicarboxylate and carbon 1,3-diene'),
    'ARYLNITROSO_DA':('42',('CCCCNO',),'neutral arylnitroso and carbon 1,3-diene'),
}
ACTIVATED={
    'AZODICARBOXYLATE_DA':'[N;X2;+0;$([N]-C(=O)-O-[#6])]=[N;X2;+0;$([N]-C(=O)-O-[#6])]',
    'ARYLNITROSO_DA':'[N;X2;+0;$([N]-[c])]=[O;X1;+0]',
}


def execute(source,event):
    from .execution import transact
    from .ledger import validate_view
    from .ast import RULE_VERSION
    if (event.template_version!=VERSION or event.event_type!='PERICYCLIC' or event.direction!='forward'
        or event.execution_level!='electron_actions' or event.mechanism_assumption!='concerted' or not event.evidence):
        raise Unsupported('Heterodienophile requires an evidence-bound forward concerted event')
    if (event.electronic_outcome or event.electron_transfers or event.stereo_effects or event.stereo_constraints
        or any((event.general_event,event.rule_template,event.element_template,event.reactive_template,event.metal_template))):
        raise Unsupported('Heterodienophile leaf does not authorize additional transformations or stereo assignments')
    validate_view(source)
    ids=event.bindings.get('ordered_atoms',[])
    if set(event.bindings)!={'ordered_atoms'} or len(ids)!=6 or len(set(ids))!=6 or not set(ids)<=atoms_by_id(source).keys():
        raise FormalError('Six distinct mapped atoms required')
    active=set(ids);atoms=atoms_by_id(source)
    if (any(atoms[i].formal_charge or atoms[i].radical_electrons!=0 for i in active)
        or any(active&set(b.donor_atom_ids) for b in source.metal_bindings)
        or any(active&{b.donor_atom_id,b.acceptor_atom_id} for b in source.lewis_adducts)
        or any(active&set(b.atom_ids) and b.direction!='inside' for b in source.material_boundaries)
        or any(c.atom_id in active for c in source.reactive_centers)):
        raise Unsupported('Registered heterodienophile requires free neutral closed-shell participants')
    mol=to_rdkit(source);Chem.SanitizeMol(mol)
    match=Chem.MolFromSmarts(ACTIVATED[event.template_id])
    if not any([mol.GetAtomWithIdx(i).GetAtomMapNum() for i in hit]==ids[4:]
               for hit in mol.GetSubstructMatches(match,uniquify=False,maxMatches=2*mol.GetNumBonds()+1)):
        raise Unsupported('Bound dienophile lacks the registered activation/aryl environment')
    actions,report=compile_profile(source,event,profiles=PROFILES,template_version=VERSION,
                                   rule_version=VERSION,neutral_linkers=set())
    signature=lambda x:json.dumps(x.model_dump(mode='json'),sort_keys=True)
    if event.actions is not None and Counter(map(signature,event.actions))!=Counter(map(signature,actions)):
        raise FormalError('Submitted action group differs from the registered [4+2] cycle')
    target=transact(source,actions,event.product_state_id,event.product_view_id)
    event.actions=actions
    report.update(event_id=event.event_id,event_type='PERICYCLIC',rule_version=RULE_VERSION,profile_rule_version=VERSION,
                  atomic_commit=True,formal_status='passed',execution_mode='electron_actions',participating_electrons=6,
                  compiled_actions=[x.model_dump(mode='json') for x in actions],chemical_feasibility='unknown',
                  source=[e.model_dump(mode='json') if hasattr(e,'model_dump') else e for e in event.evidence])
    return target,report
