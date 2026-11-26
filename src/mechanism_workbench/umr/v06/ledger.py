from ..chemistry import atoms_by_id,bond_valences,FormalError
from ..v05.elements import local_profile as profile05,charge
from ..v05.ledger import (validate_view as validate05,validate_ledger as ledger05,compare_views,signature)


def local_profile(view,atom):
    b=bond_valences(view)[atom.atom_id]
    if atom.element!='Si' or b not in {5,6}:return profile05(view,atom)
    # Only explicitly mapped fluoride coordination; not general hypervalence.
    neighbors=[(atoms_by_id(view)[x.b if x.a==atom.atom_id else x.a],x.order) for x in view.bonds if atom.atom_id in {x.a,x.b}]
    if atom.hydrogens or len(neighbors)!=b or any(a.element!='F' or o!=1 for a,o in neighbors):
        raise FormalError('Expanded Si profile is restricted to five/six monovalent fluoride ligands')
    if (atom.nonbonding_electrons,atom.radical_electrons,atom.formal_charge)!=(0,0,4-b):
        raise FormalError('Fluorosilicate charge/occupation mismatch')
    return {'rule':f'Si_fluoride_{b}_coordinate/0.6','atom_id':atom.atom_id,'element':'Si','valence_category':f'Si_Lewis_{2*b}e',
        'bond_order_sum':b,'coordination_number':b,'N':0,'r':0,'formal_charge':4-b,'shell_electrons':2*b,
        'oxidation_state':None,'chemical_feasibility':'unknown','source':'https://doi.org/10.1021/ic00178a016' if b==5 else 'https://doi.org/10.1039/D2RA04270C'}


def validate_ledger(view):return ledger05(view,_profile=local_profile)
def validate_view(view,state=None):return validate05(view,state,_profile=local_profile)
