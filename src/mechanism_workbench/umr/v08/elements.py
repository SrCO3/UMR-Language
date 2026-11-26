"""Element identity is universal; occupation and reaction capabilities are finite."""
from ..v07.elements import (VALENCE as OLD_V, MAIN, METAL_GROUP as OLD_M, TRANSITION as OLD_T,
                            METAL_OXIDATION as OLD_OS, local_profile as old_profile, oxidation)
from ..v05.coordination import checks
from ..chemistry import FormalError, Unsupported, atoms_by_id, bond_valences

NEW_METALS={'Ag':11,'Ti':4,'Mn':7,'Os':8,'Rh':9,'Co':9,'Pt':10,'Ba':2,'La':3,'Ce':4,'Sm':8}
METAL_GROUP={**OLD_M,**NEW_METALS};VALENCE={**OLD_V,**NEW_METALS}
OS={**OLD_OS,'Ag':{1},'Ti':{4},'Mn':{2,4,7},'Os':{4,8},'Rh':{1,3},'Co':{2,3},'Pt':{0,2,4},
    'Ba':{2},'La':{3},'Ce':{2,3,4},'Sm':{2,3}}
F_BLOCK={'La','Ce','Sm'}
TRANSITION=OLD_T|{'Ag','Ti','Mn','Os','Rh','Co','Pt'}


def charge(atom,b):
    if atom.element not in VALENCE or atom.nonbonding_electrons is None:raise Unsupported('Unregistered valence ownership')
    return VALENCE[atom.element]-atom.nonbonding_electrons-b


def profile(view,atom):
    if atom.element not in NEW_METALS:return old_profile(view,atom)
    b=bond_valences(view)[atom.atom_id];n=atom.nonbonding_electrons;r=atom.radical_electrons;q=atom.formal_charge
    if q!=charge(atom,b) or n<0 or r is not None and (r>n or (n-r)%2):raise FormalError('Extended element charge/occupation mismatch')
    if atom.hydrogens:raise Unsupported('New metal hydrides require explicit atom IDs and a dedicated profile')
    if atom.element in F_BLOCK and b:
        raise Unsupported('Initial lanthanide profile uses ionic center plus explicit donor-owned coordination, not covalent f-electron partition')
    if atom.element=='Ce' and q==2:
        evidence=view.annotations.get('specific_CeII_complex_evidence')
        if not isinstance(evidence,dict) or not all(evidence.get(k) for k in ('source','complex_identity','ligand_redox_assignment')):
            raise Unsupported('Ce(II) needs a cited specific complex and explicit ligand-redox assignment; no generic reductant')
    if n>VALENCE[atom.element]:raise Unsupported('Extended negative metal occupation not registered')
    return {'atom_id':atom.atom_id,'element':atom.element,'bond_order_sum':b,'N':n,'r':r,'formal_charge':q,
            'rule':'finite_ionic_ownership/0.8','shell_electrons':n+2*b,'chemical_feasibility':'unknown',
            'f_shell_is_not_d_count':atom.element in F_BLOCK}


def coordination(view):
    atoms=atoms_by_id(view)
    for m in view.metal_centers:
        if atoms[m.atom_id].element in F_BLOCK and (m.d_electrons is not None or m.valence_electrons is not None):
            raise Unsupported('Do not apply d-count/18-electron declarations to lanthanide adapter')
    reserved,rows=checks(view,_metal_group=METAL_GROUP,_metal_oxidation=OS,_transition_metals=TRANSITION)
    for row in rows:
        if row['element'] in NEW_METALS:
            maximum=12 if row['element'] in F_BLOCK else 10 if row['element']=='Ba' else 6
            if row['coordination_number']>maximum:raise Unsupported('Coordination capacity exceeds finite profile')
        if row['element'] in F_BLOCK:
            row.update(d_electrons=None,valence_electrons=None,eighteen_electron_rule='not_applicable',
                       local_nonbonding_electrons_including_f=atoms[row['metal_atom_id']].nonbonding_electrons,
                       counting_method='explicit ionic electron ownership; coupling unknown')
    return reserved,rows
