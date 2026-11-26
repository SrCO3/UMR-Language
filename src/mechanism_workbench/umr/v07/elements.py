"""Finite element profiles. Shared chemistry roles do not equate electron states."""
from functools import lru_cache
import json
from pathlib import Path
from ..chemistry import FormalError, Unsupported, atoms_by_id, bond_valences
from ..v05.elements import (METAL_GROUP as GROUP05,METAL_OXIDATION as OS05,
    VALENCE_ELECTRONS as V05,MAIN as MAIN05,local_profile as profile05)
from ..v05.coordination import checks as checks05

METAL_GROUP={**GROUP05,'Fe':8,'Cu':11}
METAL_OXIDATION={**OS05,'Fe':{2,3},'Cu':{0,1,2}}
VALENCE={**V05,**METAL_GROUP,'Sn':4}
MAIN=MAIN05|{'Sn'}
TRANSITION=frozenset({'Cr','Ru','Pd','Zn','Hg','Fe','Cu'})


def charge(atom,b):
    if atom.element not in VALENCE or atom.nonbonding_electrons is None:raise Unsupported('Unknown local electron ownership')
    return VALENCE[atom.element]-atom.nonbonding_electrons-b


def coordination(view):
    return checks05(view,_metal_group=METAL_GROUP,_metal_oxidation=METAL_OXIDATION,_transition_metals=TRANSITION)


def local_profile(view,atom):
    b=bond_valences(view)[atom.atom_id];n=atom.nonbonding_electrons;r=atom.radical_electrons;q=atom.formal_charge
    if q!=charge(atom,b):raise FormalError(f'Charge ledger mismatch at {atom.atom_id}')
    if n<0 or (r is not None and (n<r or (n-r)%2)):raise FormalError('N includes r; incompatible occupation')
    el=atom.element;ns=[(atoms_by_id(view)[e.b if e.a==atom.atom_id else e.a],e.order) for e in view.bonds if atom.atom_id in {e.a,e.b}]
    if el in {'Fe','Cu'}:
        if atom.hydrogens or len(ns)>6 or n>VALENCE[el]:raise Unsupported('Finite single-center Fe/Cu ownership scope exceeded')
        name=el+'_ionic_ownership'
    elif el=='Sn':
        if (b,n,r,q) not in {(4,0,0,0),(3,1,1,0)} or any(o!=1 for _,o in ns):
            raise Unsupported('Sn scope: neutral four-coordinate stannane or neutral R3Sn radical')
        name='Sn_tetravalent' if b==4 else 'Sn_stannyl_radical'
    elif el=='Si' and b in {5,6}:
        if n!=0 or r!=0 or q!=4-b or atom.hydrogens or len(ns)!=b or any(o!=1 for _,o in ns):
            raise FormalError('Hypercoordinate Si requires explicit single ligands and consistent charge')
        if b==6 and any(a.element!='F' for a,_ in ns):raise Unsupported('Six-coordinate Si scope is fluoride only')
        if any(a.element not in {'C','O','F','Cl'} for a,_ in ns):raise Unsupported('Five-coordinate Si ligand class outside scope')
        name=f'Si_{b}_coordinate'
    elif el=='Se' and r:
        if (b,n,r,q)!=(1,5,1,0):raise Unsupported('Se radical scope is neutral monosubstituted selenyl')
        name='Se_selenyl_radical'
    else:return profile05(view,atom)
    return {'rule':name+'/0.7','atom_id':atom.atom_id,'element':el,'bond_order_sum':b,
            'coordination_number':len(ns)+atom.hydrogens,'N':n,'r':r,'formal_charge':q,
            'shell_electrons':n+2*b,'chemical_feasibility':'unknown'}


@lru_cache(maxsize=1)
def electronegativities():
    bundled = Path(__file__).resolve().parents[2] / 'resources' / 'allen-electronegativity.json'
    return json.loads(bundled.read_text('utf8'))['elements']


def oxidation(view,atom_id):
    atom=atoms_by_id(view)[atom_id];table=electronegativities()
    if atom.formal_charge is None or atom.element not in table:return None
    if any(atom_id in b.donor_atom_ids or atom_id==b.metal_atom_id for b in view.metal_bindings if b.mode!='sigma_covalent'):return None
    value=atom.formal_charge
    pairs=[(atoms_by_id(view)[b.b if b.a==atom_id else b.a].element,b.order) for b in view.bonds if atom_id in {b.a,b.b}]
    if atom.hydrogens:pairs.append(('H',atom.hydrogens))
    for other,n in pairs:
        if other==atom.element:continue
        if other not in table or table[other]['value'] is None or table[atom.element]['value'] is None or table[other]['value']==table[atom.element]['value']:return None
        value+=n*(1 if table[other]['value']>table[atom.element]['value'] else -1)
    return value
