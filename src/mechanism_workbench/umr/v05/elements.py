"""Finite Lewis profiles; oxidation state is a separate ionic approximation."""
from ..chemistry import VALENCE, FormalError, Unsupported, atoms_by_id, bond_valences

RULE_VERSION='umr-elements/0.5'
METAL_GROUP={'Li':1,'Na':1,'K':1,'Mg':2,'Al':3,'Cr':6,'Ru':8,'Pd':10,'Zn':12,'Hg':12}
VALENCE_ELECTRONS={**VALENCE,**METAL_GROUP,'Si':4,'Se':6}
MAIN={'B','S','Se','P','Cl','Br','I','Si'}
METAL_OXIDATION={'Li':{0,1},'Na':{0,1},'K':{0,1},'Mg':{0,2},'Al':{0,3},
                 'Cr':{0,2,3,6},'Ru':{0,2,3,4,6,8},'Pd':{0,2,4},'Zn':{0,2},'Hg':{0,1,2}}
SOURCES={'rdkit':'https://www.rdkit.org/docs/RDKit_Book.html#valence-calculation-and-allowed-valences',
         'ionic_counting':'https://www.ilpi.com/organomet/electroncount.html',
         'oxidation_state':'https://doi.org/10.1515/pac-2015-1204'}


def charge(atom,b):
    if atom.element not in VALENCE_ELECTRONS or atom.nonbonding_electrons is None:
        raise Unsupported(f'Unknown electron ownership at {atom.atom_id}:{atom.element}')
    return VALENCE_ELECTRONS[atom.element]-atom.nonbonding_electrons-b


def local_profile(view,atom):
    b=bond_valences(view)[atom.atom_id];n=atom.nonbonding_electrons;r=atom.radical_electrons;q=atom.formal_charge
    if n is None or q is None:raise Unsupported(f'Unknown occupancy/charge at {atom.atom_id}')
    if q!=charge(atom,b):raise FormalError(f'Charge ledger mismatch at atom {atom.atom_id}')
    if r is not None and (n<r or (n-r)%2):raise FormalError(f'N includes r; N-r must be nonnegative and even at {atom.atom_id}')
    neighbors=[(atoms_by_id(view)[e.b if e.a==atom.atom_id else e.a],e.order) for e in view.bonds if atom.atom_id in {e.a,e.b}]
    cn=len(neighbors)+atom.hydrogens;shell=n+2*b;el=atom.element
    if el in METAL_GROUP:
        if atom.hydrogens:raise Unsupported('Metal hydrogens must be individually numbered atoms')
        if n>METAL_GROUP[el]+2 or cn>8:raise Unsupported('Metal ownership/coordination outside finite element scope')
        if el in {'Li','Na','K'} and (b,n,r,q) not in {(0,1,1,0),(0,0,0,1),(1,0,0,0)}:
            raise Unsupported('Registered alkali scopes: atomic metal, spectator cation, one sigma X ligand')
        name=f'{el}-ionic-ownership';category='metal_coordination'
    else:
        if r is None:raise Unsupported('Unknown main-group occupation')
        if r:raise Unsupported('Expanded-valence radical profile not registered; ordinary radicals use legacy rules')
        if el in {'S','Se'}:
            valid=shell in {8,10,12} and b<=6 and n<=8 and all(o<=2 for _,o in neighbors)
            category=f'{el}_Lewis_{shell}e'
        elif el=='P':
            valid=shell in {8,10} and b<=5 and n<=8 and all(o<=2 for _,o in neighbors)
            category=f'P_Lewis_{shell}e'
        elif el in {'Cl','Br','I'}:
            valid=shell in {8,10,12,14} and b<=7 and n<=8 and all(o<=2 for _,o in neighbors)
            if shell>8 and el in {'Cl','Br'}:
                valid &= all(a.element in {'O','F'} for a,_ in neighbors) and not atom.hydrogens
            category=f'{el}_Lewis_{shell}e'
        elif el=='B':
            valid=(b,n,q,cn) in {(3,0,0,3),(4,0,-1,4)}
            category='B_three_coordinate_deficient' if b==3 else 'B_four_coordinate_borate'
        elif el=='Si':
            valid=(b,n,q)==(4,0,0);category='Si_tetravalent'
        else:raise Unsupported('No new element rule')
        if not valid:raise FormalError(f'Unregistered {el} occupation: b={b}, N={n}, r={r}, q={q}, CN={cn}, shell={shell}')
        name=category
    return {'rule':name+'/0.5','atom_id':atom.atom_id,'element':el,'valence_category':category,
            'bond_order_sum':b,'coordination_number':cn,'N':n,'r':r,'formal_charge':q,
            'shell_electrons':shell,'oxidation_state':None,'chemical_feasibility':'unknown'}


def oxo_oxidation(view,atom_id):
    """O/F-only ionic approximation; ambiguous hydrides are not guessed."""
    a=atoms_by_id(view)[atom_id];neighbors=[]
    if a.hydrogens:return None
    for b in view.bonds:
        if atom_id in {b.a,b.b}:
            other=atoms_by_id(view)[b.b if b.a==atom_id else b.a]
            if other.element not in {'O','F'}:return None
            neighbors.append((other,b.order))
    if not neighbors:return a.formal_charge
    return a.formal_charge+sum(o if other.element!='H' else -o for other,o in neighbors)-a.hydrogens
