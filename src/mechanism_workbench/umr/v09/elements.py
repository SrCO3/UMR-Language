from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from typing import Callable, Mapping

from ..v08.elements import VALENCE as OLD, METAL_GROUP as METALS08,OS as OS08,MAIN,profile as profile08
from ..chemistry import atoms_by_id,bond_valences,FormalError,Unsupported
from .ast import LewisView
METALS={**METALS08,'Ni':10,'Mo':6,'W':6,'Cs':1,'Pb':4,'Ca':2,'Zr':4}
VALENCE={**OLD,**METALS}
OXIDATION={**OS08,'Cu':{0,1,2,3},'Ni':{0,1,2,3},'Cr':{0,2,3,4,6},'Co':{0,1,2,3},
 'Rh':{0,1,2,3},'Ti':{0,2,3,4},'Mn':{2,3,4,5,6,7},'Os':{4,6,8},'Ru':{0,2,3,4,5,6,7,8},'W':{6},
 'Mo':{0,2,4,6},'Zn':{0,2},'Mg':{2},'Ca':{2},'Li':{0,1},'Cs':{1},'Sm':{0,2,3},'Ce':{2,3,4},'Pb':{2,4},'Zr':{4}}
def legacy_metal_counts(element, oxidation):
    f_block = element in {'La', 'Ce', 'Sm'}
    sp_block = element in {'Cs', 'Ba', 'Mg', 'Ca', 'Li', 'Na', 'Al', 'Pb'}
    residual = None if f_block else METALS[element] - oxidation
    return (None if sp_block else residual), residual


@dataclass(frozen=True)
class ElementRules:
    version: str
    valences: Mapping
    metals: Mapping
    main: frozenset
    oxidation: Mapping
    metal_counts: Callable
    view_type: type
    local_profile: Callable | None = None
    new_local: frozenset = frozenset()
    supplemental_profile: Callable | None = None


LEGACY = ElementRules('umr-elements/0.9', VALENCE, METALS, frozenset(MAIN), OXIDATION, legacy_metal_counts, LewisView)
_rules = ContextVar('umr_element_rules', default=LEGACY)
current_elements = _rules.get


@contextmanager
def use_elements(rules):
    token = _rules.set(rules)
    try:
        yield
    finally:
        _rules.reset(token)


def charge(a,b):
    valences = current_elements().valences
    if a.element not in valences or a.nonbonding_electrons is None:raise Unsupported('Electron ownership unknown; no zero default')
    return valences[a.element]-a.nonbonding_electrons-b
def profile(view,a):
    rules = current_elements()
    if a.element in rules.new_local:
        return rules.local_profile(view, a)
    if rules.supplemental_profile is not None:
        extra = rules.supplemental_profile(view, a)
        if extra is not None:
            return extra
    if a.element=='Sn':
        b=bond_valences(view)[a.atom_id];ats=atoms_by_id(view)
        ns=[ats[x.b if x.a==a.atom_id else x.a].element for x in view.bonds if a.atom_id in {x.a,x.b}]
        if (b,a.nonbonding_electrons,a.radical_electrons,a.formal_charge,a.hydrogens)==(2,2,0,0,0) and len(ns)==2 and set(ns)<={'Cl','Br','I'}:
            return {'atom_id':a.atom_id,'element':'Sn','bond_order_sum':2,'N':2,'r':0,'formal_charge':0,'shell_electrons':6,'rule':'SnII_dihalide_paired/0.9-r10'}
        if (b,a.nonbonding_electrons,a.radical_electrons,a.formal_charge,a.hydrogens)==(4,0,0,0,0) and len(ns)==4 and set(ns)<={'C','O','Cl','Br','I'}:
            return {'atom_id':a.atom_id,'element':'Sn','bond_order_sum':4,'N':0,'r':0,'formal_charge':0,'shell_electrons':8,'rule':'SnIV_tetrahedral_ownership/0.9-r10'}
    if a.element=='Ca' and rules is LEGACY:
        b=bond_valences(view)[a.atom_id]
        if (b,a.nonbonding_electrons,a.formal_charge,a.hydrogens)!=(0,0,2,0) or a.radical_electrons not in {0,None}:
            raise Unsupported('Ca scope is explicit Ca(II) ionic ownership; no Ca(0), Ca(I), covalent or surface model')
        return {'atom_id':a.atom_id,'element':'Ca','bond_order_sum':0,'N':0,'r':a.radical_electrons,'formal_charge':2,
            'shell_electrons':0,'rule':'CaII_ionic/0.9-r5','counting_scope':'valence electrons only; closed argon core excluded'}
    if a.element=='Si' and a.radical_electrons:
        b=bond_valences(view)[a.atom_id]
        if (b,a.nonbonding_electrons,a.radical_electrons,a.formal_charge)!=(3,1,1,0) or any(e.order!=1 for e in view.bonds if a.atom_id in {e.a,e.b}):raise Unsupported('Finite neutral three-coordinate silyl radical only')
        return {'atom_id':a.atom_id,'element':'Si','bond_order_sum':3,'N':1,'r':1,'formal_charge':0,'shell_electrons':7,'rule':'silyl_radical/0.9'}
    if a.element=='Si' and bond_valences(view)[a.atom_id]==5:
        ats=atoms_by_id(view);ns=[ats[b.b if b.a==a.atom_id else b.a].element for b in view.bonds if a.atom_id in {b.a,b.b}]
        if (a.nonbonding_electrons,a.radical_electrons,a.formal_charge,a.hydrogens)!=(0,0,-1,0) or len(ns)!=5:raise FormalError('Five-coordinate silicon ownership inconsistent')
        if not set(ns)<={'C','N','O','F','Cl'}:raise Unsupported('Five-coordinate silicon donor elements outside finite profile')
        return {'atom_id':a.atom_id,'element':'Si','bond_order_sum':5,'N':0,'r':0,'formal_charge':-1,'shell_electrons':10,'rule':'Si_5_coordinate_N_extension/0.9'}
    if a.element not in rules.metals:return profile08(view,a)
    b=bond_valences(view)[a.atom_id];n=a.nonbonding_electrons;r=a.radical_electrons
    if a.formal_charge!=charge(a,b) or n<0 or (r is not None and (r>n or (n-r)%2)):
        raise FormalError('Metal charge/local electron ownership mismatch')
    if a.hydrogens:raise Unsupported('Metal hydride must have explicit atom ID')
    if a.element=='Ce' and a.formal_charge==2:return profile08(view,a)
    if a.element in {'La','Ce','Sm'} and b:raise Unsupported('Covalent f-shell partition remains outside this adapter')
    return {'atom_id':a.atom_id,'element':a.element,'bond_order_sum':b,'N':n,'r':r,'formal_charge':a.formal_charge,
            'shell_electrons':n+2*b,'rule':'explicit_ownership/'+rules.version.split('/')[-1],'ground_state':'unknown'}
