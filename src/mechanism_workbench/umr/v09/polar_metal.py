"""Finite element-role transactions inferred from typed actions, never reaction names.

A match is a candidate classification, not chemical validation. Bond preconditions,
finite oxidation profiles, donor reservations and the shared ledger remain mandatory.
"""
from dataclasses import dataclass
from itertools import permutations
from collections import Counter
import json
from ..ast import Pair,Split,Shift
from ..chemistry import FormalError,Unsupported,atoms_by_id,order,bond_valences,identity,electron_count
from .coordination import checks
from .elements import current_elements, METALS
from .ast import RULE_VERSION

@dataclass(frozen=True)
class Recipe:
    elements: dict
    actions: tuple
    bonds: tuple
    metal_states: dict
    release_roles: tuple=()
    scope: str='finite localized formal channel; geometry and feasibility unknown'

# P=PAIR; S=SPLIT with target specified; T=SHIFT. Every recipe has exact roles.
SPECS={
 'OSMATE_CYCLOADDITION':Recipe({'m':{'Os'},'o1':{'O'},'a':{'C'},'b':{'C'},'o2':{'O'}},
    (('T','m','o1','a'),('T','a','b','o2'),('S','o2','m','m')),
    (('m','o1',2),('a','b',2),('m','o2',2),('o1','a',0),('b','o2',0)),{'m':{('Os',8)}},
    scope='OsO4/alkene cyclic osmate formation; paired grouped [3+2] channel, syn geometry unresolved in 2D; template-repair/1.0'),
 'OSMATE_WATER_EXCHANGE':Recipe({'m':{'Os'},'o':{'O'},'w':{'O'},'h':{'H'}},
    (('P','w','m'),('S','m','o','o'),('P','o','h'),('S','h','w','w')),
    (('m','o',1),('w','h',1),('m','w',0),('o','h',0)),{'m':{('Os',6)}},
    scope='finite hydrolysis of an OsVI osmate ester with explicit water; retains Os and hydroxyl ligands, no redox or catalyst-regeneration claim; template-repair/1.0'),
 'TIN_LITHIUM_SIGMA_EXCHANGE':Recipe({'m':{'Sn'},'r':{'C'},'li':{'Li'},'alkyl':{'C'}},
    (('T','m','r','li'),('T','li','alkyl','m')),
    (('m','r',1),('li','alkyl',1),('r','li',0),('m','alkyl',0)),{'li':{('Li',1)}},
    scope='finite covalent RLi / neutral R4Sn paired exchange model; aggregation and kinetics unknown'),
 'LITHIUM_ZINC_SIGMA_EXCHANGE':Recipe({'li':{'Li'},'r':{'C'},'z':{'Zn'},'x':{'Cl','Br','I'}},
    (('T','li','r','z'),('T','z','x','li')),
    (('li','r',1),('z','x',1),('r','z',0),('li','x',0)),{'li':{('Li',1)},'z':{('Zn',2)}},
    scope='finite RLi / ZnX2 ligand exchange; no metal oxidation-state change'),
 'CU_ENONE_CAPTURE':Recipe({'m':{'Cu'},'b':{'C'},'a':{'C'},'c':{'C'},'o':{'O'}},
    (('P','m','b'),('T','b','a','c'),('S','c','o','o')),(('a','b',2),('a','c',1),('c','o',2),('m','b',0)),{'m':{('Cu',1)}}),
 'CU_C_C_ELIMINATION':Recipe({'m':{'Cu'},'r':{'C','Si'},'b':{'C'}},
    (('T','m','r','b'),('S','b','m','m')),(('m','r',1),('m','b',1),('r','b',0)),{'m':{('Cu',3)}}),
 'VINYL_MG_CU_EXCHANGE':Recipe({'m':{'Mg'},'r':{'C'},'u':{'Cu'},'x':{'Br','Cl'}},
    (('T','m','r','u'),('T','u','x','m')),(('m','r',1),('u','x',1),('r','u',0),('m','x',0)),{'m':{('Mg',2)},'u':{('Cu',1)}}),
 'PD_C_X_OXIDATIVE_ADDITION':Recipe({'m':{'Pd'},'c':{'C'},'x':{'Cl','Br','I'}},
    (('P','m','c'),('T','c','x','m')),(('c','x',1),('m','c',0),('m','x',0)),{'m':{('Pd',0)}}),
 'PD_AMIDO_FORMATION':Recipe({'m':{'Pd'},'n':{'N'},'h':{'H'},'x':{'Cl','Br','I'},'base':{'O'}},
    (('P','n','m'),('S','m','x','x'),('P','base','h'),('S','h','n','n')),
    (('n','h',1),('m','x',1),('m','n',0),('base','h',0)),{'m':{('Pd',2)}}),
 'PD_C_N_REDUCTIVE_ELIMINATION':Recipe({'m':{'Pd'},'c':{'C'},'n':{'N'}},
    (('S','c','m','m'),('T','m','n','c')),(('m','c',1),('m','n',1),('c','n',0)),{'m':{('Pd',2)}}),
 'PD_C_H_REDUCTIVE_ELIMINATION':Recipe({'m':{'Pd'},'c':{'C'},'h':{'H'}},
    (('S','c','m','m'),('T','m','h','c')),(('m','c',1),('m','h',1),('c','h',0)),{'m':{('Pd',2)}}),
 'MN_HYPOCHLORITE_OXO':Recipe({'m':{'Mn'},'o':{'O'},'x':{'Cl'}},
    (('S','o','x','x'),('P','m','o')),(('m','o',1),('o','x',1)),{'m':{('Mn',3)}},
    scope='declared MnIII bound-hypochlorite paired oxidation channel; occupation and kinetics not inferred'),
 'SIGMA_X_HETEROLYSIS':Recipe({'m':set(METALS),'x':{'C','N','O','S','Cl','Br','I','H'}},
    (('S','m','x','x'),),(('m','x',1),),{'m':{('Li',1),('Mg',2),('Zn',2),('Pd',2),('Cr',3),('Cu',1),('Cu',2),('Al',3),('Hg',2)}}),
 'HALIDE_METAL_ASSOCIATION':Recipe({'m':{'Cr','Al','Mg','Zn','Cu'},'x':{'Cl','Br','I'}},
    (('P','x','m'),),(('x','m',0),),{'m':{('Cr',3),('Al',3),('Mg',2),('Zn',2),('Cu',1),('Cu',2)}}),
 'METAL_O_X_EXCHANGE':Recipe({'m':{'Mn','Cr','Ti','Pd'},'d':{'O'},'x':{'O','Cl','Br','I'}},
    (('P','d','m'),('S','m','x','x')),(('m','d',0),('m','x',1)),
    {'m':{('Mn',3),('Cr',6),('Ti',4),('Pd',2)}}),
 'COORDINATED_AMINE_DEPROTONATION_EXCHANGE':Recipe({'m':{'Cu'},'n':{'N'},'x':{'O'},'h':{'H'}},
    (('P','n','m'),('S','m','x','x'),('P','x','h'),('S','h','n','n')),
    (('n','h',1),('m','x',1),('n','m',0),('x','h',0)),{'m':{('Cu',2)}},('n',)),
 'AL_HYDRIDE_ALKYNE_ADDITION':Recipe({'m':{'Al'},'h':{'H'},'a':{'C'},'b':{'C'}},
    (('T','m','h','a'),('S','a','b','b')),(('m','h',1),('a','b',3),('h','a',0)),{'m':{('Al',3)}}),
 'ALH4_CARBONYL_TRANSFER':Recipe({'m':{'Al'},'h':{'H'},'c':{'C'},'o':{'O'}},
    (('T','m','h','c'),('S','c','o','o')),(('m','h',1),('c','o',2),('h','c',0)),{}),
 'AL_ESTER_HYDRIDE_CAPTURE':Recipe({'m':{'Al'},'h':{'H'},'c':{'C'},'o':{'O'}},
    (('T','m','h','c'),('S','c','o','o'),('P','o','m')),
    (('m','h',1),('c','o',2),('h','c',0),('o','m',0)),{'m':{('Al',3)}}),
 'AL_ESTER_HEMIACETAL_COLLAPSE':Recipe({'m':{'Al'},'o':{'O'},'c':{'C'},'leave':{'O'}},
    (('T','m','o','c'),('S','c','leave','leave'),('P','leave','m')),
    (('m','o',1),('o','c',1),('c','leave',1),('m','leave',0)),{'m':{('Al',3)}}),
 'AL_ESTER_ALKOXIDE_WORKUP':Recipe({'m':{'Al'},'o':{'O'},'w':{'O'},'h':{'H'}},
    (('P','w','m'),('S','m','o','o'),('P','o','h'),('S','h','w','w')),
    (('m','o',1),('w','h',1),('m','w',0),('h','o',0)),{'m':{('Al',3)}}),
 'HG_CARBON_SUBSTITUTION':Recipe({'m':{'Hg'},'c':{'C'},'a':{'C'},'x':{'O'}},
    (('T','m','c','a'),('S','a','x','x')),(('m','c',1),('a','x',1),('c','a',0)),{'m':{('Hg',2)}}),
 'SIGMA_BOND_METATHESIS':Recipe({'a':{'Mg','Pd','Sn','Si'},'b':{'C','Sn','S','H','Cl','Br','I'},'c':{'Mg','Pd','Sn','Si'},'d':{'C','Sn','S','H','Cl','Br','I'}},
    (('T','a','b','c'),('T','c','d','a')),(('a','b',1),('c','d',1),('b','c',0),('d','a',0)),{}),
 'HALOGEN_LITHIUM_EXCHANGE':Recipe({'li':{'Li'},'c':{'C'},'x':{'I','Br'},'r':{'C'}},
    (('T','li','c','x'),('T','x','r','li')),(('li','c',1),('x','r',1),('c','x',0),('r','li',0)),{'li':{('Li',1)}}),
 'ACTIVATED_B_SI_PD_TRANSFER':Recipe({'d':{'B','Si'},'r':{'B','C'},'m':{'Pd'},'x':{'Cl','Br','I'}},
    (('T','d','r','m'),('S','m','x','x')),(('d','r',1),('m','x',1),('r','m',0)),{'m':{('Pd',2)}}),
 'BASE_ASSISTED_PD_ELIMINATION':Recipe({'base':{'O'},'h':{'H'},'a':{'C'},'b':{'C'},'m':{'Pd'},'x':{'Cl','Br','I'}},
    (('P','base','h'),('T','h','a','b'),('S','b','m','m'),('S','m','x','x')),
    (('a','h',1),('a','b',1),('b','m',1),('m','x',1),('base','h',0)),{'m':{('Pd',2)}}),
 'CHROMATE_ESTER_ELIMINATION':Recipe({'base':{'O'},'h':{'H'},'c':{'C'},'o':{'O'},'m':{'Cr'}},
    (('P','base','h'),('T','h','c','o'),('S','o','m','m')),
    (('c','h',1),('c','o',1),('o','m',1),('base','h',0)),{'m':{('Cr',6)}}),
 'PHOSPHINE_PD_O_TRANSFER':Recipe({'p':{'P'},'o':{'O'},'m':{'Pd'}},
    (('P','p','o'),('S','o','m','m')),(('p','o',0),('o','m',1)),{'m':{('Pd',2)}},('p',)),
 'TI_OXAMETALLACYCLE_CLEAVAGE':Recipe({'m':{'Ti'},'a':{'C'},'b':{'C'},'o':{'O'}},
    (('T','m','a','b'),('T','b','o','m')),(('m','a',1),('a','b',1),('b','o',1),('o','m',1)),{'m':{('Ti',4)}}),
 'OS_IMIDO_OXO_ADDITION':Recipe({'m':{'Os'},'n':{'N'},'a':{'C'},'b':{'C'},'o':{'O'}},
    (('T','m','n','a'),('T','a','b','o'),('S','o','m','m')),
    (('m','n',2),('a','b',2),('o','m',2),('n','a',0),('b','o',0)),{'m':{('Os',8)}}),
 'MG_N_O_FRAGMENT_EXCHANGE':Recipe({'n':{'N'},'o':{'O'},'x':{'O'},'m':{'Mg'}},
    (('P','o','n'),('S','n','x','x'),('S','m','o','o'),('P','x','m')),
    (('n','o',1),('n','x',1),('o','m',1),('x','m',0)),{'m':{('Mg',2)}}),
 'MG_NITROSO_O_ADDITION':Recipe({'m':{'Mg'},'r':{'C'},'o':{'O'},'n':{'N'}},
    (('T','m','r','o'),('S','o','n','n'),('P','n','m')),
    (('m','r',1),('o','n',2),('r','o',0),('n','m',0)),{'m':{('Mg',2)}}),
}

# The arrow contracts are shared with aluminate. Neutral dialkylaluminium
# has a different substrate profile and requires an explicitly named event.
for _step in ('HYDRIDE_CAPTURE','HEMIACETAL_COLLAPSE','ALKOXIDE_WORKUP'):
    SPECS['DIBAL_'+_step]=SPECS['AL_ESTER_'+_step]


def dibal_center(view, metal, atoms):
    def heavy(i):
        return {j for b in view.bonds if i in (b.a,b.b)
                for j in (b.b if b.a==i else b.a,) if atoms[j].element!='H'}

    def hydrogen_count(i):
        return atoms[i].hydrogens+sum(atoms[b.b if b.a==i else b.a].element=='H'
                                      for b in view.bonds if i in (b.a,b.b))

    ligands=[(b.b if b.a==metal else b.a,b.order) for b in view.bonds if metal in (b.a,b.b)]
    carbon=[i for i,_ in ligands if atoms[i].element=='C']
    if (atoms[metal].formal_charge!=0 or len(ligands)!=3 or len(carbon)!=2 or
            any(n!=1 or atoms[i].element not in {'C','O','H'} for i,n in ligands) or
            any(b.mode!='sigma_covalent' for b in view.metal_bindings if b.metal_atom_id==metal)):
        raise Unsupported('DIBAL requires neutral three-coordinate Al(III), two iBu ligands and one H or O')
    # Each ligand is Al–CH2–CH(CH3)2, not an arbitrary dialkylaluminium.
    for alpha in carbon:
        rest=heavy(alpha)-{metal}
        if len(rest)!=1:raise Unsupported('DIBAL requires two isobutyl ligands')
        beta=next(iter(rest));ends=heavy(beta)-{alpha}
        fragment={alpha,beta,*ends}
        if (len(ends)!=2 or any(atoms[i].element!='C' or atoms[i].formal_charge!=0 or atoms[i].radical_electrons!=0 for i in fragment) or
                hydrogen_count(alpha)!=2 or hydrogen_count(beta)!=1 or any(hydrogen_count(i)!=3 for i in ends) or
                any(heavy(i)!={beta} for i in ends) or
                any(b.order!=1 for b in view.bonds if fragment.intersection((b.a,b.b)))):
            raise Unsupported('Al-bound alkyl fragment is not isobutyl')


def ordinary_alkoxy(view, atoms, oxygen, acyl):
    neighbors=[(b.b if b.a==oxygen else b.a,b.order) for b in view.bonds if oxygen in (b.a,b.b) and acyl not in (b.a,b.b)]
    if len(neighbors)!=1 or atoms[oxygen].formal_charge!=0:return False
    other,n=neighbors[0]
    return n==1 and atoms[other].element=='C' and not any(
        z.order==2 and other in (z.a,z.b) and atoms[z.b if z.a==other else z.a].element=='O'
        for z in view.bonds)

def action_tuple(a):
    if a.op=='PAIR':return ('P',a.a,a.b)
    if a.op=='SPLIT':return ('S',a.b if a.target==a.a else a.a,a.target,a.target)
    if a.op=='SHIFT':return ('T',a.a,a.b,a.c)
    return None

def instantiate(spec,r):
    result=[]
    for t in spec.actions:
        x=[r[k] for k in t[1:]]
        result.append(Pair(a=x[0],b=x[1]) if t[0]=='P' else Split(a=x[0],b=x[1],target=x[2]) if t[0]=='S' else Shift(a=x[0],b=x[1],c=x[2]))
    return result

def role_match(v,spec,actions):
    if len(actions)!=len(spec.actions):return None
    actual=[action_tuple(a) for a in actions]
    if None in actual:return None
    ats=atoms_by_id(v)
    for seq in permutations(actual):
        r={};valid=True
        for template,item in zip(spec.actions,seq):
            if template[0]=='S' and template[1]==template[3]:template=('S',template[2],template[1],template[3])
            if template[0]!=item[0]:valid=False;break
            for role,i in zip(template[1:],item[1:]):
                if i not in ats or (role in r and r[role]!=i):valid=False;break
                r[role]=i
            if not valid:break
        if valid and len(set(r.values()))==len(r) and all(ats[r[k]].element in elements for k,elements in spec.elements.items()):return r
    return None

def candidates(v,e):
    if e.event_type not in {'POLAR','HETEROLYSIS'} or not e.actions:return []
    return [(name,r) for name,s in SPECS.items() if (not name.startswith(('AL_ESTER_','DIBAL_')) or e.general_event==name)
            and (r:=role_match(v,s,e.actions)) is not None]

def validate_roles(v,e,name,r):
    spec=SPECS[name];ats=atoms_by_id(v)
    if e.event_type not in {'POLAR','HETEROLYSIS'} or e.electron_transfers or e.electronic_outcome:
        raise Unsupported('Finite pair transaction cannot silently perform SET or a spin transition')
    if set(r)!=set(spec.elements) or len(set(r.values()))!=len(r) or not set(r.values())<=ats.keys():raise FormalError('Exact distinct finite transaction roles required')
    if any(ats[r[k]].element not in els for k,els in spec.elements.items()):raise FormalError('Element does not match finite transaction role')
    for boundary in v.material_boundaries:
        if set(r.values())&set(boundary.atom_ids) and (boundary.direction!='inside' or boundary.domain in {'bulk_metal','surface','external_electron_reservoir'}):
            raise Unsupported('Finite transaction cannot replace a bulk/surface or external exchange model')
    if not e.evidence or not(e.mechanism_assumption or e.explanation):raise FormalError('Explicit source/model evidence and explanation required')
    for a,b,n in spec.bonds:
        if order(v,r[a],r[b])!=n:raise FormalError(f'{name}: expected input bond {r[a]}-{r[b]} order {n}; preceding state does not supply it')
    counts={c['metal_atom_id']:c for c in checks(v)[1]}
    for role,allowed in spec.metal_states.items():
        i=r[role]
        if i not in counts or (ats[i].element,counts[i]['oxidation_state']) not in allowed:raise Unsupported(name+': metal oxidation profile outside finite scope')
    vals=bond_valences(v)
    if name=='TIN_LITHIUM_SIGMA_EXCHANGE':
        if (vals[r['m']],ats[r['m']].formal_charge,ats[r['m']].nonbonding_electrons)!=(4,0,0):
            raise Unsupported('Neutral tetracoordinate organotin required')
        if vals[r['li']]!=1 or ats[r['li']].formal_charge!=0:
            raise Unsupported('Explicit one-coordinate covalent RLi model required')
        if any((ats[r[k]].formal_charge,ats[r[k]].nonbonding_electrons,ats[r[k]].radical_electrons,vals[r[k]])!=(0,0,0,4) for k in ('r','alkyl')):
            raise Unsupported('Paired neutral carbon sigma ligands required')
    if name=='LITHIUM_ZINC_SIGMA_EXCHANGE':
        if vals[r['li']]!=1 or vals[r['z']]!=2 or ats[r['z']].formal_charge!=0:
            raise Unsupported('Finite RLi / two-coordinate neutral ZnX2 channel required')
        if any(ats[b.b if b.a==r['z'] else b.a].element not in {'Cl','Br','I'} for b in v.bonds if r['z'] in (b.a,b.b)):
            raise Unsupported('Zinc capture starts from a molecular zinc dihalide')
    if name.startswith(('AL_ESTER_','DIBAL_')):
        m=r['m']
        ns=[ats[b.b if b.a==m else b.a] for b in v.bonds if m in (b.a,b.b)]
        if name.startswith('DIBAL_'):
            dibal_center(v,m,ats)
        elif (ats[m].formal_charge!=-1 or len(ns)!=4 or vals[m]!=4 or
            any(a.element not in {'H','O'} for a in ns) or any(
                b.mode!='sigma_covalent' for b in v.metal_bindings if b.metal_atom_id==m)):
            raise Unsupported('Finite mononuclear tetrahydrido/alkoxy aluminate(III) required')
        if name.endswith('HYDRIDE_CAPTURE'):
            c=r['c'];neighbors=[ats[b.b if b.a==c else b.a] for b in v.bonds if c in (b.a,b.b)]
            if (ats[c].formal_charge!=0 or ats[r['o']].formal_charge!=0 or
                any(a.element not in {'C','H','O'} for a in neighbors) or
                any(a.element=='O' and a.atom_id!=r['o'] and not any(
                    a.atom_id in (b.a,b.b) and c not in (b.a,b.b) and ats[b.b if b.a==a.atom_id else b.a].element=='C'
                    for b in v.bonds) for a in neighbors)):
                raise Unsupported('Neutral ester, aldehyde or ketone required')
            if name.startswith('DIBAL_'):
                alkoxy=[a for a in neighbors if a.element=='O' and a.atom_id!=r['o']]
                if len(alkoxy)>1 or any(not ordinary_alkoxy(v,ats,a.atom_id,c) for a in alkoxy):
                    raise Unsupported('DIBAL ester profile requires one ordinary alkoxy leaving group')
        elif name.endswith('COLLAPSE'):
            if not any(r['c'] in (b.a,b.b) and ats[b.b if b.a==r['c'] else b.a].element=='H' for b in v.bonds):
                raise Unsupported('An ester-reduction hemiacetal with an explicit carbon hydrogen required')
        elif ats[r['w']].formal_charge!=0 or sum(order(v,r['w'],i)==1 and a.element=='H' for i,a in ats.items())!=2:
            raise Unsupported('An independently inventoried water with two numbered hydrogens is required')
        elif name.startswith('DIBAL_') and not any(
                b.order==1 and r['o'] in (b.a,b.b) and ats[b.b if b.a==r['o'] else b.a].element=='C' for b in v.bonds):
            raise Unsupported('DIBAL workup requires an aluminium-bound alkoxide')
    if name=='CU_ENONE_CAPTURE':
        m=r['m']
        bound=[ats[b.b if b.a==m else b.a] for b in v.bonds if m in (b.a,b.b)]
        if (not any(a.element in {'C','Si'} for a in bound) or ats[m].nonbonding_electrons!=10
                or ats[m].radical_electrons not in (0,None) or
                any(ats[r[k]].formal_charge or ats[r[k]].radical_electrons for k in ('a','b','c','o'))):
            raise Unsupported('Finite Cu(I) d10 organocopper and neutral closed-shell enone required')
    if name=='VINYL_MG_CU_EXCHANGE' and not any(
            b.order==2 and r['r'] in (b.a,b.b) and ats[b.b if b.a==r['r'] else b.a].element=='C' for b in v.bonds):
        raise Unsupported('Finite exchange requires vinylmagnesium and copper(I) halide')
    if name=='ALH4_CARBONYL_TRANSFER':
        m=r['m']
        ns=[ats[b.b if b.a==m else b.a] for b in v.bonds if m in (b.a,b.b)]
        if (ats[m].formal_charge!=-1 or vals[m]!=4 or len(ns)!=4 or any(a.element!='H' or vals[a.atom_id]!=1 for a in ns)
                or ats[r['c']].formal_charge!=0 or ats[r['o']].formal_charge!=0
                or any(ats[b.b if b.a==r['c'] else b.a].element not in {'C','H','O'} for b in v.bonds if r['c'] in (b.a,b.b))
                or sum(ats[b.b if b.a==r['c'] else b.a].element=='O' for b in v.bonds if r['c'] in (b.a,b.b))!=1):
            raise Unsupported('Finite tetrahydridoaluminate and neutral aldehyde/ketone required')
    if name=='PD_C_X_OXIDATIVE_ADDITION':
        if (ats[r['m']].nonbonding_electrons,ats[r['m']].radical_electrons)!=(10,0):
            raise Unsupported('Declare a filled Pd0 pair channel before using this finite profile')
        if any(ats[r[k]].formal_charge or ats[r[k]].radical_electrons for k in ('c','x')):
            raise Unsupported('Finite Pd oxidative addition requires neutral closed-shell C-X endpoints')
    if name=='PD_AMIDO_FORMATION':
        if (ats[r['n']].formal_charge,ats[r['n']].radical_electrons,vals[r['n']])!=(0,0,3) or (ats[r['base']].formal_charge,ats[r['base']].radical_electrons)!=(-1,0):
            raise Unsupported('Neutral secondary/primary amine and explicit oxygen anion base required')
    if name=='PD_C_N_REDUCTIVE_ELIMINATION' and (ats[r['n']].formal_charge,ats[r['n']].radical_electrons,vals[r['n']])!=(0,0,3):
        raise Unsupported('Covalent paired amido ligand required')
    if name in {'PD_C_H_REDUCTIVE_ELIMINATION','PD_C_N_REDUCTIVE_ELIMINATION'}:
        other='h' if name=='PD_C_H_REDUCTIVE_ELIMINATION' else 'n'
        if any(ats[r[k]].radical_electrons!=0 for k in ('c',other)):raise Unsupported('Paired carbon/hydride or amido ligands required')
        ids={b.binding_id for b in v.metal_bindings if b.metal_atom_id==r['m'] and set(b.donor_atom_ids)&{r['c'],r[other]}}
        center=next(c for c in v.metal_centers if c.atom_id==r['m'])
        if any(p.relation=='trans' and {p.first_binding_id,p.second_binding_id}<=ids for p in center.relative_positions):
            raise FormalError('Declared trans carbon/hydride requires isomerization before elimination')
    if name=='MN_HYPOCHLORITE_OXO':
        if e.general_event!=name or e.conditions.get('paired_donor_assumption')!='MnIII_two_electron_channel':
            raise Unsupported('Mn oxo formation requires explicit paired_donor_assumption; unknown spin is not a free pair')
        if (ats[r['o']].formal_charge,ats[r['o']].nonbonding_electrons,ats[r['o']].radical_electrons)!=(0,4,0):
            raise FormalError('Neutral paired bound hypochlorite required')
    if name=='SIGMA_BOND_METATHESIS':
        a,b,c,d=[r[k] for k in 'abcd'];ea,eb,ec,ed=[ats[i].element for i in [a,b,c,d]]
        valid=((ea,ec)==('Mg','Pd') and eb=='C' and ed in {'Cl','Br','I'} or
               (ea,ec)==('Pd','Mg') and eb in {'Cl','Br','I'} and ed=='C' or
               (ea,ec)==('Sn','Pd') and eb in {'C','Sn'} and ed in {'Cl','Br','I'} or
               (ea,ec)==('Pd','Sn') and eb in {'Cl','Br','I'} and ed in {'C','Sn'} or
               (ea,ec)==('Si','Pd') and eb=='H' and ed=='S' or
               (ea,ec)==('Pd','Si') and eb=='S' and ed=='H')
        if not valid:raise Unsupported('Only Mg/Pd, Sn/Pd transmetalation or Si-H/Pd-S metathesis registered')
        for i in [a,c]:
            if ats[i].element in current_elements().metals and counts[i]['oxidation_state']!=2:raise Unsupported('Metathesis requires MgII/PdII')
    if name=='ACTIVATED_B_SI_PD_TRANSFER':
        d=r['d']
        if (ats[d].element,vals[d],ats[d].formal_charge) not in {('B',4,-1),('Si',5,-1)}:raise Unsupported('Activated tetra-borate or penta-silicate required')
    if name=='HALIDE_METAL_ASSOCIATION' and (ats[r['x']].formal_charge,vals[r['x']])!=(-1,0):raise FormalError('Free explicit halide anion required')
    if name=='HG_CARBON_SUBSTITUTION' and (vals[r['a']]!=4 or any(b.order!=1 for b in v.bonds if r['a'] in {b.a,b.b})):raise Unsupported('Hg organic substitution limited to saturated electrophilic carbon')
    if name.startswith('MG_N') and any(ats[r[x]].formal_charge!=0 or ats[r[x]].radical_electrons!=0 for x in ('n','o')):
        raise Unsupported('Mg nitrogen/oxygen transaction requires the neutral paired input profile')
    for role in spec.release_roles:
        a=ats[r[role]]
        if a.radical_electrons!=0:raise Unsupported('Resolved paired donor occupation required')
    return counts

def execute(v,e,name=None,roles=None):
    from .execution import transact
    from .ledger import validate_view
    from .metal_events import repartition
    from ..v07.ledger import project_symbols
    validate_view(v)
    if name is None:
        matches=candidates(v,e)
        if len(matches)!=1:raise Unsupported('Finite metal transaction does not have one unambiguous role match')
        name,roles=matches[0]
    r=roles if roles is not None else e.general_roles
    before=validate_roles(v,e,name,r);spec=SPECS[name];actions=instantiate(spec,r)
    key=action_tuple
    if e.actions is not None and Counter(map(key,e.actions))!=Counter(map(key,actions)):raise FormalError('Actions disagree with finite transaction')
    metals={i for i in r.values() if atoms_by_id(v)[i].element in current_elements().metals}
    release={b.binding_id for b in v.metal_bindings if b.mode=='dative' and b.metal_atom_id in metals and set(b.donor_atom_ids)&{r[x] for x in spec.release_roles}}
    if any(len(b.donor_atom_ids)>1 for b in v.metal_bindings if b.binding_id in release):raise Unsupported('Partial chelate release needs an explicit donor-wise binding')
    def post(t):
        t.metal_bindings=[b for b in t.metal_bindings if b.binding_id not in release]
        repartition(t)
        replacements=[]
        if name=='TIN_LITHIUM_SIGMA_EXCHANGE':replacements=[(r['r'],r['m'],r['li']),(r['alkyl'],r['li'],r['m'])]
        if name=='LITHIUM_ZINC_SIGMA_EXCHANGE':replacements=[(r['r'],r['li'],r['z'])]
        for old in v.bonds:
            if old.order!=2 or old.stereo.get('status')!='specified':continue
            refs=list(old.stereo.get('reference_atoms',[]))
            if len(refs)!=2:continue
            affected=False
            for center,leaving,incoming in replacements:
                if center not in (old.a,old.b):continue
                j=0 if center==old.a else 1
                if refs[j]==leaving:refs[j]=incoming;affected=True
            if not affected:continue
            new=next((b for b in t.bonds if {b.a,b.b}=={old.a,old.b} and b.order==2),None)
            if new is not None and all(order(t,end,ref)==1 for end,ref in zip((old.a,old.b),refs)):
                new.stereo={**old.stereo,'reference_atoms':refs if new.a==old.a else refs[::-1],
                            'basis':'declared paired sigma-exchange vinyl retention; no product-derived assignment'}
    declared={r['m']} if name in {'MN_HYPOCHLORITE_OXO','CU_ENONE_CAPTURE'} else set()
    target=transact(v,actions,e.product_state_id,e.product_view_id,after=post,_released_binding_ids=release,_declared_pair_donors=declared)
    s,t=project_symbols(v),project_symbols(target)
    if identity(v)!=identity(target) or electron_count(s)!=electron_count(t) or sum(a.formal_charge for a in s.atoms)!=sum(a.formal_charge for a in t.atoms):raise FormalError('Finite metal transaction conservation failure')
    after={c['metal_atom_id']:c for c in checks(target)[1]}
    if name.startswith('DIBAL_') and after[r['m']]['oxidation_state']!=3:
        raise FormalError('DIBAL formal Al(III) changed during hydride transfer or hydrolysis')
    if name in {'CU_ENONE_CAPTURE','CU_C_C_ELIMINATION'}:
        if after[r['m']]['oxidation_state']!=({'CU_ENONE_CAPTURE':3,'CU_C_C_ELIMINATION':1}[name]):
            raise FormalError('Copper event oxidation-state closure failed')
    if name=='VINYL_MG_CU_EXCHANGE' and any(before[i]['oxidation_state']!=after[i]['oxidation_state'] for i in metals):
        raise FormalError('Vinyl transmetalation changed a metal oxidation state')
    if name in {'SIGMA_BOND_METATHESIS','HALOGEN_LITHIUM_EXCHANGE','ACTIVATED_B_SI_PD_TRANSFER'} and any(before[i]['oxidation_state']!=after[i]['oxidation_state'] for i in metals):raise FormalError('Transfer profile unexpectedly changed metal oxidation state')
    e.actions=actions
    return target,{'event_id':e.event_id,'event_type':e.event_type,'resolved_profile':name,'resolved_roles':r,'rule_version':RULE_VERSION,
      'atomic_commit':True,'execution_mode':'electron_actions','formal_status':'passed','compiled_actions':[a.model_dump(mode='json') for a in actions],
      'released_binding_ids':sorted(release),'chemical_feasibility':'unknown','geometry_status':'unverified',
      'occupation_assumption':e.conditions.get('paired_donor_assumption') if declared else None,
      'occupation_verified':False if declared else None,
      'source_scope':'typed topology/profile match; not a determination of the experimental pathway','ledger':validate_view(target)}
