"""Finite ordered-topology families. Atomic commit uses the existing ledger."""
from collections import Counter
from ..ast import Pair, Split, Shift
from ..chemistry import FormalError, Unsupported, atoms_by_id, order
from ..v02.pericyclic import compile_pericyclic as legacy

SOURCE = 'https://goldbook.iupac.org/terms/view/S05660'
PROFILES = {
 'COPE_3_3_EXTENDED': ('33', ('CCCCCC',), 'carbon [3,3]; O/O-/OSi or epoxide substituents remain explicitly in the ledger'),
 'CLAISEN_AROMATIC_3_3': ('33', ('CCOCCC',), 'O-linked [3,3], including cyclic vinyl and a specified aromatic Kekule view'),
 'OVERMAN_3_3': ('33', ('NCOCCC',), 'imidate ordered N=C-O-C-C=C; evidence-bound formal profile'),
 'FISCHER_3_3': ('33', ('CCNNCC',), 'enehydrazine N-N chain, actual protonation must pass the ledger'),
 'BARTOLI_3_3': ('33', ('CCNOCC',), 'ordered N-O-linked [3,3] hypothesis; no inference of full named mechanism'),
 'YLIDE_N_2_3': ('23', ('CNCCC','ONCCC'), 'ammonium ylide or amine oxide, negative terminal donor and positive N'),
 'YLIDE_S_2_3': ('23', ('CSCCC',), 'sulfonium ylide, negative carbon and positive S'),
 'WITTIG_O_2_3': ('23', ('COCCC',), 'oxygen-linked carbanion; product alkoxide retained'),
 'MISLOW_EVANS_2_3': ('23oxo', ('OSCCC',), 'sulfoxide S=O convention; sulfenate product'),
 'HETERO_DA_4_2': ('42', ('CCCCCO','CCCCOC','OCCCCC','CCCOCC'), 'one O heterodiene or O dienophile; no metal or ketene [2+2] inference'),
 'OZONE_3_2': ('32', ('OOOCC',), 'O(-)-O(+)=O + carbon alkene; primary ozonide/reverse only'),
 'CARBONYL_OXIDE_3_2': ('32', ('OOCCO','OOCOC'), 'O(-)-O(+)=C with carbonyl; orientation must be declared'),
 'HETERO_ELECTROCYCLIC_6PI': ('ec6', ('NCCCCC','CCCCCN'), 'one terminal N in a six-electron chain; not a generic Zincke completion'),
 'SULFOLENE_SO2_EXTRUSION': ('cheletropic', ('SCCCC',), 'one 3-sulfolene with two explicit S=O bonds; SO2 retained; formal thermal extrusion/reverse topology only'),
}


def inverse(a):
    if a.op == 'PAIR': return Split(a=a.b,b=a.a,target=a.a)
    if a.op == 'SPLIT': return Pair(a=a.target,b=a.b if a.target==a.a else a.a)
    return Shift(a=a.c,b=a.b,c=a.a)


def compile_pericyclic(source,event):
    if event.template_id not in PROFILES:
        return legacy(source,event)
    if event.template_version != '0.8': raise Unsupported('v0.8 family requires template_version 0.8')
    return compile_profile(source,event,profiles=PROFILES,template_version='0.8',
        rule_version='umr-pericyclic/0.8-r1',neutral_linkers={'WITTIG_O_2_3'})


def compile_profile(source,event,*,profiles,template_version,rule_version,neutral_linkers):
    """Shared ordered-topology compiler; registry/version ownership stays with callers."""
    if event.template_id not in profiles or event.template_version != template_version:
        raise Unsupported('Profile/version not registered by this compiler')
    family, patterns, scope = profiles[event.template_id]
    ids = event.bindings.get('ordered_atoms')
    if not isinstance(ids,list) or len(ids) != len(patterns[0]) or len(ids) != len(set(ids)):
        raise FormalError('Exact distinct ordered_atoms roles required')
    atoms = atoms_by_id(source)
    if not set(ids) <= atoms.keys() or set(event.region_atom_ids) != set(ids) or len(event.region_atom_ids)!=len(ids):
        raise FormalError('Participant IDs/region references mismatch')
    if ''.join(atoms[i].element for i in ids) not in patterns:
        raise Unsupported('Element ordering outside this registered family profile')
    if any(atoms[i].radical_electrons != 0 for i in ids):
        raise Unsupported('Closed-shell family does not authorize radicals or unknown local occupation')
    if event.mechanism_assumption not in {None,'concerted'}:
        raise FormalError('A stepwise assertion cannot execute as a concerted event')
    if event.conditions.get('geometry_status') == 'incompatible' or event.conditions.get('s_cis') is False:
        raise FormalError('Explicitly incompatible reaction geometry')
    activation = event.activation
    if activation in {'photochemical','hν'} and family != 'ec6':
        raise Unsupported('Photochemical mechanism is outside this thermal family')
    if family != 'ec6' and event.facial_mode not in {'unknown','supra_supra'}:
        raise Unsupported('Only registered supra/supra topology')
    a,b,c,d,e = ids[:5]
    if family == '33':
        f=ids[5]; required=[(a,b,2),(b,c,1),(c,d,1),(d,e,1),(e,f,2),(a,f,0)]
        actions=[Shift(a=b,b=a,c=f),Shift(a=f,b=e,c=d),Shift(a=d,b=c,c=b)]
    elif family in {'23','23oxo'}:
        required=[(a,b,2 if family=='23oxo' else 1),(b,c,1),(c,d,1),(d,e,2),(a,e,0)]
        actions=[Shift(a=b,b=a,c=e) if family=='23oxo' else Pair(a=a,b=e),Shift(a=e,b=d,c=c),Split(a=c,b=b,target=b)]
        if event.direction == 'forward' and family=='23':
            expected_b = 0 if event.template_id in neutral_linkers else 1
            if atoms[a].formal_charge != -1 or atoms[b].formal_charge != expected_b:
                raise Unsupported('Ylide donor/onium charge profile mismatch')
    elif family=='cheletropic':
        oxygens=[x.b if x.a==a else x.a for x in source.bonds if a in {x.a,x.b} and x.order==2 and atoms[x.b if x.a==a else x.a].element=='O']
        if len(oxygens)!=2 or any(atoms[i].formal_charge for i in [a,*oxygens]):raise Unsupported('Sulfolene adapter requires neutral explicit S(=O)2 representation')
        expected_n=0 if event.direction=='forward' else 2
        if atoms[a].nonbonding_electrons!=expected_n:raise FormalError('Sulfolene/SO2 sulfur electron occupancy mismatch')
        required=[(a,b,1),(b,c,1),(c,d,2),(d,e,1),(e,a,1)]
        actions=[Shift(a=a,b=b,c=c),Shift(a=c,b=d,c=e),Split(a=e,b=a,target=a)]
    elif family == '42':
        f=ids[5]; required=[(a,b,2),(b,c,1),(c,d,2),(e,f,2),(a,e,0),(d,f,0)]
        actions=[Shift(a=b,b=a,c=e),Shift(a=e,b=f,c=d),Shift(a=d,b=c,c=b)]
    elif family == '32':
        required=[(a,b,1),(b,c,2),(d,e,2),(a,e,0),(c,d,0)]
        if event.direction=='forward' and (atoms[a].formal_charge,atoms[b].formal_charge,atoms[c].formal_charge)!=(-1,1,0):
            raise Unsupported('Explicit 1,3-dipole charge arrangement required')
        actions=[Pair(a=a,b=e),Shift(a=e,b=d,c=c),Split(a=c,b=b,target=b)]
    else:
        f=ids[5]; required=[(a,b,2),(b,c,1),(c,d,2),(d,e,1),(e,f,2),(a,f,0)]
        actions=[Shift(a=b,b=a,c=f),Shift(a=f,b=e,c=d),Shift(a=d,b=c,c=b)]
        rotation='conrotatory' if activation in {'photochemical','hν'} else 'disrotatory'
        if activation != 'unknown' and event.rotation_mode not in {'unknown',rotation}:
            raise FormalError('Rotation contradicts the registered thermal/photo six-electron rule')
    delta=Counter()
    for x in actions:
        if x.op in {'SHIFT','SPLIT'}: delta[tuple(sorted((x.a,x.b)))]-=1
        if x.op=='SHIFT':delta[tuple(sorted((x.b,x.c)))]+=1
        if x.op=='PAIR':delta[tuple(sorted((x.a,x.b)))]+=1
    if event.direction not in {'forward','reverse'}:raise Unsupported('Unknown family direction')
    for x,y,n in required:
        expected=n+(delta[tuple(sorted((x,y)))] if event.direction=='reverse' else 0)
        if order(source,x,y)!=expected:raise FormalError(f'Family {event.direction} requires bond {x}-{y} order {expected}')
    if event.direction=='reverse':actions=[inverse(x) for x in actions]
    constraints=[]
    for constraint in event.stereo_constraints:
        if constraint.get('kind')=='facial' and constraint.get('mode') not in {'unknown','supra_supra'}:
            raise FormalError('Explicit facial assertion contradicts registered family')
        if constraint.get('kind')=='rotation' and family=='ec6' and activation!='unknown' and constraint.get('mode') not in {'unknown',rotation}:
            raise FormalError('Explicit rotation constraint contradicts registered six-electron rule')
        constraints.append({'declaration':constraint,'status':'unverified','not_a_product_assignment':True})
    return actions, {'family':family, 'profile':event.template_id, 'scope':scope, 'direction':event.direction,
        'topology':'matched', 'orbital_symmetry':'unknown' if activation=='unknown' else 'allowed_in_declared_model',
        'geometry':'unverified', 'stereo_constraints':constraints,'selectivity':'unknown',
        'source':'https://pubs.acs.org/doi/10.1021/jo00114a039' if family=='cheletropic' else SOURCE,
        'rule_version':rule_version, 'chemical_evidence':'requires case-specific source review'}
