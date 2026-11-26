"""Finite topology/occupation profiles added from inspected textbook evidence.

These selectors describe local events, not unrestricted named reactions. All
actions still go through the shared atomic ledger and source evidence gate.
"""
from ..ast import Pair, Split, Shift
from ..v07.ast import Fishhook
from ..v04.ast import ReactiveCenter
from ..chemistry import FormalError, Unsupported, atoms_by_id, order, bond_valences

SPECS = {
    'HALOCYCLOPROPYL_LI_ALPHA_ELIMINATION': dict(roles=['center','a','b','halogen','lithium'], event_type='POLAR', source='Kurti/Czako, printed p146 / PDF200; simultaneous C-Li/C-X cleavage to selected cyclopropylidene channel, not all carbenoids'),
    'CYCLOPROPYLIDENE_RING_OPENING': dict(roles=['center','a','b'], event_type='POLAR', source='Kurti/Czako, printed p146 / PDF200; free paired cyclopropylidene to allene'),
    'SELENIUM_DIOXIDE_ENE': dict(roles=['terminal','middle','allylic','hydrogen','selenium','oxo','acceptor_o'], event_type='PERICYCLIC', source='Kurti/Czako, printed p380 / PDF434; explicit SeO2 ene topology'),
    'ALLYLSELENINIC_2_3': dict(roles=['oxo','selenium','a','b','c','hydroxy_o','hydrogen'], event_type='PERICYCLIC', source='Kurti/Czako, printed p380 / PDF434; allylseleninic acid [2,3] shift'),
    'AMMONIUM_YLIDE_CAGE_FRAGMENTATION': dict(roles=['ylide','nitrogen','migrating'], event_type='RADICAL', source='Kurti/Czako, printed p434 / PDF488; ammonium ylide cage radical-pair channel; not concerted [1,2]'),
}
for spec in SPECS.values():
    spec['rule_revision'] = '0.9-r8'
    spec['geometry_status'] = 'unverified'

SPECS.update({name:dict(roles=roles.split(),event_type=kind,rule_revision='0.9-r9',
    geometry_status='unverified',source=source) for name,roles,kind,source in [
    ('CYCLOBUTENONE_OPENING','carbonyl alpha beta gamma oxygen','PERICYCLIC','Kurti/Czako PDF176: thermal cyclobutenone electrocyclic opening'),
    ('KETENE_ALKYNE_CAPTURE','carbonyl ketene_c oxygen a b','PERICYCLIC','Kurti/Czako PDF176: ketene/alkyne capture; not unrestricted thermal alkene [2+2]'),
    ('DIENYLKETENE_CLOSURE','a b c d e f oxygen','PERICYCLIC','Kurti/Czako PDF176: six-pi dienylketene closure'),
    ('DIOXIRANE_OXYGEN_TRANSFER','carbon retained_o transferred_o a b','POLAR','Kurti/Czako PDF464: concerted dioxirane oxygen transfer; transition-state selection uncomputed'),
    ('NA_PAIR_HALIDE_REDUCTION','first_na second_na carbon halogen','SET','Kurti/Czako PDF552: unresolved organosodium formation accounting; not a proven concerted elementary step'),
    ('KETYL_CARBONYL_REDUCTION','sodium radical_c carbonyl_c oxygen','SET','Kurti/Czako PDF58: selected acyloin diketone reduction branch'),
]})


def compile(source, event):
    name=event.general_event; spec=SPECS[name]; roles=event.general_roles
    ids=[roles[k] for k in spec['roles']]; atoms=atoms_by_id(source); vals=bond_valences(source)
    def elements(seq):
        if tuple(atoms[i].element for i in ids)!=seq:raise Unsupported('Element roles outside finite textbook profile')
    def bond(a,b,n=1):
        if order(source,a,b)!=n:raise FormalError(f'Required bond {a}-{b} order {n}')
    def neutral_closed(indices):
        if any(atoms[i].formal_charge or atoms[i].radical_electrons!=0 for i in indices):
            raise Unsupported('Profile requires neutral paired local occupation')
    allowed=lambda b: (name=='HALOCYCLOPROPYL_LI_ALPHA_ELIMINATION' and b.mode=='sigma_covalent'
        and b.metal_atom_id==roles['lithium'] and b.donor_atom_ids==[roles['center']])
    if any(set(ids)&set(b.donor_atom_ids) and not allowed(b) for b in source.metal_bindings):
        raise Unsupported('Coordinated reacting center requires an independent profile')
    if name in {'NA_PAIR_HALIDE_REDUCTION','KETYL_CARBONYL_REDUCTION'}:
        def f(src,dst):return Fishhook.model_validate(dict(op='FISHHOOK',source=src,destination=dst))
        def na(i):
            if (atoms[i].element,vals[i],atoms[i].nonbonding_electrons,atoms[i].radical_electrons,atoms[i].formal_charge)!=('Na',0,1,1,0):
                raise Unsupported('Isolated Na0 with one explicit ownership slot required; no bulk inference')
        if name=='NA_PAIR_HALIDE_REDUCTION':
            m,n,c,x=ids;na(m);na(n);bond(c,x)
            if atoms[c].element!='C' or atoms[x].element not in {'Cl','Br','I'}:raise Unsupported('Saturated C-halogen only')
            if (vals[c],atoms[c].nonbonding_electrons,atoms[c].radical_electrons,atoms[c].formal_charge)!=(4,0,0,0):raise FormalError('Neutral saturated carbon required')
            if any(b.order!=1 for b in source.bonds if c in (b.a,b.b)):raise Unsupported('Wurtz accounting fixture is saturated C-X only')
            neutral_closed([x])
            return [f(dict(kind='radical',atom_id=i,slot=0),dict(kind='lone_pair',atom_id=c)) for i in [m,n]]+[
                f(dict(kind='bond',a=c,b=x,pair='sigma',slot=k),dict(kind='lone_pair',atom_id=x)) for k in [0,1]],[]
        m,c,d,o=ids;na(m);bond(c,d);bond(d,o,2)
        if [atoms[i].element for i in [c,d,o]]!=['C','C','O']:raise Unsupported('Adjacent ketyl/carbonyl scope')
        if (vals[c],atoms[c].nonbonding_electrons,atoms[c].radical_electrons,atoms[c].formal_charge)!=(3,1,1,0):raise FormalError('Neutral ketyl carbon radical required')
        if not any(atoms[j].element=='O' and atoms[j].formal_charge==-1 and order(source,c,j)==1 for j in atoms):raise Unsupported('Ketyl anionic oxygen required')
        neutral_closed([d,o])
        return [f(dict(kind='radical',atom_id=i,slot=0),dict(kind='bond',a=c,b=d,pair='pi1')) for i in [m,c]]+[
            f(dict(kind='bond',a=d,b=o,pair='pi1',slot=k),dict(kind='lone_pair',atom_id=o)) for k in [0,1]],[]
    if name in {'CYCLOBUTENONE_OPENING','KETENE_ALKYNE_CAPTURE','DIENYLKETENE_CLOSURE','DIOXIRANE_OXYGEN_TRANSFER'}:
        neutral_closed(ids)
        if name=='CYCLOBUTENONE_OPENING':
            elements(('C','C','C','C','O'));c,a,b,g,o=ids
            for i,j,n in [(c,o,2),(c,a,1),(a,b,2),(b,g,1),(g,c,1)]:bond(i,j,n)
            if event.conditions.get('activation')!='thermal':raise Unsupported('Thermal activation declaration required')
            return [Shift(a=g,b=c,c=a),Shift(a=a,b=b,c=g)],[]
        if name=='KETENE_ALKYNE_CAPTURE':
            elements(('C','C','O','C','C'));c,k,o,a,b=ids
            for i,j,n in [(c,o,2),(c,k,2),(a,b,3),(k,b,0),(a,c,0)]:bond(i,j,n)
            return [Shift(a=c,b=k,c=b),Shift(a=b,b=a,c=c)],[]
        if name=='DIENYLKETENE_CLOSURE':
            elements(('C','C','C','C','C','C','O'));a,b,c,d,e,f,o=ids
            for i,j,n in [(a,o,2),(a,b,2),(b,c,1),(c,d,2),(d,e,1),(e,f,2),(a,f,0)]:bond(i,j,n)
            if event.conditions.get('activation')!='thermal':raise Unsupported('Thermal activation declaration required')
            return [Shift(a=b,b=a,c=f),Shift(a=f,b=e,c=d),Shift(a=d,b=c,c=b)],[]
        elements(('C','O','O','C','C'));c,o,t,a,b=ids
        for i,j,n in [(c,o,1),(o,t,1),(t,c,1),(a,b,2),(t,a,0),(t,b,0)]:bond(i,j,n)
        if any(vals[i]!=2 or atoms[i].nonbonding_electrons!=4 for i in [o,t]):raise FormalError('Paired dioxirane oxygens required')
        return [Shift(a=c,b=t,c=a),Shift(a=a,b=b,c=t),Shift(a=t,b=o,c=c)],[]
    if name.startswith('HALOCYCLO') or name=='CYCLOPROPYLIDENE_RING_OPENING':
        c,a,b=ids[:3]
        if any(atoms[i].element!='C' for i in [c,a,b]):raise Unsupported('Three-carbon ring required')
        for x,y in [(c,a),(c,b),(a,b)]:bond(x,y)
        neutral_closed([a,b])
        if name.startswith('HALOCYCLO'):
            x,li=ids[3:];bond(c,x);bond(c,li)
            if atoms[x].element not in {'Cl','Br'}:raise Unsupported('Only cyclopropyl chloride/bromide channel registered')
            neutral_closed([x])
            if (vals[c],atoms[c].formal_charge,atoms[c].nonbonding_electrons,atoms[c].radical_electrons)!=(4,0,0,0):
                raise FormalError('Sigma-bound halo-cyclopropyl carbon occupation required')
            if (atoms[li].element,vals[li],atoms[li].formal_charge,atoms[li].nonbonding_electrons)!=('Li',1,0,0):
                raise Unsupported('Finite sigma-bound Li(I) input required')
            return [Split(a=li,b=c,target=c),Split(a=c,b=x,target=x)],[]
        center=next((p for p in source.reactive_centers if p.atom_id==c),None)
        if center is None or center.profile!='carbene_paired' or center.environment!='free':
            raise Unsupported('Explicit free paired cyclopropylidene required')
        if (vals[c],atoms[c].formal_charge,atoms[c].nonbonding_electrons,atoms[c].radical_electrons)!=(2,0,2,0):
            raise FormalError('Cyclopropylidene occupation mismatch')
        if source.electronic_state.coupling=='high_spin':raise Unsupported('requires_spin_transition')
        return [Shift(a=a,b=b,c=c),Pair(a=c,b=a)],[]
    if name=='SELENIUM_DIOXIDE_ENE':
        elements(('C','C','C','H','Se','O','O'));a,b,c,h,se,o,x=ids
        neutral_closed(ids)
        for i,j,n in [(a,b,2),(b,c,1),(c,h,1),(se,o,2),(se,x,2),(a,se,0),(x,h,0)]:bond(i,j,n)
        if vals[se]!=4 or atoms[se].nonbonding_electrons!=2:raise FormalError('Molecular SeO2 occupation required')
        return [Shift(a=b,b=a,c=se),Shift(a=h,b=c,c=b),Shift(a=se,b=x,c=h)],[]
    if name=='ALLYLSELENINIC_2_3':
        elements(('O','Se','C','C','C','O','H'));o,se,a,b,c,x,h=ids
        neutral_closed(ids)
        for i,j,n in [(o,se,2),(se,a,1),(a,b,1),(b,c,2),(se,x,1),(x,h,1),(o,c,0)]:bond(i,j,n)
        if vals[se]!=4 or atoms[se].nonbonding_electrons!=2:raise FormalError('Allylseleninic acid occupation required')
        return [Shift(a=se,b=o,c=c),Shift(a=c,b=b,c=a),Split(a=a,b=se,target=se)],[]
    elements(('C','N','C'));c,n,m=ids
    bond(c,n);bond(n,m);bond(c,m,0)
    if [(vals[i],atoms[i].nonbonding_electrons,atoms[i].radical_electrons,atoms[i].formal_charge) for i in ids] != [(3,2,0,-1),(4,0,0,1),(4,0,0,0)]:
        raise FormalError('Exact closed-shell C-/N+/C ammonium ylide required')
    if any(a.radical_electrons!=0 for a in source.atoms) or source.electronic_state.coupling=='high_spin':
        raise Unsupported('Cage-pair scope starts from known paired inventory; no unrelated spin conversion')
    if (event.electronic_outcome is None or event.electronic_outcome.character!='open_shell_singlet'
            or event.electronic_outcome.coupling!='unknown' or not event.electronic_outcome.evidence):
        raise Unsupported('Source-declared singlet cage-pair model required; coupling dynamics not computed')
    def f(src,dst):return Fishhook.model_validate(dict(op='FISHHOOK',source=src,destination=dst))
    return [
        f(dict(kind='bond',a=n,b=m,pair='sigma',slot=0),dict(kind='atom',atom_id=m)),
        f(dict(kind='bond',a=n,b=m,pair='sigma',slot=1),dict(kind='lone_pair',atom_id=n)),
        f(dict(kind='lone_pair',atom_id=c,pair_index=0,slot=0),dict(kind='lone_pair',atom_id=n)),
        f(dict(kind='lone_pair',atom_id=c,pair_index=0,slot=1),dict(kind='atom',atom_id=c)),
    ],[]


def postprocess(target,event):
    if event.general_event=='HALOCYCLOPROPYL_LI_ALPHA_ELIMINATION':
        li=event.general_roles['lithium'];c=event.general_roles['center']
        target.metal_bindings=[b for b in target.metal_bindings if not (b.metal_atom_id==li and b.donor_atom_ids==[c])]
        from .metal_events import repartition
        repartition(target)
        target.reactive_centers.append(ReactiveCenter(center_id=event.event_id+'.cyclopropylidene',
            atom_id=event.general_roles['center'],species='carbene',profile='carbene_paired',
            environment='free',existence_basis='formal_model',reaction_context=event.mechanism_assumption,
            evidence=event.evidence))
    elif event.general_event=='CYCLOPROPYLIDENE_RING_OPENING':
        target.reactive_centers=[c for c in target.reactive_centers if c.atom_id!=event.general_roles['center']]
