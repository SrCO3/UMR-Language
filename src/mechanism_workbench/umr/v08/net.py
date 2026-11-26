"""Registered net changes use an internal ledger recipe, never microscopic arrows."""
from ..ast import Split, Shift, Pair
from ..chemistry import FormalError, Unsupported, atoms_by_id, order
from .transaction import transact

REGISTRY={
 'ALKENE_H2': {'roles':['c1','c2','h1','h2'], 'condition':'catalyst', 'scope':'one ordinary carbon alkene + explicit H2; Pd condition, no surface mechanism'},
 'PROTON_RELAY': {'roles':['donor','hydrogen','acceptor'], 'condition':'relay_medium', 'scope':'one explicit H+ donor/acceptor net pair transfer, no HAT/hydride/PCET'},
 'KETO_ENOL_RELAY': {'roles':['alpha_carbon','carbonyl_carbon','oxygen','hydrogen'], 'condition':'relay_medium', 'scope':'explicit C-H/O-H and pi-bond changes; water retained; no uncatalyzed concerted 1,3-H-shift claim'},
}


def execute(source,event):
    name=event.net_template
    if name not in REGISTRY:raise Unsupported('Unregistered net transformation; not a fallback for failed pathways')
    spec=REGISTRY[name];roles=event.net_roles
    if set(roles)!=set(spec['roles']) or len(set(roles.values()))!=len(roles):raise FormalError('Exact distinct net roles required')
    if not event.evidence or not event.conditions.get(spec['condition']):raise FormalError('Net source and explicit relevant conditions required')
    atoms=atoms_by_id(source)
    if not set(roles.values())<=atoms.keys():raise FormalError('Net role atom absent')
    if any(atoms[i].radical_electrons!=0 for i in roles.values()):raise Unsupported('Net profile is closed-shell; not HAT/PCET or ISC')
    if name in {'PROTON_RELAY','KETO_ENOL_RELAY'}:
        w=event.conditions.get('relay_water_oxygen_id')
        from ..chemistry import bond_valences
        if w not in atoms or w in roles.values() or (atoms[w].element,atoms[w].formal_charge,atoms[w].nonbonding_electrons,atoms[w].radical_electrons)!=('O',0,4,0) or bond_valences(source)[w]!=2:
            raise FormalError('Finite relay requires an explicit separate water species in the inventory')
        if any(atoms[b.b if b.a==w else b.a].element!='H' for b in source.bonds if w in {b.a,b.b}):raise Unsupported('First relay medium is explicit water only')
    if name=='ALKENE_H2':
        a,b,h,j=(roles[x] for x in spec['roles'])
        if [atoms[i].element for i in (a,b,h,j)]!=['C','C','H','H'] or order(source,a,b)!=2 or order(source,h,j)!=1:
            raise FormalError('Alkene hydrogenation requires the specified C=C and H-H inventory')
        if 'Pd' not in str(event.conditions['catalyst']):raise Unsupported('First net hydrogenation scope is declared Pd/H2')
        if not any(a.element=='Pd' for a in source.atoms):raise FormalError('Declared Pd catalyst must remain in the explicit material inventory')
        recipe=[Shift(a=b,b=a,c=h),Shift(a=h,b=j,c=b)]
    elif name=='KETO_ENOL_RELAY':
        a,c,o,h=(roles[x] for x in spec['roles'])
        if [atoms[i].element for i in (a,c,o,h)]!=['C','C','O','H'] or any(atoms[i].formal_charge for i in (a,c,o,h)):
            raise Unsupported('Neutral keto/enol C-C-O-H profile required')
        keto={(a,c):1,(c,o):2,(a,h):1,(o,h):0}
        enol={(a,c):2,(c,o):1,(a,h):0,(o,h):1}
        expected=keto if event.direction=='forward' else enol if event.direction=='reverse' else None
        if expected is None or any(order(source,*pair)!=n for pair,n in expected.items()):raise FormalError('Keto/enol roles and specified hydrogen do not match')
        recipe=[Shift(a=o,b=c,c=a),Shift(a=a,b=h,c=o)]
        if event.direction=='reverse':
            from .pericyclic import inverse
            recipe=[inverse(x) for x in recipe]
    else:
        donor,h,acceptor=(roles[x] for x in spec['roles'])
        if atoms[h].element!='H' or order(source,donor,h)!=1 or order(source,acceptor,h)!=0:
            raise FormalError('Proton relay requires explicit bonded donor H and unbonded acceptor')
        if atoms[donor].element not in {'O','N','S'} or atoms[acceptor].element not in {'O','N','S'}:
            raise Unsupported('First relay profile limited to heteroatom acid/base sites; tautomerism needs its own recipe')
        recipe=[Pair(a=acceptor,b=h),Split(a=h,b=donor,target=donor)]
    target=transact(source,recipe,event.product_state_id,event.product_view_id)
    return target, {'event_id':event.event_id,'event_type':event.event_type,'execution_mode':'net_transformation',
         'net_template':name,'scope':spec['scope'],'formal_status':'net_ledger_passed','atomic_commit':True,
         'compiled_actions':[], 'internal_balance_recipe':[x.model_dump(mode='json') for x in recipe],
         'microscopic_mechanism_verified':False,'chemical_feasibility':'unknown','human_approved':False,
         'stereochemistry':'source_declarations_separate; no automatic syn/E/Z product prediction'}
