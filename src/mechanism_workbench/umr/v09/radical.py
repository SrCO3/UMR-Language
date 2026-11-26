"""Role/topology-bound open-shell events, independent of reaction names."""
from ..chemistry import FormalError,Unsupported,atoms_by_id,order
from ..v02.electrons import atom_to_bond,bond_to_atom,bond_to_bond
SPECS={
 'ENEDIYNE_DIRADICAL':('a b c d e f','https://doi.org/10.1021/ja00757a071'),
 'DIVINYL_BRIDGING_DIRADICAL':('a b c d e','source-declared photochemical divinylmethane first event; no ground-state inference'),
 'CHARGED_BOND_HOMOLYSIS':('a b','charged N-C/N-halogen local homolysis, source-required'),
 'HAT_ROLE':('radical hydrogen donor','explicit H/D/T; local occupation, not acidity'),
 'HETERO_BETA_SCISSION':('radical adjacent departing','C/N/O/S radical adjacent to C/N; source-bound'),
 'ACYL_DENITROGENATION':('radical nitrogen carbon oxygen','acyl-N=N radical fragmentation; distinct from aryl profile'),
 'DISSOCIATIVE_NA_REDUCTION':('sodium carbon sulfur','source-declared molecular Na0 electron donor and sulfone fragmentation; not amalgam surface dynamics'),
 'NA_RADICAL_REDUCTION':('sodium carbon','known sodium unpaired slot pairs a carbon radical; no bulk substitution'),
 'LI_RADICAL_REDUCTION':('lithium carbon','known finite Li0 single-electron ownership slot; no bulk or solvation dynamics'),
 'LI_CONJUGATED_REDUCTION':('metal a b c d','finite Li0 donor, declared conjugated region and localized radical-anion contributor'),
 'METAL_PI_REDUCTION':('metal a b','known donor slot reduces C=O/C=C; no solvated-electron dynamics'),
 'METAL_CONJUGATED_REDUCTION':('metal a b c d','known Na slot and declared conjugated reduction contributor'),
 'METAL_DISSOCIATIVE_HALIDE_SET':('metal carbon halogen','known Cu0 radical or CrII pair slot; explicit radical/halide products'),
 'HETERO_RADICAL_PI_ADDITION':('radical attacked distal','C/Sn/Si radical onto explicit C=S/C=O region'),
 'AZO_INITIATOR_N2':('a n1 n2 b','source-declared azo fragmentation; N2 and two radicals retained'),
 'THIOCARBONYL_RELAY_FRAGMENTATION':('radical sulfur carbon nitrogen oxygen','Sn radical, thiocarbonyl and N-O fragmentation; original source role-bound'),
}
def compile(source,event):
    name=event.general_event;names=SPECS[name][0].split();r=event.general_roles;atoms=atoms_by_id(source)
    if set(r)!=set(names) or len(set(r.values()))!=len(r) or not set(r.values())<=atoms.keys():raise FormalError('Missing/duplicate radical roles')
    ids=[r[n] for n in names]
    def bond(a,b,n=1):
        if order(source,a,b)!=n:raise FormalError(f'Expected bond {a}-{b} order {n}')
    def rad(a):
        if atoms[a].radical_electrons!=1:raise FormalError('Occupied single-electron slot required')
    if name=='HETERO_RADICAL_PI_ADDITION':
        r,a,b=ids;rad(r);bond(a,b,2);bond(r,a,0)
        if atoms[r].element not in {'C','Sn','Si'} or {atoms[a].element,atoms[b].element} not in [{'C','S'},{'C','O'}]:raise Unsupported('Finite radical thiocarbonyl/carbonyl addition')
        return [atom_to_bond(r,r,a),bond_to_bond(a,b,'pi1',0,r,a),bond_to_atom(a,b,'pi1',1,b)],[]
    if name=='AZO_INITIATOR_N2':
        a,n,m,b=ids
        if [atoms[i].element for i in ids]!=['C','N','N','C']:raise Unsupported('Finite dialkyl azo precursor required')
        for x,y,n0 in [(a,n,1),(n,m,2),(m,b,1)]:bond(x,y,n0)
        return [bond_to_atom(a,n,'sigma',0,a),bond_to_bond(a,n,'sigma',1,n,m,'pi2'),bond_to_bond(m,b,'sigma',0,n,m,'pi2'),bond_to_atom(m,b,'sigma',1,b)],[n,m]
    if name=='THIOCARBONYL_RELAY_FRAGMENTATION':
        r,s,c,n,o=ids;rad(r)
        if [atoms[i].element for i in ids]!=['Sn','S','C','N','O']:raise Unsupported('Finite Sn/thiocarbonyl/N-O relay scope')
        for a,b,n0 in [(s,c,2),(c,n,1),(o,n,1),(r,s,0)]:bond(a,b,n0)
        return [atom_to_bond(r,r,s),bond_to_bond(c,s,'pi1',0,r,s),bond_to_bond(c,s,'pi1',1,n,c,'pi1'),bond_to_bond(o,n,'sigma',0,n,c,'pi1'),bond_to_atom(o,n,'sigma',1,o)],[]
    if name.startswith('METAL_') or name=='LI_CONJUGATED_REDUCTION':
        from ..v07.ast import Fishhook
        from .coordination import checks
        m=ids[0];counts={c['metal_atom_id']:c for c in checks(source)[1]}
        def lp(src,i):return Fishhook.model_validate({'op':'FISHHOOK','source':src,'destination':{'kind':'lone_pair','atom_id':i}})
        donor={'kind':'radical','atom_id':m,'slot':0}
        if name=='METAL_DISSOCIATIVE_HALIDE_SET':
            m,c,x=ids;bond(c,x)
            if atoms[c].element!='C' or atoms[x].element not in {'Cl','Br','I'}:raise Unsupported('Explicit C-halogen SET target required')
            if atoms[m].element=='Cu' and counts[m]['oxidation_state']==0:rad(m)
            elif atoms[m].element=='Cr' and counts[m]['oxidation_state']==2 and atoms[m].radical_electrons is not None and atoms[m].nonbonding_electrons-atoms[m].radical_electrons>=2:
                donor={'kind':'lone_pair','atom_id':m,'pair_index':0,'slot':0}
            else:raise Unsupported('Only declared Cu0 radical or known CrII paired-slot donor supported')
            return [lp(donor,x),bond_to_atom(c,x,'sigma',0,c),lp({'kind':'bond','a':c,'b':x,'pair':'sigma','slot':1},x)],[]
        allowed='Li' if name=='LI_CONJUGATED_REDUCTION' else 'Na'
        if atoms[m].element!=allowed or atoms[m].formal_charge!=0 or atoms[m].nonbonding_electrons!=1 or any(m in (z.a,z.b) for z in source.bonds):raise Unsupported('Finite isolated '+allowed+'0 ownership donor only; not bulk/solvated electron')
        rad(m);a,b=ids[1:3];bond(a,b,2)
        if atoms[a].element!='C' or atoms[b].element not in {'C','O'}:raise Unsupported('Finite C=C/C=O reduction')
        result=[lp(donor,b),lp({'kind':'bond','a':a,'b':b,'pair':'pi1','slot':0},b)]
        if name=='METAL_PI_REDUCTION':return result+[bond_to_atom(a,b,'pi1',1,a)],[]
        c,d=ids[3:];bond(a,d);bond(c,d,2)
        if any(atoms[i].element!='C' for i in [a,b,c,d]):raise Unsupported('Conjugated reduction is carbon pi region only')
        return result+[bond_to_bond(a,b,'pi1',1,a,d,'pi1'),bond_to_bond(c,d,'pi1',0,a,d,'pi1'),bond_to_atom(c,d,'pi1',1,c)],[]
    if name in {'DISSOCIATIVE_NA_REDUCTION','NA_RADICAL_REDUCTION','LI_RADICAL_REDUCTION'}:
        from ..v07.ast import Fishhook
        m,c=ids[:2];rad(m)
        allowed='Li' if name=='LI_RADICAL_REDUCTION' else 'Na'
        if atoms[m].element!=allowed or atoms[m].formal_charge!=0 or atoms[m].nonbonding_electrons!=1 or atoms[c].element!='C':raise Unsupported('Explicit molecular '+allowed+'0 single-electron donor and carbon required')
        if any(b.a==m or b.b==m for b in source.bonds):raise Unsupported('Bound/surface sodium requires a separate external model')
        def to_lp(src,d):return Fishhook.model_validate({'op':'FISHHOOK','source':src,'destination':{'kind':'lone_pair','atom_id':d}})
        donor={'kind':'radical','atom_id':m,'slot':0}
        if name in {'NA_RADICAL_REDUCTION','LI_RADICAL_REDUCTION'}:
            rad(c)
            return [to_lp({'kind':'radical','atom_id':c,'slot':0},c),to_lp(donor,c)],[]
        s=ids[2];bond(c,s)
        if atoms[s].element!='S' or atoms[s].formal_charge or len([a for a in atoms if atoms[a].element=='O' and order(source,s,a)==2])!=2:raise Unsupported('Finite neutral sulfone C-S dissociative reduction')
        return [bond_to_atom(c,s,'sigma',0,c),to_lp({'kind':'bond','a':c,'b':s,'pair':'sigma','slot':1},s),to_lp(donor,s)],[]
    if name=='ENEDIYNE_DIRADICAL':
        a,b,c,d,e,f=ids
        if any(atoms[i].element!='C' for i in ids):raise Unsupported('Enediyne profile requires six carbon centers')
        for x,y,n in [(a,b,3),(b,c,1),(c,d,2),(d,e,1),(e,f,3),(a,f,0)]:bond(x,y,n)
        return [bond_to_bond(a,b,'pi2',0,a,f),bond_to_atom(a,b,'pi2',1,b),bond_to_bond(e,f,'pi2',0,a,f),bond_to_atom(e,f,'pi2',1,e)],[]
    if name=='DIVINYL_BRIDGING_DIRADICAL':
        a,b,c,d,e=ids
        if any(atoms[i].element!='C' for i in ids):raise Unsupported('Divinyl first-event profile requires five carbons')
        if event.conditions.get('activation') not in {'photochemical','excited_state'}:raise Unsupported('Explicit excitation declaration required; photodynamics not computed')
        for x,y,n in [(a,b,2),(b,c,1),(c,d,1),(d,e,2),(b,d,0)]:bond(x,y,n)
        return [bond_to_bond(a,b,'pi1',0,b,d),bond_to_atom(a,b,'pi1',1,a),bond_to_bond(d,e,'pi1',0,b,d),bond_to_atom(d,e,'pi1',1,e)],[]
    if name=='CHARGED_BOND_HOMOLYSIS':
        a,b=ids;bond(a,b)
        if atoms[a].element!='N' or atoms[b].element not in {'C','Cl','Br','I'} or atoms[a].formal_charge!=1 or any(atoms[i].radical_electrons!=0 for i in ids):raise Unsupported('Charged homolysis profile is saturated ammonium N-C/N-halogen')
        return [bond_to_atom(a,b,'sigma',0,a),bond_to_atom(a,b,'sigma',1,b)],[]
    if name=='HAT_ROLE':
        r,h,d=ids;rad(r);bond(d,h);bond(r,h,0)
        if atoms[h].element!='H':raise FormalError('Explicit transferring hydrogen required')
        if atoms[r].element not in {'C','N','O','S','Sn','Se'} or atoms[d].element not in {'C','N','O','S','Sn','Si','Se'}:raise Unsupported('HAT element capabilities not registered')
        return [atom_to_bond(r,r,h),bond_to_bond(d,h,'sigma',0,r,h),bond_to_atom(d,h,'sigma',1,d)],[]
    if name=='HETERO_BETA_SCISSION':
        r,a,b=ids;rad(r);bond(r,a);bond(a,b)
        if atoms[r].element not in {'C','N','O','S'} or atoms[a].element not in {'C','N','O'}:raise Unsupported('Hetero beta-scission element profile missing')
        return [atom_to_bond(r,r,a,'pi1'),bond_to_bond(a,b,'sigma',0,r,a,'pi1'),bond_to_atom(a,b,'sigma',1,b)],[]
    r,n,c,o=ids;rad(r)
    if [atoms[i].element for i in ids]!=['N','N','C','O']:raise Unsupported('Acyl diazenyl roles required')
    bond(r,n,2);bond(n,c);bond(c,o,2)
    return [atom_to_bond(r,r,n,'pi2'),bond_to_bond(n,c,'sigma',0,r,n,'pi2'),bond_to_atom(n,c,'sigma',1,c)],[r,n]
