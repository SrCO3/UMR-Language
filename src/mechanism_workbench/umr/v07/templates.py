"""Role-bound finite recipes. Sources support families, not these synthetic atom maps."""
from ..ast import Pair,Split,Shift
from ..chemistry import FormalError,Unsupported,atoms_by_id,order,bond_valences
from ..v02.electrons import atom_to_bond,bond_to_atom,bond_to_bond
from ..v05.ast import Ligand,MetalBinding,MetalCenter
from .ast import Fishhook,LonePairSource,LonePairDestination
from .elements import METAL_GROUP

SPECS={
 'LEWIS_ASSOCIATION':('donor metal','single neutral connected L donor; AlIII/FeIII/ZnII/HgII/PdII','https://doi.org/10.1515/pac-2015-1204'),
 'LEWIS_DISSOCIATION':('donor metal','inverse of a stored L dative binding','https://doi.org/10.1515/pac-2015-1204'),
 'LEWIS_ACYL_ACTIVATION':('carbon leaving metal donor','O-bound acyl chloride/AlCl3 localized acylium/AlCl4 model','https://openstax.org/books/organic-chemistry/pages/16-3-alkylation-and-acylation-of-aromatic-rings-the-friedel-crafts-reaction'),
 'LEWIS_BROMINATION':('adjacent ipso bromine leaving metal','aryl pi capture of FeIII-bound Br2; complete FeBr4 ledger','https://openstax.org/books/organic-chemistry/pages/16-1-electrophilic-aromatic-substitution-reactions-bromination'),
 'LEWIS_REAROMATIZATION':('metal leaving hydrogen ipso adjacent','aryl H transfer to a bound halide regenerates Lewis acid','https://openstax.org/books/organic-chemistry/pages/16-1-electrophilic-aromatic-substitution-reactions-bromination'),
 'AMMONIUM_YLIDE_2_3':('ylide nitrogen first middle last','allylic ammonium ylide thermal [2,3]; not [1,2] Stevens','https://doi.org/10.1021/jacs.6b11851'),
 'RADICAL_DENITROGENATION':('radical nitrogen carbon','Ar-N=N radical to Ar radical + retained N2','https://pmc.ncbi.nlm.nih.gov/articles/PMC12117561/'),
 'HAT_AROMATIZATION':('radical hydrogen first second third donor','formal cyclohexadienyl H transfer/rearomatization; source supports related radical disproportionation, not this exact partner','https://doi.org/10.1039/P29720001256'),
 'SN_HAT':('radical hydrogen tin','carbon radical + explicit R3Sn-H to C-H + R3Sn radical','https://doi.org/10.1021/ja00416a008'),
 'CU_DIAZONIUM_SET':('copper carbon nitrogen terminal','Cu0/I donation grouped with aryl diazonium C-N cleavage; formal net model, microscopic timing unverified','https://doi.org/10.1021/acs.orglett.0c04209'),
 'SAEGUSA_PALLADATION':('metal alpha carbon oxygen silicon leaving','formal grouped desilylation/palladation; no claim about microscopic timing','https://pubmed.ncbi.nlm.nih.gov/23256839/'),
 'PD_BETA_H_ELIMINATION':('metal alpha beta hydrogen','mononuclear PdII sigma alkyl; free alkene product; vacancy/syn geometry required or unknown','https://www.ilpi.com/organomet/betahydride.html'),
 'PD_H_ACETATE_ELIMINATION':('metal hydrogen oxygen','PdII hydride/acetate to Pd0 and acetic acid; stoichiometric endpoint','https://pubmed.ncbi.nlm.nih.gov/23256839/'),
 'P_SE_ACTIVATION':('phosphorus selenium cyanocarbon nitrogen','phosphine/selenocyanate activation with cyanide retained; formal candidate; original Grieco substrate scope is aryl','https://doi.org/10.1021/jo00870a052'),
 'SELENOXIDE_ELIMINATION':('selenium oxygen alpha beta hydrogen','five-center syn elimination of alkyl selenoxide; explicit H','https://doi.org/10.1016/S0040-4039(01)96098-8'),
}


def rebuild_bindings(view):
    """Derive connected ligand partitions after registered graph changes only."""
    atoms=atoms_by_id(view);mids={i for i,a in atoms.items() if a.element in METAL_GROUP}
    datives=[(b.donor_atom_ids[0],b.metal_atom_id) for b in view.metal_bindings if b.mode=='dative']
    if any(b.mode=='eta' for b in view.metal_bindings):raise Unsupported('Binding rebuild excludes eta/bridging changes')
    todo=set(atoms)-mids;parts=[]
    while todo:
        stack=[todo.pop()];seen=set(stack)
        while stack:
            i=stack.pop()
            for b in view.bonds:
                if i not in {b.a,b.b}:continue
                j=b.b if b.a==i else b.a
                if j in todo:todo.remove(j);seen.add(j);stack.append(j)
        parts.append(seen)
    ligands=[];bindings=[]
    for members in parts:
        contacts=[(d,m,'dative',0) for d,m in datives if d in members]
        contacts += [(b.a if b.a in members else b.b,b.b if b.a in members else b.a,'sigma_covalent',b.order)
                     for b in view.bonds if {b.a,b.b}&mids and {b.a,b.b}&members]
        if not contacts:continue
        if len(contacts)!=1:raise Unsupported('Registered rebuild requires monodentate, nonbridging ligands')
        d,m,mode,o=contacts[0];iq=sum(atoms[i].formal_charge for i in members)-o
        cls='L' if mode=='dative' and iq==0 else 'X' if o==1 and iq==-1 else 'X2' if o==2 and iq==-2 and len(members)==1 and atoms[d].element=='O' else 'unknown'
        if cls=='unknown':raise Unsupported('Ligand ionic class is outside the registered partition')
        lid=f'ligand.{min(members)}';ligands.append(Ligand(ligand_id=lid,atom_ids=sorted(members),donor_atom_ids=[d],ionic_class=cls,ionic_charge=iq,assignment_basis='registered event product graph'))
        bindings.append(MetalBinding(binding_id=f'bind.{m}.{d}',metal_atom_id=m,ligand_id=lid,donor_atom_ids=[d],mode=mode))
    view.ligands=ligands;view.metal_bindings=bindings
    for m in view.metal_centers:
        # Input assertions remain in the input. Product derived fields must be recomputed.
        m.oxidation_state=None;m.d_electrons=None;m.valence_electrons=None
        m.relative_positions=[];m.vacant_site_ids=None


def compile_template(view,event):
    name=event.rule_template
    if name not in SPECS:raise Unsupported('Unregistered v0.7 template')
    names,scope,source=SPECS[name];names=names.split();r=event.bindings;atoms=atoms_by_id(view)
    if set(r)!=set(names) or any(type(r[k]) is not int or r[k] not in atoms for k in names) or len(set(r.values()))!=len(r):
        raise FormalError('Template requires distinct explicit roles: '+','.join(names))
    if not event.evidence:raise Unsupported('Template requires source and assumption evidence')
    expected_kind=('COORDINATION' if name in {'LEWIS_ASSOCIATION','LEWIS_DISSOCIATION'} else
        'RADICAL' if name in {'RADICAL_DENITROGENATION','HAT_AROMATIZATION','SN_HAT'} else
        'SET' if name=='CU_DIAZONIUM_SET' else
        'PERICYCLIC' if name in {'AMMONIUM_YLIDE_2_3','SELENOXIDE_ELIMINATION'} else 'POLAR')
    if event.event_type!=expected_kind:raise FormalError('Template/event classification mismatch: '+expected_kind)
    ids=[r[k] for k in names];after=None
    def bond(a,b,n=1):
        if order(view,a,b)!=n:raise FormalError(f'{name} requires bond {a}-{b} order {n}')
    def el(i,*symbols):
        if atoms[i].element not in symbols:raise FormalError(f'{name} role atom {i} must be {symbols}')
    def rad(i):
        if atoms[i].radical_electrons!=1:raise FormalError('Specified singly occupied center required')
    def metal(i,element=None,os=None):
        if i not in {m.atom_id for m in view.metal_centers}:raise FormalError('Explicit MetalCenter required')
        if element:el(i,element)
        if os is not None:
            from .elements import coordination
            if next(c for c in coordination(view)[1] if c['metal_atom_id']==i)['oxidation_state']!=os:raise Unsupported('Metal oxidation state outside template scope')
    if name.startswith('LEWIS_'):
        m=r['metal'];metal(m)
        allowed={'Al':3,'Fe':3,'Zn':2,'Hg':2,'Pd':2}
        if atoms[m].element not in allowed:raise Unsupported('Lewis acid family not registered for this element')
        metal(m,os=allowed[atoms[m].element])
        if name in {'LEWIS_ASSOCIATION','LEWIS_DISSOCIATION'}:
            d=r['donor'];el(d,'N','O','P','S','Se','Br','Cl');bond(d,m,0)
            old=next((b for b in view.metal_bindings if b.metal_atom_id==m and b.donor_atom_ids==[d] and b.mode=='dative'),None)
            if name.endswith('ASSOCIATION'):
                if old:raise FormalError('Donor is already coordinated')
                limit={'Al':4,'Fe':6,'Zn':4,'Hg':4,'Pd':4}[atoms[m].element]
                if sum(b.metal_atom_id==m for b in view.metal_bindings)>=limit:raise Unsupported('Registered coordination capacity reached')
                def after(v):
                    v.metal_bindings.append(MetalBinding(binding_id=f'new.{event.event_id}',metal_atom_id=m,ligand_id='temporary',donor_atom_ids=[d],mode='dative'))
                    rebuild_bindings(v)
            else:
                if not old:raise FormalError('Specified coordination is absent')
                def after(v):
                    v.metal_bindings=[b for b in v.metal_bindings if b.binding_id!=old.binding_id];rebuild_bindings(v)
            actions=[]
        elif name=='LEWIS_ACYL_ACTIVATION':
            c,x,m,d=ids;el(c,'C');el(d,'O');el(x,'Cl','Br');bond(c,d,2);bond(c,x);el(m,'Al')
            if not any(b.metal_atom_id==m and b.donor_atom_ids==[d] and b.mode=='dative' for b in view.metal_bindings):raise FormalError('Acyl oxygen must first be coordinated')
            actions=[Shift(a=c,b=x,c=m)]
            def after(v):
                v.metal_bindings=[b for b in v.metal_bindings if not(b.metal_atom_id==m and b.donor_atom_ids==[d])];rebuild_bindings(v)
        elif name=='LEWIS_BROMINATION':
            a,c,br,x,m=ids;el(a,'C');el(c,'C');el(br,'Br');el(x,'Br');el(m,'Fe');bond(a,c,2);bond(br,x)
            if not any(b.metal_atom_id==m and b.donor_atom_ids==[x] and b.mode=='dative' for b in view.metal_bindings):raise FormalError('Br2 donor must first be coordinated')
            actions=[Shift(a=a,b=c,c=br),Shift(a=br,b=x,c=m)]
            def after(v):
                v.metal_bindings=[b for b in v.metal_bindings if not(b.metal_atom_id==m and b.donor_atom_ids==[x])];rebuild_bindings(v)
        else:
            m,x,h,c,a=ids;el(h,'H');el(x,'Cl','Br');bond(m,x);bond(h,c);bond(c,a)
            actions=[Pair(a=x,b=h),Shift(a=h,b=c,c=a),Split(a=m,b=x,target=x)];after=rebuild_bindings
    elif name=='AMMONIUM_YLIDE_2_3':
        y,n,a,b,c=ids;el(y,'C');el(n,'N')
        for i in (a,b,c):el(i,'C')
        for i,j,k in [(y,n,1),(n,a,1),(a,b,1),(b,c,2),(y,c,0)]:bond(i,j,k)
        if (atoms[y].formal_charge,atoms[n].formal_charge)!=(-1,1) or any(atoms[i].radical_electrons for i in ids):raise FormalError('Closed-shell ammonium ylide occupation required')
        if event.activation=='photochemical':raise Unsupported('Only thermal ylide [2,3] profile registered')
        actions=[Pair(a=y,b=c),Shift(a=c,b=b,c=a),Split(a=a,b=n,target=n)]
    elif name=='RADICAL_DENITROGENATION':
        n2,n1,c=ids;el(n2,'N');el(n1,'N');el(c,'C');rad(n2);bond(n2,n1,2);bond(n1,c)
        from .substitution import aromatic_ring
        if not aromatic_ring(view,c):raise Unsupported('Denitrogenation profile requires an aryl carbon')
        if bond_valences(view)[n2]!=2 or bond_valences(view)[n1]!=3:raise Unsupported('Unsubstituted terminal diazenyl N required')
        actions=[atom_to_bond(n2,n2,n1,'pi2'),bond_to_bond(n1,c,'sigma',0,n2,n1,'pi2'),bond_to_atom(n1,c,'sigma',1,c)]
    elif name=='SN_HAT':
        c,h,sn=ids;el(c,'C');el(h,'H');el(sn,'Sn');rad(c);bond(sn,h);bond(c,h,0)
        if bond_valences(view)[sn]!=4:raise Unsupported('Tin donor must be tetravalent R3SnH')
        others=[atoms[b.b if b.a==sn else b.a].element for b in view.bonds if sn in {b.a,b.b} and h not in {b.a,b.b}]
        if sorted(others)!=['C','C','C']:raise Unsupported('Sn HAT requires three explicit carbon substituents')
        actions=[atom_to_bond(c,c,h),bond_to_bond(sn,h,'sigma',0,c,h),bond_to_atom(sn,h,'sigma',1,sn)]
    elif name=='HAT_AROMATIZATION':
        r0,h,a,b,c,d=ids;rad(r0);rad(a);el(h,'H')
        for i in (a,b,c,d):el(i,'C')
        for i,j,k in [(a,b,1),(b,c,2),(c,d,1),(d,h,1),(r0,h,0)]:bond(i,j,k)
        remainder=set(atoms)-set(ids)
        if not any(atoms[e].element==atoms[f].element=='C' and order(view,d,e)==1 and
                   order(view,e,f)==2 and order(view,f,a)==1 for e in remainder for f in remainder if e!=f):
            raise Unsupported('HAT/aromatization requires the specified cyclohexadienyl six-membered ring')
        actions=[atom_to_bond(r0,r0,h),bond_to_bond(d,h,'sigma',0,r0,h),
                 bond_to_bond(d,h,'sigma',1,c,d,'pi1'),bond_to_bond(b,c,'pi1',0,c,d,'pi1'),
                 bond_to_bond(b,c,'pi1',1,a,b,'pi1'),atom_to_bond(a,a,b,'pi1')]
    elif name=='CU_DIAZONIUM_SET':
        cu,c,n,t=ids;el(cu,'Cu');el(c,'C');el(n,'N');el(t,'N');bond(c,n);bond(n,t,3)
        from .substitution import aromatic_ring
        if not aromatic_ring(view,c) or bond_valences(view)[t]!=3:raise Unsupported('Only aryl diazonium SET registered')
        if (atoms[n].formal_charge,atoms[n].nonbonding_electrons)!=(1,0):raise FormalError('Aryl diazonium N+ input required')
        if (atoms[cu].nonbonding_electrons,atoms[cu].radical_electrons,atoms[cu].formal_charge) not in {(11,1,0),(10,0,1)} or bond_valences(view)[cu]:
            raise Unsupported('SET scope: explicitly declared atomic Cu0 or CuI ion; no guessed spin/ligands')
        from ..v02.ast import RadicalSource
        src=RadicalSource(atom_id=cu) if atoms[cu].radical_electrons==1 else LonePairSource(atom_id=cu,pair_index=0,slot=0)
        actions=[Fishhook(source=src,destination=LonePairDestination(atom_id=n)),
                 Fishhook(source={'kind':'bond','a':c,'b':n,'pair':'sigma','slot':0},destination=LonePairDestination(atom_id=n)),
                 bond_to_atom(c,n,'sigma',1,c)]
        after=rebuild_bindings
    elif name=='SAEGUSA_PALLADATION':
        m,a,c,o,si,x=ids;metal(m,'Pd',2);el(si,'Si');el(o,'O');el(x,'O');el(a,'C');el(c,'C')
        for i,j,k in [(a,c,2),(c,o,1),(o,si,1),(m,x,1)]:bond(i,j,k)
        actions=[Split(a=si,b=o,target=o),Pair(a=o,b=c),Shift(a=c,b=a,c=m),Split(a=m,b=x,target=x),Pair(a=x,b=si)]
        after=rebuild_bindings
    elif name=='PD_BETA_H_ELIMINATION':
        m,a,b,h=ids;metal(m,'Pd',2);el(a,'C');el(b,'C');el(h,'H');bond(m,a);bond(a,b);bond(b,h);bond(m,h,0)
        if event.conditions.get('geometry_status')=='incompatible' or event.conditions.get('vacant_site_available') is False:
            raise FormalError('Explicit geometry/vacancy contradicts beta-H template')
        if event.conditions.get('product_alkene_binding','free')!='free':raise Unsupported('This beta-H profile exports free alkene; eta-bound variant not registered')
        actions=[Shift(a=m,b=a,c=b),Shift(a=b,b=h,c=m)];after=rebuild_bindings
    elif name=='PD_H_ACETATE_ELIMINATION':
        m,h,o=ids;metal(m,'Pd',2);el(h,'H');el(o,'O');bond(m,h);bond(m,o)
        actions=[Shift(a=m,b=h,c=o),Split(a=o,b=m,target=m)];after=rebuild_bindings
    elif name=='P_SE_ACTIVATION':
        p,se,c,n=ids;el(p,'P');el(se,'Se');el(c,'C');el(n,'N');bond(se,c);bond(c,n,3);bond(p,se,0)
        if (bond_valences(view)[p],atoms[p].nonbonding_electrons)!=(3,2):raise FormalError('Trivalent phosphine donor required')
        actions=[Pair(a=p,b=se),Split(a=se,b=c,target=c)]
    else:
        se,o,a,b,h=ids;el(se,'Se');el(o,'O');el(a,'C');el(b,'C');el(h,'H')
        for i,j,k in [(se,o,2),(se,a,1),(a,b,1),(b,h,1)]:bond(i,j,k)
        if event.conditions.get('geometry_status')=='incompatible':raise FormalError('Explicit syn geometry contradicted')
        actions=[Shift(a=se,b=o,c=h),Shift(a=h,b=b,c=a),Split(a=a,b=se,target=se)]
    report={'template_id':name,'scope':scope,'source':source,'geometry_status':'unverified',
            'mechanism_evidence':'assumed_by_template','chemical_feasibility':'unknown','atomic_commit':True}
    return actions,after,report
