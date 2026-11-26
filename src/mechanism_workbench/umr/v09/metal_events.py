"""Finite topology and ownership recipes, with complete product ligand partitions."""
from ..ast import Pair,Split,Shift
from ..chemistry import FormalError,Unsupported,atoms_by_id,bond_orders,order
from .ast import Ligand,MetalBinding
from .elements import current_elements
from .coordination import checks
SPECS={
 'ZR_CP2_H_ALKENE_ADDITION':['metal','hydrogen','a','b'],
 'PD_ASSOCIATIVE_ALKYNE_INSERTION':['metal','migrating','a','b'],
 'RH_ACYL_H_OXIDATIVE_ADDITION':['metal','carbonyl','oxygen','hydrogen'],
 'RH_ACYL_DECARBONYLATION':['metal','carbonyl','oxygen','organic','hydrogen'],
 'RH_C_H_REDUCTIVE_ELIMINATION':['metal','organic','hydrogen'],
 'TI_ENONE_IONIZING_ACTIVATION':['metal','oxygen','carbonyl','alpha','beta','halogen'],
 'TI_ENOLATE_VINYL_CLOSURE':['metal','oxygen','carbonyl','alpha','vinyl','halogen'],
 'AL_CYANIDE_ASSOCIATION':['metal','carbon','nitrogen'],
 'AL_CYANIDE_CONJUGATE_TRANSFER':['metal','cyanide_c','nitrogen','beta','alpha','carbonyl','oxygen','hydrogen'],
 'ACYL_CS_OXIDATIVE_ADDITION':['metal','carbon','sulfur','oxygen'],
 'METAL_RX_INSERTION':['metal','organic','halogen'],
 'METAL_CARBONYL_ADDITION':['metal','nucleophile','carbon','oxygen'],
 'METAL_OXO_TRANSFER':['metal','oxygen','a','b'],
 'ALKENE_MIGRATORY_INSERTION':['metal','migrating','a','b'],
 'CO_MIGRATORY_INSERTION':['metal','migrating','carbon','oxygen'],
 'METALLACYCLE_FORMATION':['metal','alkylidene','a','b'],
 'METALLACYCLE_CLEAVAGE':['metal','first','second','third'],
 'METAL_SOLVENT_X_EXCHANGE':['metal','leaving','donor'],
 'AQUEOUS_METAL_QUENCH':['metal','organic','water_o','hydrogen'],
 'AMINE_HYDRIDE_TRANSFER':['metal','nitrogen','carbon','hydrogen'],
 'AMINOPALLADATION':['metal','nitrogen','a','b'],
 'ORGANO_MG_NITRO_ADDITION':['metal','organic','nitro_o','nitrogen','second_o'],
 'IONIC_C_LI_ASSOCIATION':['metal','carbon'],
 'ANION_METAL_ASSOCIATION':['metal','oxygen'],
 'O_DONOR_X_SUBSTITUTION':['metal','donor','leaving'],
 'PB_ALCOHOL_ACETATE_EXCHANGE':['metal','oxygen','leaving','hydrogen'],
 'AL_ALKOXIDE_HYDRIDE_TRANSFER':['metal','donor_o','donor_c','hydrogen','acceptor_c','acceptor_o'],
 'ZN_CARBENOID_TRANSFER':['metal','carbenoid','halogen','a','b'],
 'TI_ALPHA_H_ELIMINATION':['metal','alkylidene','departing','hydrogen'],
 'TI_OXAMETALLACYCLE':['metal','alkylidene','carbon','oxygen'],
 'HG_ENOL_ETHER_CAPTURE':['metal','a','b','oxygen'],
 'ASSOCIATIVE_ALKENE_INSERTION':['metal','migrating','a','b'],
}

def repartition(v):
    """No inference of haptic/bridge electronic assignments after graph edits."""
    ats=atoms_by_id(v);mids={i for i,a in ats.items() if a.element in current_elements().metals};edges=bond_orders(v)
    if any(l.internal_metal_ids for l in v.ligands):raise Unsupported('Internal-metal ligand repartition needs explicit dedicated model')
    # Untouched Cp/allyl regions remain authoritative haptic objects. The
    # transaction guard already prohibits unregistered edits to their donors.
    frozen={b.ligand_id for b in v.metal_bindings if b.mode=='eta'}
    ls=[l for l in v.ligands if l.ligand_id in frozen];bs=[b for b in v.metal_bindings if b.ligand_id in frozen]
    fixed={i for l in ls for i in l.atom_ids}
    ds=[(i,b.metal_atom_id) for b in v.metal_bindings if b.mode=='dative' and b.ligand_id not in frozen for i in b.donor_atom_ids]
    old_ligands=v.ligands;todo=set(ats)-mids-fixed
    while todo:
        start=todo.pop();ids={start};stack=[start]
        while stack:
            i=stack.pop()
            for a,b in edges:
                if i not in {a,b}:continue
                j=b if a==i else a
                if j in todo:todo.remove(j);ids.add(j);stack.append(j)
        contacts=[(i,m,'dative',0) for i,m in ds if i in ids]
        contacts += [(i,m,'sigma_covalent',o) for (a,b),o in edges.items() for i,m in [(a,b),(b,a)] if i in ids and m in mids]
        if not contacts:
            # A declared incoming pi ligand may be used by a later event.
            # Preserve it only while its full component is still unchanged
            # and it was unbound before this transaction.
            old=next((l for l in old_ligands if set(l.atom_ids)==ids and not any(b.ligand_id==l.ligand_id for b in v.metal_bindings)),None)
            if old is not None:ls.append(old.model_copy(deep=True))
            continue
        if len({m for i,m,mode,o in contacts})>1:raise Unsupported('Repartition of an active bridging ligand needs explicit allocations')
        q=sum(ats[i].formal_charge for i in ids)-sum(o for i,m,mode,o in contacts)
        neutral=sum(mode=='dative' for _,_,mode,_ in contacts)
        cls={(0,0):'L',(0,1):'L',(0,2):'L2',(-1,0):'X',(-1,1):'LX',(-1,2):'L2X',(-2,0):'X2',(-2,1):'LX2',(-2,2):'L2X2',(-3,0):'X3'}.get((q,neutral))
        if cls is None:raise Unsupported('Product ligand ionic class not registered')
        lid='ligand.'+str(min(ids));donors=sorted({i for i,_,_,_ in contacts})
        previous=next((l for l in old_ligands if set(l.atom_ids)==ids),None)
        typ=previous.ligand_type if previous else 'generic'
        neutral_carbene=bool(previous and previous.ionic_class=='L' and previous.ligand_type in {'Ru_alkylidene','Fischer_carbene'}
            and len(contacts)==1 and contacts[0][2:] == ('sigma_covalent',2) and previous.donor_atom_ids==donors)
        if neutral_carbene:cls='L';q=0
        if cls=='X3' and all(ats[m].element in {'W','Mo'} for _,m,_,_ in contacts):
            typ='Schrock_alkylidyne' if len(contacts)==1 and contacts[0][3]==3 else 'metallacyclobutadiene'
        ls.append(Ligand(ligand_id=lid,atom_ids=sorted(ids),donor_atom_ids=donors,ionic_class=cls,ionic_charge=q,
          ligand_type=typ,assignment_basis='v0.9 registered event ionic partition'))
        for k,(i,m,mode,o) in enumerate(contacts):bs.append(MetalBinding(binding_id=f'{lid}.{m}.{k}',ligand_id=lid,metal_atom_id=m,donor_atom_ids=[i],mode=mode,donated_electrons=2 if mode=='dative' or neutral_carbene else 2*o))
    v.ligands=ls;v.metal_bindings=bs
    for m in v.metal_centers:
        m.oxidation_state=None;m.d_electrons=None;m.valence_electrons=None;m.relative_positions=[];m.vacant_site_ids=None

def compile(v,e):
    names=SPECS[e.general_event];r=e.general_roles;ats=atoms_by_id(v)
    if set(r)!=set(names) or len(set(r.values()))!=len(r) or not set(r.values())<=ats.keys():raise FormalError('Exact distinct metal event roles required')
    if not e.evidence or not e.mechanism_assumption:raise FormalError('Metal event needs evidence and a declared channel')
    ids=[r[k] for k in names];m=ids[0];counts={x['metal_atom_id']:x for x in checks(v)[1]}
    if m not in counts:raise FormalError('Missing metal declaration')
    el=ats[m].element;os=counts[m]['oxidation_state']
    def bond(a,b,o):
        if order(v,a,b)!=o:raise FormalError(f'Expected bond {a}-{b} order {o}')
    if e.general_event=='ZR_CP2_H_ALKENE_ADDITION':
        m,h,a,b=ids
        if (el,os)!=('Zr',4) or [ats[i].element for i in ids[1:]]!=['H','C','C']:raise Unsupported('ZrIV hydride/alkene scope')
        ligs={l.ligand_id:l for l in v.ligands}
        cps=[x for x in v.metal_bindings if x.metal_atom_id==m and x.mode=='eta' and x.hapticity==5 and ligs[x.ligand_id].ligand_type=='Cp']
        if len(cps)!=2 or any(ligs[x.ligand_id].ionic_class!='L2X' or x.donated_electrons!=6 for x in cps):raise Unsupported('Two explicitly assigned Cp(-) six-electron donors required')
        if not any(ats[j].element=='Cl' and order(v,m,j)==1 for j in ats):raise Unsupported('Chlorozirconocene hydride only')
        for i,j,n in [(m,h,1),(a,b,2),(h,a,0),(b,m,0)]:bond(i,j,n)
        if any(ats[i].formal_charge or ats[i].radical_electrons for i in [a,b,h]):raise FormalError('Neutral paired hydride/alkene required')
        return [Shift(a=m,b=h,c=a),Shift(a=a,b=b,c=m)]
    if e.general_event=='PD_ASSOCIATIVE_ALKYNE_INSERTION':
        m,r,a,b=ids
        if (el,os)!=('Pd',2) or any(ats[i].element!='C' for i in ids[1:]):raise Unsupported('PdII sigma-carbon/alkyne input only')
        if e.conditions.get('association_model')!='unresolved_within_event':raise Unsupported('Explicit unisolated association declaration required')
        for i,j,n in [(m,r,1),(a,b,3),(r,a,0),(b,m,0)]:bond(i,j,n)
        if not any(set(l.donor_atom_ids)=={a,b} and l.ionic_class=='L' for l in v.ligands):raise FormalError('Full incoming pi-ligand declaration required')
        return [Shift(a=m,b=r,c=a),Shift(a=a,b=b,c=m)]
    if e.general_event.startswith('RH_'):
        if el!='Rh':raise Unsupported('Rhodium-only acyl/hydride profile')
        if e.general_event=='RH_ACYL_H_OXIDATIVE_ADDITION':
            m,c,o,h=ids
            if os!=1 or [ats[i].element for i in ids[1:]]!=['C','O','H']:raise Unsupported('RhI aldehyde activation only')
            if ats[m].radical_electrons!=0 or ats[m].nonbonding_electrons<2:raise Unsupported('Explicit paired RhI donor occupation required; no automatic low-spin assignment')
            for i,j,n in [(c,o,2),(c,h,1),(m,c,0),(m,h,0)]:bond(i,j,n)
            return [Pair(a=m,b=c),Shift(a=c,b=h,c=m)]
        if os!=3:raise Unsupported('RhIII intermediate required')
        if e.general_event=='RH_ACYL_DECARBONYLATION':
            m,c,o,r,h=ids
            if [ats[i].element for i in ids[1:]]!=['C','O','C','H']:raise Unsupported('RhIII acyl/hydride topology only')
            for i,j,n in [(m,c,1),(c,o,2),(c,r,1),(m,h,1),(m,r,0)]:bond(i,j,n)
            return [Shift(a=c,b=r,c=m),Split(a=m,b=c,target=c),Pair(a=o,b=c)]
        m,c,h=ids
        if [ats[i].element for i in ids[1:]]!=['C','H']:raise Unsupported('RhIII C/H reductive elimination only')
        for i,j,n in [(m,c,1),(m,h,1),(c,h,0)]:bond(i,j,n)
        return [Split(a=c,b=m,target=m),Shift(a=m,b=h,c=c)]
    if e.general_event in {'TI_ENONE_IONIZING_ACTIVATION','TI_ENOLATE_VINYL_CLOSURE'}:
        m,o,k,a,b,x=ids
        if (el,os)!=('Ti',4) or [ats[i].element for i in ids[1:]]!=['O','C','C','C','Cl']:raise Unsupported('TiIV/O/C/C/C/Cl finite enone scope')
        if e.general_event=='TI_ENONE_IONIZING_ACTIVATION':
            for i,j,n in [(m,x,1),(o,k,2),(k,a,1),(a,b,2),(m,o,0)]:bond(i,j,n)
            if any(ats[i].formal_charge or ats[i].radical_electrons for i in ids[1:]):raise FormalError('Neutral paired enone/chloride input required')
            return [Pair(a=o,b=m),Split(a=m,b=x,target=x),Split(a=k,b=o,target=o),Shift(a=b,b=a,c=k)]
        for i,j,n in [(m,o,1),(o,k,1),(k,a,2),(a,b,0),(m,x,0)]:bond(i,j,n)
        if (ats[b].formal_charge,ats[b].nonbonding_electrons,ats[b].radical_electrons)!=(1,0,0) or ats[x].formal_charge!=-1:raise FormalError('Declared vinyl cation and retained chloride required')
        if not any(order(v,b,j)==2 and ats[j].element=='C' for j in ats):raise Unsupported('Vinyl cation scope only')
        return [Pair(a=o,b=k),Shift(a=k,b=a,c=b),Split(a=m,b=o,target=o),Pair(a=x,b=m)]
    if e.general_event=='AL_CYANIDE_ASSOCIATION':
        m,c,n=ids
        if (el,os)!=('Al',3) or [ats[i].element for i in [c,n]]!=['C','N']:raise Unsupported('AlIII/cyanide scope')
        if (ats[c].formal_charge,ats[c].nonbonding_electrons,ats[c].radical_electrons)!=(-1,2,0):raise FormalError('Cyanide paired carbon anion required')
        bond(c,n,3);bond(m,c,0)
        return [Pair(a=c,b=m)]
    if e.general_event=='AL_CYANIDE_CONJUGATE_TRANSFER':
        m,c,n,b,a,k,o,h=ids
        if (el,os)!=('Al',3) or [ats[i].element for i in ids[1:]]!=['C','N','C','C','C','O','H']:raise Unsupported('AlIII cyanide/protonated enone scope')
        for i,j,z in [(m,c,1),(c,n,3),(b,a,2),(a,k,1),(k,o,2),(o,h,1),(c,b,0)]:bond(i,j,z)
        if ats[m].formal_charge!=-1 or ats[o].formal_charge!=1 or any(ats[i].radical_electrons for i in ids[1:]):raise FormalError('Explicit cyanidoaluminate/protonated paired enone required')
        return [Shift(a=m,b=c,c=b),Shift(a=b,b=a,c=k),Split(a=k,b=o,target=o)]
    if e.general_event=='ACYL_CS_OXIDATIVE_ADDITION':
        m,c,s,o=ids
        if (el,os)!=('Pd',0) or [ats[i].element for i in [c,s,o]]!=['C','S','O']:
            raise Unsupported('Finite mononuclear Pd0 thioester acyl C-S oxidative addition')
        if any(ats[i].formal_charge or ats[i].radical_electrons for i in [c,s,o]):
            raise Unsupported('Neutral closed-shell thioester required')
        for i,j,n in [(c,s,1),(c,o,2),(m,c,0),(m,s,0)]:bond(i,j,n)
        if ats[s].hydrogens or sum(order(v,s,i) for i in ats)!=2:
            raise Unsupported('Thioester S must have two single covalent contacts')
        return [Pair(a=m,b=c),Shift(a=c,b=s,c=m)]
    if e.general_event=='ANION_METAL_ASSOCIATION':
        m,o=ids
        if (el,os)!=('Al',3) or ats[o].element!='O' or (ats[o].formal_charge,ats[o].radical_electrons)!=(-1,0):raise Unsupported('Finite AlIII/alkoxide covalent Lewis association')
        bond(m,o,0);return [Pair(a=o,b=m)]
    if e.general_event=='O_DONOR_X_SUBSTITUTION':
        m,o,x=ids
        if (el,os) not in {('Cr',6),('Ti',4)} or ats[o].element!='O' or ats[x].element not in {'O','Cl'}:raise Unsupported('Finite CrVI/TiIV oxygen donor substitution')
        bond(m,x,1);bond(m,o,0)
        return [Pair(a=o,b=m),Split(a=m,b=x,target=x)]
    if e.general_event=='PB_ALCOHOL_ACETATE_EXCHANGE':
        m,o,x,h=ids
        if (el,os)!=('Pb',4) or [ats[i].element for i in [o,x,h]]!=['O','O','H']:raise Unsupported('Finite PbIV alcohol/acetate ligand exchange')
        for a,b,n in [(m,x,1),(o,h,1),(m,o,0),(x,h,0)]:bond(a,b,n)
        return [Pair(a=o,b=m),Split(a=m,b=x,target=x),Pair(a=x,b=h),Split(a=h,b=o,target=o)]
    if e.general_event=='AL_ALKOXIDE_HYDRIDE_TRANSFER':
        m,o,c,h,a,x=ids
        if (el,os)!=('Al',3) or [ats[i].element for i in [o,c,h,a,x]]!=['O','C','H','C','O']:raise Unsupported('Finite monomeric AlIII six-center alkoxide/carbonyl hydride transfer')
        for i,j,n in [(m,o,1),(o,c,1),(c,h,1),(a,x,2),(h,a,0),(x,m,0)]:bond(i,j,n)
        if not released_bindings(v,e):raise FormalError('Carbonyl oxygen must already be donor-bound to Al')
        return [Shift(a=c,b=h,c=a),Split(a=a,b=x,target=x),Pair(a=o,b=c),Split(a=m,b=o,target=o),Pair(a=x,b=m)]
    if e.general_event=='ZN_CARBENOID_TRANSFER':
        m,c,x,a,b=ids
        if (el,os)!=('Zn',2) or [ats[i].element for i in [c,x,a,b]] not in [['C','I','C','C'],['C','Br','C','C']]:raise Unsupported('Finite ZnII halomethyl/alkene transfer')
        for i,j,n in [(m,c,1),(c,x,1),(a,b,2),(c,a,0),(c,b,0),(m,x,0)]:bond(i,j,n)
        return [Shift(a=m,b=c,c=a),Shift(a=a,b=b,c=c),Shift(a=c,b=x,c=m)]
    if e.general_event=='TI_ALPHA_H_ELIMINATION':
        m,c,r,h=ids
        if (el,os)!=('Ti',4) or [ats[i].element for i in [c,r,h]]!=['C','C','H']:raise Unsupported('Finite TiIV dialkyl alpha-H elimination')
        for i,j,n in [(m,r,1),(m,c,1),(c,h,1),(r,h,0)]:bond(i,j,n)
        return [Shift(a=m,b=r,c=h),Shift(a=h,b=c,c=m)]
    if e.general_event=='TI_OXAMETALLACYCLE':
        m,c,a,o=ids
        if (el,os)!=('Ti',4) or [ats[i].element for i in [c,a,o]]!=['C','C','O']:raise Unsupported('Finite TiIV alkylidene/carbonyl metallacycle; not ordinary organic [2+2]')
        if e.event_type!='POLAR':raise FormalError('Metal-assisted event classified POLAR, not ordinary pericyclic rule')
        for i,j,n in [(m,c,2),(a,o,2),(c,a,0),(o,m,0)]:bond(i,j,n)
        return [Shift(a=m,b=c,c=a),Shift(a=a,b=o,c=m)]
    if e.general_event=='HG_ENOL_ETHER_CAPTURE':
        m,a,b,o=ids
        if (el,os)!=('Hg',2) or [ats[i].element for i in [a,b,o]]!=['C','C','O']:raise Unsupported('Finite HgII capture of enol ether')
        for i,j,n in [(a,b,2),(o,a,1),(b,m,0)]:bond(i,j,n)
        return [Shift(a=a,b=b,c=m),Pair(a=o,b=a)]
    if e.general_event=='IONIC_C_LI_ASSOCIATION':
        m,c=ids
        if (el,os)!=('Li',1) or ats[c].element!='C' or ats[c].formal_charge!=-1 or ats[c].radical_electrons!=0:raise Unsupported('Finite ionic Li+/carbanion association to covalent Lewis convention')
        bond(m,c,0)
        return [Pair(a=c,b=m)]
    if e.general_event=='ORGANO_MG_NITRO_ADDITION':
        m,c,o,n,x=ids
        if (el,os)!=('Mg',2) or [ats[i].element for i in [c,o,n,x]]!=['C','O','N','O'] or (ats[n].formal_charge,ats[x].formal_charge)!=(1,-1):raise Unsupported('Explicit nitro ionic topology and monomeric organomagnesium required')
        for a,b,z in [(m,c,1),(o,n,2),(n,x,1),(c,o,0),(x,m,0)]:bond(a,b,z)
        return [Shift(a=m,b=c,c=o),Split(a=o,b=n,target=n),Pair(a=x,b=m)]
    if e.general_event=='AMINOPALLADATION':
        m,n,a,b=ids
        if (el,os)!=('Pd',2) or [ats[i].element for i in [n,a,b]]!=['N','C','C']:raise Unsupported('Finite PdII/N-donor alkene addition')
        if ats[n].formal_charge!=0 or ats[n].radical_electrons!=0:raise Unsupported('Neutral closed-shell amine nucleophile required')
        for x,y,z in [(a,b,2),(n,a,0),(b,m,0)]:bond(x,y,z)
        if not released_bindings(v,e):raise FormalError('Aminopalladation requires explicit eta2 substrate coordination')
        return [Pair(a=n,b=a),Shift(a=a,b=b,c=m)]
    if e.general_event=='AMINE_HYDRIDE_TRANSFER':
        m,n,c,h=ids
        if (el,os)!=('Pd',2) or [ats[i].element for i in [n,c,h]]!=['N','C','H']:raise Unsupported('Finite PdII amine alpha-H transfer scope')
        if ats[n].formal_charge!=0 or ats[n].nonbonding_electrons!=2 or ats[n].radical_electrons!=0:raise FormalError('Neutral closed-shell coordinated amine required')
        for a,b,o in [(n,c,1),(c,h,1),(m,h,0)]:bond(a,b,o)
        if not released_bindings(v,e):raise FormalError('Amine transfer requires explicit donor-owned N coordination')
        return [Pair(a=n,b=c),Shift(a=c,b=h,c=m)]
    if e.general_event=='AQUEOUS_METAL_QUENCH':
        m,x,o,h=ids
        if (el,os) not in {('Mg',2),('Zn',2)} or ats[x].element not in {'C','O'}:
            raise Unsupported('Aqueous quench limited to monomeric Mg/Zn C/O sigma ligands')
        if ats[o].element!='O' or ats[h].element!='H' or ats[o].formal_charge or ats[o].hydrogens:
            raise FormalError('Quench water and transferred H must be explicit')
        hs=[i for edge in bond_orders(v) if o in edge for i in edge if i!=o]
        if len(hs)!=2 or any(ats[i].element!='H' for i in hs):raise FormalError('Specified donor must be water')
        for a,b,n in [(m,x,1),(o,h,1),(x,h,0),(o,m,0)]:bond(a,b,n)
        if e.conditions.get('quench_medium')!='water':raise FormalError('Explicit aqueous assumption required')
        # Both arrows read the same input. M-X electrons form X-H, while
        # O-H electrons form O-M. No isolated proton or missing hydroxide.
        return [Shift(a=m,b=x,c=h),Shift(a=h,b=o,c=m)]
    if e.general_event=='METAL_SOLVENT_X_EXCHANGE':
        m,x,d=ids
        if (el,os) not in {('Mg',2),('Zn',2)} or ats[d].element!='O' or ats[x].element not in {'O','Cl','Br','I'}:raise Unsupported('Finite Mg/Zn oxygen-donor X replacement')
        bond(m,x,1);bond(d,m,0)
        if ats[d].formal_charge or ats[d].nonbonding_electrons!=4:raise FormalError('Neutral closed-shell oxygen donor required')
        return [Pair(a=d,b=m),Split(a=m,b=x,target=x)]
    if e.general_event in {'METALLACYCLE_FORMATION','METALLACYCLE_CLEAVAGE'}:
        m,a,b,c=ids
        if el not in {'Ru','Mo'} or any(ats[i].element!='C' for i in [a,b,c]):raise Unsupported('Finite Ru/Mo alkylidene metallacycle profile')
        if e.event_type!='POLAR':raise FormalError('Metal-assisted metallacycle is distinct from ordinary pericyclic templates')
        if e.general_event=='METALLACYCLE_FORMATION':
            # a = original alkylidene; b,c = specified alkene ends.
            for x,y,n in [(m,a,2),(b,c,2),(a,b,0),(c,m,0)]:bond(x,y,n)
            ls={l.ligand_id:l for l in v.ligands}
            lig=next((ls[x.ligand_id] for x in v.metal_bindings if x.metal_atom_id==m and x.donor_atom_ids==[a] and x.mode=='sigma_covalent'),None)
            if lig is None or lig.ionic_class!='X2' or lig.ligand_type not in {'Ru_alkylidene','Schrock_alkylidene'}:raise Unsupported('Explicit X2 alkylidene convention required; neutral-L representation is not silently converted')
            return [Shift(a=m,b=a,c=b),Shift(a=b,b=c,c=m)]
        for x,y in [(m,a),(a,b),(b,c),(c,m)]:bond(x,y,1)
        return [Shift(a=a,b=m,c=c),Shift(a=c,b=b,c=a)]
    if e.general_event in {'ALKENE_MIGRATORY_INSERTION','ASSOCIATIVE_ALKENE_INSERTION'}:
        m,r,a,b=ids
        if (el,os) not in {('Pd',2),('Ni',2),('Rh',1),('Ru',2)} or ats[r].element not in {'C','H'} or any(ats[i].element!='C' for i in [a,b]):raise Unsupported('Finite cis/vacancy-dependent migratory insertion scope')
        for x,y,n in [(m,r,1),(a,b,2),(r,a,0),(b,m,0)]:bond(x,y,n)
        if e.general_event=='ALKENE_MIGRATORY_INSERTION' and not released_bindings(v,e):raise FormalError('Insertion requires explicit preceding eta2 association')
        if e.general_event=='ASSOCIATIVE_ALKENE_INSERTION':
            if (el,os)!=('Pd',2) or e.conditions.get('association_model')!='unresolved_within_event':raise Unsupported('Associative PdII insertion requires explicit unisolated-association assumption')
            if not any(set(l.donor_atom_ids)=={a,b} and l.ionic_class=='L' for l in v.ligands):raise FormalError('Declared full incoming pi-ligand required')
        return [Shift(a=m,b=r,c=a),Shift(a=a,b=b,c=m)]
    if e.general_event=='CO_MIGRATORY_INSERTION':
        m,r,c,o=ids
        if el not in {'Pd','Ni','Rh','Co','Cr','Ru'} or [ats[i].element for i in [r,c,o]]!=['C','C','O']:raise Unsupported('Finite organometallic carbonyl insertion scope')
        for x,y,n in [(m,r,1),(c,o,3),(r,c,0),(m,c,0)]:bond(x,y,n)
        if not released_bindings(v,e):raise FormalError('CO must already be explicitly C-bound')
        return [Shift(a=m,b=r,c=c),Split(a=c,b=o,target=o),Pair(a=c,b=m)]
    if e.general_event=='METAL_RX_INSERTION':
        m,c,x=ids
        if (el,os)!=('Zn',0) or ats[c].element!='C' or ats[x].element not in {'Cl','Br','I'}:raise Unsupported('Finite RX insertion is Zn(0) molecular ownership model, not bulk surface kinetics')
        bond(c,x,1);bond(m,c,0);bond(m,x,0)
        return [Pair(a=m,b=c),Shift(a=c,b=x,c=m)]
    if e.general_event=='METAL_CARBONYL_ADDITION':
        m,r,c,o=ids
        if (el,os) not in {('Mg',2),('Zn',2)} or [ats[i].element for i in [r,c,o]]!=['C','C','O']:raise Unsupported('Finite organo-Mg/Zn carbonyl addition scope')
        for a,b,n in [(m,r,1),(c,o,2),(r,c,0),(o,m,0)]:bond(a,b,n)
        return [Shift(a=m,b=r,c=c),Split(a=c,b=o,target=o),Pair(a=o,b=m)]
    m,o,a,b=ids
    if (el,os) not in {('Mn',5),('Ru',6),('Mo',6)} or [ats[i].element for i in [o,a,b]]!=['O','C','C']:raise Unsupported('Finite terminal oxo epoxidation scope requires MnV/RuVI/MoVI')
    for x,y,n in [(m,o,2),(a,b,2),(o,a,0),(o,b,0)]:bond(x,y,n)
    return [Split(a=o,b=m,target=m),Shift(a=m,b=o,c=a),Shift(a=a,b=b,c=o)]

def released_bindings(v,e):
    r=e.general_roles;m=r.get('metal')
    if e.general_event in {'ALKENE_MIGRATORY_INSERTION','AMINOPALLADATION'}:return {b.binding_id for b in v.metal_bindings if b.metal_atom_id==m and b.mode=='eta' and set(b.donor_atom_ids)=={r['a'],r['b']} and b.hapticity==2}
    if e.general_event=='CO_MIGRATORY_INSERTION':return {b.binding_id for b in v.metal_bindings if b.metal_atom_id==m and b.mode=='dative' and b.donor_atom_ids==[r['carbon']]}
    if e.general_event=='AMINE_HYDRIDE_TRANSFER':return {b.binding_id for b in v.metal_bindings if b.metal_atom_id==m and b.mode=='dative' and b.donor_atom_ids==[r['nitrogen']]}
    if e.general_event=='AL_ALKOXIDE_HYDRIDE_TRANSFER':return {b.binding_id for b in v.metal_bindings if b.metal_atom_id==m and b.mode=='dative' and b.donor_atom_ids==[r['acceptor_o']]}
    return set()
