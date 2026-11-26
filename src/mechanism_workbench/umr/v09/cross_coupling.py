"""Source-conditioned molecular coupling channels; never a name-only bypass.

This optional r10 extension uses the existing atomic transaction and ionic
ligand partition. It does not assign a global spin, surface model, or rate.
"""
import json
from collections import Counter
from ..ast import Pair, Split, Shift
from ..chemistry import FormalError, Unsupported, atoms_by_id, order, bond_valences
from .coordination import checks

VERSION = 'umr-cross-coupling/1.0'
SPECS = {name: dict(roles=roles.split(), event_type='POLAR', rule_revision=VERSION,
    source=source, library='organometallic') for name, roles, source in [
 ('CC_OXIDATIVE_ADDITION', 'metal carbon leaving', 'Kurti/Czako: Pd0/Ni0 or CuI, declared molecular two-electron channel'),
 ('CC_TRANSMETALATION', 'metal leaving donor carbon', 'Kurti/Czako Stille/Negishi/Kumada/Sonogashira; explicit donor-ligand exchange'),
 ('CC_SUZUKI_OH_TRANSFER', 'metal oxygen boron carbon', 'doi:10.1021/ja1108326; hydroxo-Pd/neutral boronic partner'),
 ('CC_PD_ACETATE_BORON_TRANSFER', 'metal oxygen boron carbon', 'doi:10.1039/C1CY00241D; finite proposed acetate/boronic acid exchange in PdII reductive activation; microscopic sequence not established'),
 ('CC_SUZUKI_BORONATE_TRANSFER', 'metal leaving boron carbon', 'Kurti/Czako PDF502: separately declared tetrahedral boronate path; not the universal aqueous weak-base mechanism'),
 ('CC_HIYAMA_F_TRANSFER', 'metal leaving silicon carbon fluorine', 'Hiyama fluoride-activated pentacoordinate silicate, finite molecular exchange'),
 ('CC_REDUCTIVE_ELIMINATION', 'metal first second', 'Kurti/Czako molecular cis C-C/C-N/C-O/C-Sn reductive elimination'),
 ('CC_BASE_LIGAND_EXCHANGE', 'metal leaving donor hydrogen base', 'Base-assisted X-H substitution: retains H, leaving ion and proton acceptor'),
 ('CC_ANION_LIGAND_EXCHANGE', 'metal leaving donor', 'Explicit anionic ligand exchange, no oxidation-state overwrite'),
 ('CC_HECK_INSERTION', 'metal carbon first second', 'Kurti/Czako Heck: unresolved association within migratory insertion'),
 ('CC_BETA_H_ELIMINATION', 'metal alpha beta hydrogen', 'ILPI beta-hydride elimination; syn accessibility unverified in 2D'),
 ('CC_HX_BASE_RELEASE', 'metal hydrogen leaving base', 'Heck PdII hydride/halide plus base gives Pd0 and retained BH+/X-'),
 ('CC_ALLYL_ACETATE_ACTIVATION', 'metal first middle last oxygen', 'Kurti/Czako PDF512: localized sigma+pi allyl-PdII convention, simultaneous ionization/activation'),
 ('CC_SOFT_ALLYL_CAPTURE', 'metal first middle last nucleophile', 'Kurti/Czako PDF512: soft malonate carbon attack, same-terminus double-inversion branch'),
]}
SPECS['CC_PI_DISSOCIATION']=dict(roles=['metal','first','second'],event_type='COORDINATION',rule_revision=VERSION,
    source='doi:10.1021/ja403259c; doi:10.1021/ja9637098; finite dba ligand substitution hypothesis',library='organometallic')

def compile(v, e):
    spec=SPECS[e.general_event]; r=e.general_roles; ats=atoms_by_id(v)
    if set(r)!=set(spec['roles']) or len(set(r.values()))!=len(r) or not set(r.values())<=ats.keys():
        raise FormalError('Exact distinct coupling role references required')
    if e.event_type!='POLAR' or not e.evidence or not e.mechanism_assumption:
        raise FormalError('Coupling profile requires evidence and a declared molecular channel')
    ids=[r[k] for k in spec['roles']]; m=ids[0]
    counts={x['metal_atom_id']:x for x in checks(v)[1]}
    if m not in counts: raise FormalError('Coupling metal declaration missing')
    el=ats[m].element; os=counts[m]['oxidation_state']; name=e.general_event
    def bonds(rows):
        for a,b,n in rows:
            if order(v,a,b)!=n: raise FormalError(f'Coupling requires bond {a}-{b} order {n}, actual {order(v,a,b)}')
    def elements(items,allowed):
        if any(ats[i].element not in choices.split('/') for i,choices in zip(items,allowed.split())):
            raise Unsupported('Element roles outside registered coupling channel')
    def acceptor():
        if (el,os) not in {('Pd',2),('Ni',2),('Cu',3)}:
            raise Unsupported('Molecular PdII/NiII/CuIII coupling acceptor required; SET is a separate channel')
    if name=='CC_OXIDATIVE_ADDITION':
        m,c,x=ids;elements([c,x],'C Cl/Br/I/O')
        if ats[x].element=='O':
            triflate=False
            for s in ats:
                if ats[s].element!='S' or order(v,x,s)!=1:continue
                ox=sum(ats[j].element=='O' and order(v,s,j)==2 for j in ats)
                cf=[j for j in ats if ats[j].element=='C' and order(v,s,j)==1 and sum(ats[f].element=='F' and order(v,j,f)==1 for f in ats)==3]
                triflate=bool(ox==2 and cf)
            if not triflate or el not in {'Pd','Ni'}:raise Unsupported('C-O activation is restricted to explicit Pd/Ni aryl/vinyl triflate input')
        if (el,os) not in {('Pd',0),('Ni',0),('Cu',1)} or counts[m]['d_electrons']!=10:
            raise Unsupported('Declared d10 two-electron oxidative addition only')
        if ats[m].radical_electrons not in (None,0):raise Unsupported('Open-shell local donor needs a different channel')
        bonds([(c,x,1),(m,c,0),(m,x,0)])
        # Alkyl Ni/aryl fluoride are not silently included in this paired branch.
        if el=='Ni' and bond_valences(v)[c]==4 and not any(order(v,c,j)>1 for j in ats):
            raise Unsupported('Saturated nickel activation needs a separately evidenced radical/paired branch')
        return [Pair(a=m,b=c),Shift(a=c,b=x,c=m)]
    if name=='CC_ALLYL_ACETATE_ACTIVATION':
        m,a,b,c,o=ids;elements([a,b,c,o],'C C C O')
        if (el,os)!=('Pd',0) or counts[m]['d_electrons']!=10 or ats[m].radical_electrons not in (None,0):raise Unsupported('Paired Pd0 allyl activation only')
        bonds([(a,b,1),(b,c,2),(a,o,1),(m,a,0)])
        if not any(ats[k].element=='C' and order(v,o,k)==1 and any(ats[z].element=='O' and order(v,k,z)==2 for z in ats) for k in ats):
            raise Unsupported('Allylic carboxylate leaving group required')
        return [Pair(a=m,b=a),Split(a=a,b=o,target=o)]
    if name=='CC_SOFT_ALLYL_CAPTURE':
        m,a,b,c,n=ids;elements([a,b,c,n],'C C C C');acceptor()
        if el!='Pd' or (ats[n].formal_charge,ats[n].nonbonding_electrons,ats[n].radical_electrons)!=(-1,2,0):raise Unsupported('Paired malonate carbanion/PdII branch only')
        if sum(ats[k].element=='C' and order(v,n,k)==1 and any(ats[z].element=='O' and order(v,k,z)==2 for z in ats) for k in ats)!=2:
            raise Unsupported('Two carbonyl-stabilized soft-carbon substituents required')
        bonds([(m,a,1),(a,b,1),(b,c,2),(a,n,0)])
        if not any(x.metal_atom_id==m and x.mode=='eta' and set(x.donor_atom_ids)=={b,c} for x in v.metal_bindings):
            raise FormalError('Localized sigma+pi allyl binding required')
        return [Pair(a=n,b=a),Split(a=a,b=m,target=m)]
    if name=='CC_TRANSMETALATION':
        acceptor();m,x,d,c=ids;elements([x,c],'Cl/Br/I/O C/Sn')
        de=ats[d].element
        if de not in {'Zn','Mg','Sn','Cu'}:raise Unsupported('Distinct organozinc/magnesium/tin/copper donor required')
        if de in {'Zn','Mg','Cu'} and (d not in counts or counts[d]['oxidation_state']!={'Zn':2,'Mg':2,'Cu':1}[de]):
            raise Unsupported('Donor oxidation state outside coupling channel')
        if de=='Sn' and (bond_valences(v)[d],ats[d].formal_charge,ats[d].nonbonding_electrons)!=(4,0,0):
            raise Unsupported('Neutral tetracoordinate tin donor required')
        if ats[c].element=='Sn' and de!='Sn':raise Unsupported('Distannane stannylation is a separate donor subcase')
        bonds([(m,x,1),(d,c,1),(m,c,0),(d,x,0)])
        return [Shift(a=m,b=x,c=d),Shift(a=d,b=c,c=m)]
    if name in {'CC_SUZUKI_OH_TRANSFER','CC_PD_ACETATE_BORON_TRANSFER'}:
        acceptor();m,o,b,c=ids;elements([o,b,c],'O B C')
        if (bond_valences(v)[b],ats[b].formal_charge,ats[b].nonbonding_electrons)!=(3,0,0):
            raise Unsupported('Neutral tricoordinate boronic acid/ester only; not every B reagent')
        if sum(ats[j].element=='O' and order(v,b,j)==1 for j in ats)!=2:raise Unsupported('Two B-O substituents required')
        if name=='CC_SUZUKI_OH_TRANSFER':
            if not any(ats[j].element=='H' and order(v,o,j)==1 for j in ats):raise FormalError('Pd hydroxo hydrogen must be explicit')
        else:
            if (el,os)!=('Pd',2):raise Unsupported('Acetate activation branch requires PdII')
            acetates=[j for j in ats if ats[j].element=='C' and order(v,o,j)==1 and
                any(ats[k].element=='O' and order(v,j,k)==2 for k in ats) and
                any(ats[k].element=='C' and order(v,j,k)==1 and bond_valences(v)[k]==4 for k in ats)]
            if not acetates:raise Unsupported('Explicit bound acetate required; other Pd salts need separate activation')
        bonds([(m,o,1),(b,c,1),(m,c,0),(b,o,0)])
        return [Shift(a=m,b=o,c=b),Shift(a=b,b=c,c=m)]
    if name=='CC_SUZUKI_BORONATE_TRANSFER':
        acceptor();m,x,b,c=ids;elements([x,b,c],'Cl/Br/I B C')
        if (bond_valences(v)[b],ats[b].formal_charge,ats[b].nonbonding_electrons)!=(4,-1,0):raise Unsupported('Explicit tetrahedral boronate anion required')
        if sum(ats[j].element=='O' and order(v,b,j)==1 for j in ats)!=3:raise Unsupported('Three B-O bonds required for this boronate channel')
        bonds([(m,x,1),(b,c,1),(m,c,0)])
        return [Split(a=m,b=x,target=x),Shift(a=b,b=c,c=m)]
    if name=='CC_HIYAMA_F_TRANSFER':
        acceptor();m,x,s,c,f=ids;elements([x,s,c,f],'Cl/Br/I Si C F')
        if (bond_valences(v)[s],ats[s].formal_charge)!=(5,-1):raise Unsupported('Explicit fluoride-activated pentacoordinate silicate required')
        bonds([(m,x,1),(s,c,1),(s,f,1),(m,c,0)])
        return [Split(a=m,b=x,target=x),Shift(a=s,b=c,c=m)]
    if name=='CC_REDUCTIVE_ELIMINATION':
        acceptor();m,a,b=ids;elements([a,b],'C C/N/O/Sn')
        bonds([(m,a,1),(m,b,1),(a,b,0)])
        center=next(x for x in v.metal_centers if x.atom_id==m)
        bindings={x.donor_atom_ids[0]:x.binding_id for x in v.metal_bindings if x.metal_atom_id==m and len(x.donor_atom_ids)==1}
        if any({p.first_binding_id,p.second_binding_id}=={bindings.get(a),bindings.get(b)} and p.relation=='trans' for p in center.relative_positions):
            raise FormalError('Declared trans ligands cannot undergo cis reductive elimination without isomerization')
        return [Split(a=a,b=m,target=m),Shift(a=m,b=b,c=a)]
    if name in {'CC_BASE_LIGAND_EXCHANGE','CC_ANION_LIGAND_EXCHANGE'}:
        if (el,os) not in {('Pd',2),('Ni',2),('Cu',1),('Cu',3)}:raise Unsupported('Registered ligand exchange metal/OS required')
        m,x,d=ids[:3];elements([x,d],'Cl/Br/I/O C/N/O')
        bonds([(m,x,1),(m,d,0)])
        if name=='CC_ANION_LIGAND_EXCHANGE':
            if ats[d].formal_charge!=-1 or ats[d].radical_electrons!=0:raise FormalError('Paired anionic donor required')
            return [Pair(a=d,b=m),Split(a=m,b=x,target=x)]
        h,b=ids[3:];elements([h,b],'H N/O');bonds([(d,h,1),(b,h,0)])
        if ats[d].element=='C' and not any(order(v,d,j)==3 and ats[j].element=='C' for j in ats):
            raise Unsupported('Base-assisted C-H metalation restricted to terminal acetylene')
        if ats[b].radical_electrons!=0 or ats[b].nonbonding_electrons<2:raise FormalError('Explicit paired proton acceptor required')
        return [Pair(a=b,b=h),Shift(a=h,b=d,c=m),Split(a=m,b=x,target=x)]
    if name=='CC_HECK_INSERTION':
        acceptor();m,c,a,b=ids;elements([c,a,b],'C C C')
        bonds([(m,c,1),(a,b,2),(c,a,0),(m,b,0)])
        if e.conditions.get('association_model')!='unresolved_within_event':raise Unsupported('Unresolved alkene association must be declared')
        return [Shift(a=m,b=c,c=a),Shift(a=a,b=b,c=m)]
    if name=='CC_BETA_H_ELIMINATION':
        acceptor();m,a,b,h=ids;elements([a,b,h],'C C H')
        bonds([(m,a,1),(a,b,1),(b,h,1),(m,h,0)])
        if e.conditions.get('syn_accessibility')=='excluded':raise FormalError('Explicitly inaccessible beta-H geometry')
        return [Shift(a=m,b=a,c=b),Shift(a=b,b=h,c=m)]
    if name=='CC_HX_BASE_RELEASE':
        acceptor();m,h,x,b=ids;elements([h,x,b],'H Cl/Br/I/O N/O')
        bonds([(m,h,1),(m,x,1),(b,h,0)])
        if ats[b].nonbonding_electrons<2 or ats[b].radical_electrons!=0:raise FormalError('Paired base required')
        return [Pair(a=b,b=h),Split(a=h,b=m,target=m),Split(a=m,b=x,target=x)]
    raise Unsupported('Unimplemented coupling channel')

def execute(v,e):
    from .execution import transact
    from .ledger import validate_view
    from .metal_events import repartition
    from ..v06.stereo import apply_event
    validate_view(v)
    if e.general_event=='CC_PI_DISSOCIATION':
        r=e.general_roles
        if set(r)!=set(SPECS[e.general_event]['roles']) or len(set(r.values()))!=3 or e.event_type!='COORDINATION' or e.actions or not e.evidence or not e.mechanism_assumption:
            raise FormalError('Typed evidenced pi dissociation without covalent actions required')
        m,a,b=[r[k] for k in SPECS[e.general_event]['roles']]
        hits=[x for x in v.metal_bindings if x.metal_atom_id==m and x.mode=='eta' and x.hapticity==2 and set(x.donor_atom_ids)=={a,b}]
        if len(hits)!=1 or atoms_by_id(v)[m].element!='Pd' or order(v,a,b)!=2:raise Unsupported('Existing Pd/alkene pi pool required')
        t=v.model_copy(deep=True);t.state_id=e.product_state_id;t.view_id=e.product_view_id
        t.metal_bindings=[x for x in t.metal_bindings if x.binding_id!=hits[0].binding_id]
        for metal in t.metal_centers:
            if metal.atom_id==m:metal.valence_electrons=None;metal.relative_positions=[];metal.vacant_site_ids=None
        e.actions=[]
        return t,dict(event_id=e.event_id,event_type=e.event_type,general_event=e.general_event,rule_version=VERSION,
            execution_mode='electron_actions',formal_status='passed',atomic_commit=True,compiled_actions=[],ledger=validate_view(t),
            chemical_feasibility='unknown',coordination_change={'released_pi_pool':[a,b],'metal':m},electron_ownership_changed=False)
    actions=compile(v,e)
    key=lambda a:json.dumps(a.model_dump(mode='json'),sort_keys=True)
    if e.actions is not None and Counter(map(key,e.actions))!=Counter(map(key,actions)):
        raise FormalError('Actions disagree with registered coupling profile')
    donors={e.general_roles['metal']} if e.general_event in {'CC_OXIDATIVE_ADDITION','CC_ALLYL_ACETATE_ACTIVATION'} else set()
    replacements=[];r=e.general_roles
    if e.general_event=='CC_OXIDATIVE_ADDITION':replacements=[(r['carbon'],r['leaving'],r['metal'])]
    if e.general_event=='CC_TRANSMETALATION':replacements=[(r['carbon'],r['donor'],r['metal'])]
    if e.general_event in {'CC_SUZUKI_OH_TRANSFER','CC_PD_ACETATE_BORON_TRANSFER'}:replacements=[(r['carbon'],r['boron'],r['metal'])]
    if e.general_event=='CC_SUZUKI_BORONATE_TRANSFER':replacements=[(r['carbon'],r['boron'],r['metal'])]
    if e.general_event=='CC_HIYAMA_F_TRANSFER':replacements=[(r['carbon'],r['silicon'],r['metal'])]
    if e.general_event=='CC_REDUCTIVE_ELIMINATION':replacements=[(r['first'],r['metal'],r['second']),(r['second'],r['metal'],r['first'])]
    stereo=[]
    released={b.binding_id for b in v.metal_bindings if e.general_event=='CC_SOFT_ALLYL_CAPTURE' and b.metal_atom_id==r['metal'] and b.mode=='eta' and set(b.donor_atom_ids)=={r['middle'],r['last']}}
    def after(x):
        x.metal_bindings=[b for b in x.metal_bindings if b.binding_id not in released]
        repartition(x)
        if e.general_event=='CC_ALLYL_ACETATE_ACTIVATION':
            from .ast import MetalBinding
            ligand=next(l for l in x.ligands if r['first'] in l.atom_ids)
            ligand.ionic_class='LX';ligand.donor_atom_ids=sorted(set(ligand.donor_atom_ids)|{r['middle'],r['last']})
            ligand.ligand_type='allyl';ligand.assignment_basis='localized sigma+pi representation of allyl-Pd; whole event atomic'
            x.metal_bindings.append(MetalBinding(binding_id=e.event_id+'.pi',ligand_id=ligand.ligand_id,metal_atom_id=r['metal'],
                donor_atom_ids=[r['middle'],r['last']],mode='eta',hapticity=2,ionic_allocation=0,donated_electrons=2,pi_pool_slot=0))
        if e.general_event in {'CC_ALLYL_ACETATE_ACTIVATION','CC_SOFT_ALLYL_CAPTURE'}:
            c=r['first'];old=atoms_by_id(v)[c].stereo
            departing=r['oxygen'] if e.general_event=='CC_ALLYL_ACETATE_ACTIVATION' else r['metal']
            arriving=r['metal'] if e.general_event=='CC_ALLYL_ACETATE_ACTIVATION' else r['nucleophile']
            if old.get('status')=='specified' and departing in old.get('neighbor_order',[]) and old.get('rdkit_tag') in {'CHI_TETRAHEDRAL_CW','CHI_TETRAHEDRAL_CCW'}:
                atoms_by_id(x)[c].stereo={**old,'neighbor_order':[arriving if i==departing else i for i in old['neighbor_order']],
                   'rdkit_tag':'CHI_TETRAHEDRAL_CCW' if old['rdkit_tag']=='CHI_TETRAHEDRAL_CW' else 'CHI_TETRAHEDRAL_CW',
                   'basis':'declared same-terminus soft allylation double-inversion channel; not universal stereoselectivity'}
        # These registered paired channels retain the relative geometry of an
        # UNBROKEN vinyl pi bond. This is not a product-guided E/Z assignment,
        # nor a claim that every metal catalyst forbids isomerization.
        for old in v.bonds:
            if old.order!=2 or old.stereo.get('status')!='specified':continue
            if not any(c in (old.a,old.b) for c,_,_ in replacements):continue
            new=next((b for b in x.bonds if {b.a,b.b}=={old.a,old.b} and b.order==2),None)
            if new is None:continue
            refs=list(old.stereo.get('reference_atoms',[]))
            if len(refs)!=2:continue
            for c,depart,arrive in replacements:
                if c in (old.a,old.b):
                    k=0 if c==old.a else 1
                    if refs[k]==depart:refs[k]=arrive
            if any(order(x,end,ref)!=1 for end,ref in zip((old.a,old.b),refs)):continue
            new.stereo={**old.stereo,'reference_atoms':refs if new.a==old.a else refs[::-1],
                        'basis':'registered stereoretentive vinyl coupling channel; CIP independently recomputed'}
            stereo.append(dict(bond=[old.a,old.b],reference_atoms=refs,effect='relative_vinyl_retention'))
    t=transact(v,actions,e.product_state_id,e.product_view_id,after=after,_declared_pair_donors=donors,_released_binding_ids=released)
    e.actions=actions
    return t,dict(event_id=e.event_id,event_type=e.event_type,general_event=e.general_event,
        rule_version=VERSION,execution_mode='electron_actions',atomic_commit=True,formal_status='passed',
        compiled_actions=[a.model_dump(mode='json') for a in actions],stereochemistry=apply_event(v,t,e),
        ledger=validate_view(t),chemical_feasibility='unknown',geometry_status='unverified',
        microscopic_association='not_resolved',vinyl_geometry=stereo,declared_pair_donor_ids=sorted(donors),
        pair_occupation_basis='declared molecular two-electron channel; total spin unchanged')
