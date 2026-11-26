"""Role-bound r10 metal-layer transactions over existing PAIR/SPLIT/SHIFT.

Not a named-reaction dispatcher. Grouped stoichiometric redox channels are
explicitly excluded from complete elementary-path certification.
"""
import json
from collections import Counter
from ..ast import Pair, Split, Shift
from ..chemistry import FormalError, Unsupported, atoms_by_id, order, bond_valences, identity, electron_count
from .coordination import checks
from .ledger import validate_view

VERSION='umr-rules/0.9-r10'
SPECS={name:dict(roles=roles.split(),event_type='POLAR',rule_revision='0.9-r10',
    library='organometallic',source=source,geometry_status='unverified') for name,roles,source in [
 ('ALKYLIDYNE_ALKYNE_CAPTURE','metal carbyne a b','Schrock alkylidyne/MCBD channel; doi:10.1021/acscatal.1c01587'),
 ('ALKYLIDYNE_CYCLOREVERSION','metal first second third','Same MCBD channel; registered Lewis redistribution is resonance, not a third real event'),
 ('RU_ALKYNE_CAPTURE','metal carbene a b','Kurti/Czako PDF206 path II; doi:10.1021/ja907155r; channel dependent'),
 ('RU_ALKENE_CAPTURE','metal carbene a b','Declared Ru alkylidene/alkene channel'),
 ('RU_METALLACYCLE_OPENING','metal first second third','Explicit Ru metallacyclobutane/butene topology; retained vinyl carbene'),
 ('CO_CO_EPOXIDE_HYDROLYSIS','donor_metal acceptor_metal hydroxy_o hydrogen carbon ring_o other_c','CoIII salen cooperative HKR; doi:10.1021/ja038590z'),
 ('CO_HYDROXO_REGENERATION','metal acceptor_o water_o hydrogen','CoIII water/proton transfer regeneration; explicit water replaces transferred OH'),
 ('METAL_ALKOXIDE_WATER_EXCHANGE','metal product_o water_o hydrogen','Finite protonolysis/ligand exchange; explicit water atom mapping'),
 ('TI_ASSISTED_ALLYLTIN_ADDITION','metal tin first middle terminal carbon oxygen leaving_o','Keck declared monomer channel, not aggregation/ee prediction; doi:10.1021/ja00071a074'),
 ('TI_SN_ALKOXIDE_EXCHANGE','metal tin product_o tin_o','Source-proposed titanium product alkoxide/tributyltin alkoxide exchange; no OS change'),
 ('SN_ALKOXIDE_HYDROLYSIS','tin product_o water_o hydrogen','Finite tetraorganotin alkoxide hydrolysis; explicit water'),
 ('MG_TI_ALKOXIDE_ALKYL_EXCHANGE','metal magnesium carbon oxygen','Kulinkovich precursor ligand exchange; original donor/acceptor retained'),
 ('TI_DIALKYL_CYCLOPROPANATION','metal alpha beta departing hydrogen','Grouped dialkyltitanium accounting; alternative transfer-epititanation exists: doi:10.1002/ejoc.200300588'),
 ('TI_CARBONYL_INSERTION','metal carbon carbonyl oxygen','Finite TiIV sigma carbon to ester/ketone, no general Grignard premise'),
 ('TI_TETRAHEDRAL_COLLAPSE','metal carbonyl oxygen leaving_o','Finite oxatitanacycle collapse with retained alkoxide'),
 ('RU_OXO_ALCOHOL_ASSOCIATION','metal oxo alcohol_o hydrogen','Molecular TPAP model; competing RuO2 chemistry retained as uncertainty'),
 ('RU_ALKOXIDE_OXIDATION','metal carbon alcohol_o hydrogen hydroxy_o','RuVII -> RuV formal paired elimination; doi:10.1039/C6CC05440D'),
 ('RU_N_OXIDE_REOXIDATION','metal nitrogen oxygen','RuV rescue by N-oxide; oxygen transfer, not arbitrary OS overwrite'),
 ('SNII_IMIDOYL_CHLORIDE_REDUCTION_GROUP','tin carbon nitrogen leaving_cl hydrogen donor_cl','Stephen protonated imidoyl chloride/SnCl2/HCl grouped reduction, not a resolved microscopic SET sequence'),
]}
GROUPED={'TI_DIALKYL_CYCLOPROPANATION','SNII_IMIDOYL_CHLORIDE_REDUCTION_GROUP'}


def compile(v,e):
    spec=SPECS[e.general_event];r=e.general_roles;ats=atoms_by_id(v)
    if set(r)!=set(spec['roles']) or len(set(r.values()))!=len(r) or not set(r.values())<=ats.keys():raise FormalError('Exact distinct catalytic role references required')
    if e.event_type!='POLAR' or not e.evidence or not e.mechanism_assumption:raise FormalError('Typed event, evidence and mechanism assumption required')
    if e.general_event in GROUPED and e.conditions.get('mechanistic_resolution')!='grouped_ledger':raise Unsupported('Unresolved microscopic channel requires explicit grouped_ledger declaration')
    ids=[r[k] for k in spec['roles']];counts={c['metal_atom_id']:c for c in checks(v)[1]}
    def metal(i,element,os):
        if i not in counts or ats[i].element not in set(element.split('/')) or counts[i]['oxidation_state']!=os:raise Unsupported('Metal/oxidation input outside catalytic profile')
    def elements(items,types):
        if [ats[i].element for i in items]!=types.split():raise Unsupported('Element roles outside catalytic profile')
    def bonds(items):
        for a,b,n in items:
            if order(v,a,b)!=n:raise FormalError(f'Catalytic profile requires bond {a}-{b} order {n}')
    def bound(m,d):
        if not any(b.metal_atom_id==m and d in b.donor_atom_ids and b.mode=='dative' for b in v.metal_bindings):raise FormalError('Pre-associated donor required')
    name=e.general_event
    if name.startswith('ALKYLIDYNE_'):
        m,c,a,b=ids;metal(m,'W/Mo',6);elements([c,a,b],'C C C')
        if sum(ats[j].element=='O' and order(v,m,j)==1 for j in ats)!=3:raise Unsupported('Three covalent alkoxide/aryloxide donors required')
        if name=='ALKYLIDYNE_ALKYNE_CAPTURE':
            bonds([(m,c,3),(a,b,3),(c,a,0),(m,b,0)])
            if not any(l.ligand_type=='Schrock_alkylidyne' and l.ionic_class=='X3' and c in l.donor_atom_ids for l in v.ligands):raise Unsupported('Declared X3 alkylidyne required')
        else:bonds([(m,c,1),(c,a,2),(a,b,1),(b,m,2)])
        return [Shift(a=m,b=c,c=a),Shift(a=a,b=b,c=m)]
    if name in {'RU_ALKYNE_CAPTURE','RU_ALKENE_CAPTURE','RU_METALLACYCLE_OPENING'}:
        m,c,a,b=ids;elements([c,a,b],'C C C')
        if name=='RU_METALLACYCLE_OPENING':
            metal(m,'Ru',4);bonds([(m,c,1),(c,a,1),(m,b,1)])
            if order(v,a,b) not in {1,2}:raise FormalError('Ru four-membered ring internal order must be one or two')
        else:
            metal(m,'Ru',2);bonds([(m,c,2),(a,b,3 if name=='RU_ALKYNE_CAPTURE' else 2),(c,a,0),(m,b,0)])
            if not any(l.ligand_type=='Ru_alkylidene' and l.ionic_class=='L' and c in l.donor_atom_ids for l in v.ligands):raise Unsupported('Explicit Ru L-type alkylidene convention required')
        return [Shift(a=m,b=c,c=a),Shift(a=a,b=b,c=m)]
    if name=='CO_CO_EPOXIDE_HYDROLYSIS':
        m,n,o,h,c,x,a=ids;metal(m,'Co',3);metal(n,'Co',3);elements([o,h,c,x,a],'O H C O C')
        for center in (m,n):
            if not any(l.ionic_class=='L2X2' and len(l.donor_atom_ids)==4 and any(b.metal_atom_id==center and b.ligand_id==l.ligand_id for b in v.metal_bindings) for l in v.ligands):raise Unsupported('Two complete tetradentate dianionic salen ligands required')
        bound(n,x);bonds([(m,o,1),(o,h,1),(c,x,1),(x,a,1),(a,c,1),(o,c,0)])
        if bond_valences(v)[o]!=2 or ats[o].nonbonding_electrons!=4:raise FormalError('Co hydroxo donor must retain its explicit H')
        return [Shift(a=m,b=o,c=c),Split(a=c,b=x,target=x)]
    if name in {'METAL_ALKOXIDE_WATER_EXCHANGE','SN_ALKOXIDE_HYDROLYSIS'}:
        m,o,w,h=ids;elements([o,w,h],'O O H');bonds([(m,o,1),(w,h,1),(w,m,0),(o,h,0)])
        if name=='SN_ALKOXIDE_HYDROLYSIS':
            if ats[m].element!='Sn' or (bond_valences(v)[m],ats[m].nonbonding_electrons,ats[m].formal_charge)!=(4,0,0):raise Unsupported('Neutral SnIV alkoxide only')
        elif m not in counts or (ats[m].element,counts[m]['oxidation_state']) not in {('Ti',4),('Co',3)}:raise Unsupported('TiIV/CoIII product alkoxide only')
        if bond_valences(v)[w]!=2 or sum(order(v,w,j)==1 and ats[j].element=='H' for j in ats)!=2:raise FormalError('Water with both mapped hydrogens required')
        return [Pair(a=w,b=m),Split(a=m,b=o,target=o),Pair(a=o,b=h),Split(a=h,b=w,target=w)]
    if name=='CO_HYDROXO_REGENERATION':
        m,o,w,h=ids;metal(m,'Co',3);elements([o,w,h],'O O H');bonds([(w,h,1),(m,w,0),(o,h,0)])
        if ats[o].formal_charge!=-1 or ats[m].formal_charge!=1 or bond_valences(v)[w]!=2:raise FormalError('Cationic CoIII and product alkoxide/water required')
        if sum(order(v,w,j)==1 and ats[j].element=='H' for j in ats)!=2:raise FormalError('Water with both mapped hydrogens required')
        return [Pair(a=w,b=m),Pair(a=o,b=h),Split(a=h,b=w,target=w)]
    if name=='TI_ASSISTED_ALLYLTIN_ADDITION':
        m,s,a,b,c,k,o,x=ids;metal(m,'Ti',4);elements([s,a,b,c,k,o,x],'Sn C C C C O O');bound(m,o)
        bonds([(s,a,1),(a,b,1),(b,c,2),(k,o,2),(c,k,0),(m,o,0),(m,x,1),(s,x,0)])
        if (bond_valences(v)[s],ats[s].formal_charge,ats[s].nonbonding_electrons)!=(4,0,0):raise FormalError('Neutral tetracoordinate allylstannane required')
        return [Shift(a=s,b=a,c=b),Shift(a=b,b=c,c=k),Split(a=k,b=o,target=o),Pair(a=o,b=m),Shift(a=m,b=x,c=s)]
    if name=='TI_SN_ALKOXIDE_EXCHANGE':
        m,s,o,x=ids;metal(m,'Ti',4);elements([s,o,x],'Sn O O');bonds([(m,o,1),(s,x,1),(s,o,0),(m,x,0)])
        if (bond_valences(v)[s],ats[s].formal_charge,ats[s].nonbonding_electrons)!=(4,0,0):raise FormalError('SnIV alkoxide exchange only')
        return [Shift(a=m,b=o,c=s),Shift(a=s,b=x,c=m)]
    if name=='MG_TI_ALKOXIDE_ALKYL_EXCHANGE':
        m,g,c,o=ids;metal(m,'Ti',4);metal(g,'Mg',2);elements([c,o],'C O');bonds([(m,o,1),(g,c,1),(m,c,0),(g,o,0)])
        if not any(ats[j].element in {'Cl','Br','I'} and order(v,g,j)==1 for j in ats):raise Unsupported('Explicit organomagnesium halide required')
        return [Shift(a=m,b=o,c=g),Shift(a=g,b=c,c=m)]
    if name=='TI_DIALKYL_CYCLOPROPANATION':
        m,a,b,c,h=ids;metal(m,'Ti',4);elements([a,b,c,h],'C C C H');bonds([(m,a,1),(a,b,1),(m,c,1),(b,h,1),(m,b,0),(c,h,0)])
        return [Shift(a=h,b=b,c=m),Shift(a=m,b=c,c=h)]
    if name=='TI_CARBONYL_INSERTION':
        m,c,k,o=ids;metal(m,'Ti',4);elements([c,k,o],'C C O');bonds([(m,c,1),(k,o,2),(c,k,0),(m,o,0)])
        return [Shift(a=m,b=c,c=k),Split(a=k,b=o,target=o),Pair(a=o,b=m)]
    if name=='TI_TETRAHEDRAL_COLLAPSE':
        m,k,o,x=ids;metal(m,'Ti',4);elements([k,o,x],'C O O');bonds([(m,o,1),(k,o,1),(k,x,1),(m,x,0)])
        if bond_valences(v)[k]!=4:raise FormalError('Tetrahedral acyl carbon required')
        return [Pair(a=o,b=k),Split(a=k,b=x,target=x),Split(a=m,b=o,target=o),Pair(a=x,b=m)]
    if name=='RU_OXO_ALCOHOL_ASSOCIATION':
        m,x,o,h=ids;metal(m,'Ru',7);elements([x,o,h],'O O H');bonds([(m,x,2),(o,h,1),(m,o,0),(x,h,0)])
        if sum(order(v,m,j)==2 and ats[j].element=='O' for j in ats)!=4:raise Unsupported('Molecular tetraoxo RuVII anion required')
        return [Pair(a=o,b=m),Split(a=m,b=x,target=x),Pair(a=x,b=h),Split(a=h,b=o,target=o)]
    if name=='RU_ALKOXIDE_OXIDATION':
        m,c,o,h,x=ids;metal(m,'Ru',7);elements([c,o,h,x],'C O H O');bonds([(m,o,1),(o,c,1),(c,h,1),(m,x,1),(x,h,0)])
        if sum(order(v,x,j)==1 and ats[j].element=='H' for j in ats)!=1:raise FormalError('Ru hydroxo leaving group needs its retained H')
        return [Shift(a=h,b=c,c=o),Split(a=o,b=m,target=m),Pair(a=x,b=h),Split(a=m,b=x,target=x)]
    if name=='RU_N_OXIDE_REOXIDATION':
        m,n,o=ids;metal(m,'Ru',5);elements([n,o],'N O');bonds([(n,o,1),(m,o,0)])
        if (ats[n].formal_charge,ats[o].formal_charge)!=(1,-1) or (ats[m].nonbonding_electrons,ats[m].radical_electrons)!=(3,1):raise Unsupported('Explicit RuV paired-donor occupation and amine N-oxide required')
        return [Pair(a=o,b=m),Pair(a=m,b=o),Split(a=o,b=n,target=n)]
    if name=='SNII_IMIDOYL_CHLORIDE_REDUCTION_GROUP':
        s,c,n,x,h,y=ids;elements(ids,'Sn C N Cl H Cl')
        if (bond_valences(v)[s],ats[s].nonbonding_electrons,ats[s].radical_electrons,ats[s].formal_charge)!=(2,2,0,0):raise Unsupported('Paired neutral SnII dihalide only')
        bonds([(c,n,2),(c,x,1),(h,y,1),(s,x,0),(s,y,0),(c,h,0)])
        if ats[n].formal_charge!=1 or sum(order(v,n,j)==1 and ats[j].element=='H' for j in ats)!=2:raise FormalError('Protonated imidoyl chloride with two mapped N-H atoms required')
        return [Pair(a=s,b=x),Shift(a=x,b=c,c=h),Shift(a=h,b=y,c=s)]
    raise Unsupported('Unimplemented catalytic profile')


def postprocess(v,e):
    from .metal_events import repartition
    repartition(v)
    if e.general_event=='RU_METALLACYCLE_OPENING':
        m=e.general_roles['metal'];c=e.general_roles['third']
        for l in v.ligands:
            if c in l.donor_atom_ids and order(v,m,c)==2:
                l.ligand_type='Ru_alkylidene';l.ionic_class='L';l.ionic_charge=0
                for b in v.metal_bindings:
                    if b.ligand_id==l.ligand_id:b.donated_electrons=2


def execute(v,e):
    from .execution import transact
    from ..v06.stereo import apply_event
    validate_view(v);actions=compile(v,e)
    key=lambda a:json.dumps(a.model_dump(mode='json'),sort_keys=True)
    if e.actions is not None and Counter(map(key,e.actions))!=Counter(map(key,actions)):raise FormalError('Actions disagree with catalytic profile')
    released={b.binding_id for b in v.metal_bindings if e.general_event=='TI_ASSISTED_ALLYLTIN_ADDITION' and b.mode=='dative' and b.metal_atom_id==e.general_roles['metal'] and b.donor_atom_ids==[e.general_roles['oxygen']]}
    def after(x):
        x.metal_bindings=[b for b in x.metal_bindings if b.binding_id not in released]
        postprocess(x,e)
    t=transact(v,actions,e.product_state_id,e.product_view_id,after=after,_released_binding_ids=released)
    if identity(t)!=identity(v) or electron_count(t)!=electron_count(v) or sum(a.formal_charge for a in t.atoms)!=sum(a.formal_charge for a in v.atoms):raise FormalError('Catalytic material/electron conservation failure')
    e.actions=actions
    report=dict(event_id=e.event_id,event_type=e.event_type,general_event=e.general_event,
        rule_version=VERSION,execution_mode='electron_actions',atomic_commit=True,formal_status='passed',
        compiled_actions=[a.model_dump(mode='json') for a in actions],stereochemistry=apply_event(v,t,e),
        chemical_feasibility='unknown',library='organometallic',metal_elements=sorted({a.element for a in v.atoms if a.element in {'W','Mo','Ru','Co','Ti','Sn','Mg'}}))
    if e.general_event in GROUPED:
        report['mechanistic_resolution']='grouped_redox_proton_ledger_not_asserted_elementary'
    return t,report
