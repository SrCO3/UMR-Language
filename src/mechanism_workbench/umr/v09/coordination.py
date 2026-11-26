"""Donor-owned ligand accounting with explicit multimetal allocation."""
from collections import Counter
from ..chemistry import atoms_by_id,bond_orders,FormalError,Unsupported
from .elements import LEGACY, current_elements
from .local_rules import ValenceRule, active_rules, require_equal
CLASSES={'L':(0,2),'X':(-1,2),'X2':(-2,4),'X3':(-3,6),'L2':(0,4),'LX':(-1,4),'L2X':(-1,6),'LX2':(-2,6),'L2X2':(-2,8)}

def neutral_ru_carbene(view,ligand):
    """Finite L2 chelated ether or alkylidene+pi scope, not a label bypass."""
    if ligand.ligand_type not in {'Ru_alkylidene','Fischer_carbene'}:return False
    if ligand.ionic_class=='L':return True  # retain historical convention
    if ligand.ionic_class!='L2' or ligand.ligand_type!='Ru_alkylidene':return False
    ats=atoms_by_id(view);edges=bond_orders(view)
    bs=[b for b in view.metal_bindings if b.ligand_id==ligand.ligand_id]
    sig=[b for b in bs if b.mode=='sigma_covalent']
    other=[b for b in bs if b.mode!='sigma_covalent']
    if len(sig)!=1 or len(other)!=1 or len(sig[0].donor_atom_ids)!=1:return False
    m=sig[0].metal_atom_id;c=sig[0].donor_atom_ids[0];b=other[0]
    if ats[m].element!='Ru' or b.metal_atom_id!=m or ats[c].element!='C' or edges.get(tuple(sorted((m,c))))!=2:return False
    if ligand.ionic_charge!=0 or sum(ats[i].formal_charge for i in ligand.atom_ids)!=0:return False
    if b.mode=='dative':
        if len(b.donor_atom_ids)!=1:return False
        o=ats[b.donor_atom_ids[0]]
        return o.element=='O' and (o.formal_charge,o.nonbonding_electrons,o.radical_electrons)==(0,4,0)
    return (b.mode=='eta' and b.hapticity==2 and len(b.donor_atom_ids)==2 and c not in b.donor_atom_ids
            and all(ats[i].element=='C' for i in b.donor_atom_ids) and edges.get(tuple(sorted(b.donor_atom_ids)))==2)
def checks(view,*,_symbolic_ids=frozenset()):
    rules = current_elements()
    METALS, OXIDATION = rules.metals, rules.oxidation
    # Legacy transactions construct v0.3 ligand nodes. Adapt a private copy at
    # the boundary, never interpreting absent v0.9 fields as new permissions.
    from ..v07.ledger import symbols,project_symbols
    symbolic=symbols(view)|set(_symbolic_ids)
    if symbols(view):view=project_symbols(view)
    try:view=rules.view_type.model_validate(view.model_dump(mode='json'))
    except ValueError as exc:raise FormalError(str(exc)) from exc
    atoms=atoms_by_id(view);edges=bond_orders(view);centers={m.atom_id:m for m in view.metal_centers}
    mids={i for i,a in atoms.items() if a.element in METALS}
    if set(centers)!=mids:raise FormalError('Each metal requires a unique state-local MetalCenter')
    ligands={l.ligand_id:l for l in view.ligands};reserved=Counter();covered=set();alloc=Counter();donations=Counter();seen=set();sites=set()
    if len(ligands)!=len(view.ligands) or len({b.binding_id for b in view.metal_bindings})!=len(view.metal_bindings):raise FormalError('Duplicate ligand/binding ID')
    for k,l in enumerate(view.ligands):
        for other in view.ligands[k+1:]:
            a,b=set(l.atom_ids),set(other.atom_ids)
            if a&b and not(a<b and other.internal_metal_ids or b<a and l.internal_metal_ids):
                raise FormalError('Overlapping ligand identities need explicit nested internal-metal representation')
    for l in ligands.values():
        ids=set(l.atom_ids)
        if len(ids)!=len(l.atom_ids) or not ids<=atoms.keys() or not set(l.donor_atom_ids)<=ids:raise FormalError('Invalid ligand membership')
        if l.ionic_class=='X3':
            bs=[b for b in view.metal_bindings if b.ligand_id==l.ligand_id]
            if any(b.metal_atom_id not in mids for b in bs):raise FormalError('X3 binding references missing metal')
            if l.ligand_type not in {'Schrock_alkylidyne','metallacyclobutadiene'} or not bs or any(atoms[b.metal_atom_id].element not in {'W','Mo'} or b.mode!='sigma_covalent' for b in bs):
                raise Unsupported('X3 scope is declared W/Mo alkylidyne or metallacyclobutadiene')
            centers_x3={b.metal_atom_id for b in bs};donors=set(l.donor_atom_ids)
            if len(centers_x3)!=1 or any(atoms[i].element!='C' for i in donors):raise Unsupported('X3 requires one W/Mo and carbon donors')
            m=next(iter(centers_x3));orders=[edges.get(tuple(sorted((m,i))),0) for i in donors]
            if l.ligand_type=='Schrock_alkylidyne' and (len(donors)!=1 or orders!=[3]):raise FormalError('Alkylidyne label requires an actual metal-carbon triple bond')
            if l.ligand_type=='metallacyclobutadiene':
                if len(donors)!=2 or sorted(orders)!=[1,2]:raise FormalError('MCBD requires one single and one double metal-carbon bond')
                a,c=sorted(donors)
                if not any(atoms[b].element=='C' and sorted([edges.get(tuple(sorted((a,b))),0),edges.get(tuple(sorted((b,c))),0)])==[1,2] for b in ids-donors):raise FormalError('MCBD requires a connected four-membered alternating ring')
        if ids&mids!=set(l.internal_metal_ids):raise FormalError('Internal ligand metal must be explicit, e.g. dppf Fe')
        if not ids or not l.donor_atom_ids:raise FormalError('Empty ligand or donor set')
        graph={i:set() for i in ids}
        for a,b in edges:
            if {a,b}<=ids:graph[a].add(b);graph[b].add(a)
        for binding in view.metal_bindings:
            for donor in binding.donor_atom_ids:
                if {binding.metal_atom_id,donor}<=ids:
                    graph[donor].add(binding.metal_atom_id);graph[binding.metal_atom_id].add(donor)
        reached={min(ids)}
        while True:
            nxt=reached|{j for i in reached for j in graph[i]}
            if nxt==reached:break
            reached=nxt
        if reached!=ids:raise FormalError('Ligand must be a connected covalent/coordination fragment')
        if l.ionic_class!='unknown' and l.ionic_charge!=CLASSES[l.ionic_class][0]+l.remote_ionic_charge and not(l.ionic_charge is None and symbolic&ids and l.ionic_class=='L'):raise FormalError('Ligand ionic class/charge mismatch')
        if l.remote_ionic_charge:
            if l.ionic_class not in {'LX','L2X'}:raise Unsupported('Remote charge partition only registered for haptic ligand spectators')
            remote=set(l.atom_ids)-set(l.donor_atom_ids)
            if sum(atoms[i].formal_charge for i in remote)!=l.remote_ionic_charge:raise FormalError('Remote charge does not match explicitly charged non-donor atoms')
        bs=[b for b in view.metal_bindings if b.ligand_id==l.ligand_id]
        if len({b.metal_atom_id for b in bs})>1 and any(b.ionic_allocation is None or b.donated_electrons is None for b in bs):
            raise Unsupported('Bridge needs explicit per-metal ionic allocation and donor count')
        if bs:
            if any(bool(a in ids)!=bool(b in ids) and not ({a,b}-ids)<=mids for a,b in edges):raise FormalError('Bound ligand omits covalently connected atoms')
            cov=sum(o for (a,b),o in edges.items() if (a in ids and b in mids-ids) or (b in ids and a in mids-ids))
            # Fischer/Ru carbene convention differs from electronegative X2 partition.
            if not neutral_ru_carbene(view,l):
                if l.ionic_charge is not None and sum(atoms[i].formal_charge for i in ids)-cov!=l.ionic_charge:
                    raise FormalError('Ligand ionic partition differs from declared charge')
    pi_reservations=set()
    for b in view.metal_bindings:
        if b.metal_atom_id not in mids or b.ligand_id not in ligands:raise FormalError('Binding references missing metal/ligand')
        l=ligands[b.ligand_id];ds=b.donor_atom_ids
        if not ds:raise FormalError('Empty coordinated donor subset')
        if b.site_id:
            if b.site_id not in centers[b.metal_atom_id].site_ids or (b.metal_atom_id,b.site_id) in sites:raise FormalError('Unknown or multiply occupied coordination site')
            sites.add((b.metal_atom_id,b.site_id))
        if len(set(ds))!=len(ds) or not set(ds)<=set(l.donor_atom_ids):raise FormalError('Invalid donor subset')
        key=(b.metal_atom_id,b.ligand_id,tuple(sorted(ds)))
        if key in seen:raise FormalError('Duplicate binding')
        seen.add(key)
        q,ne=CLASSES.get(l.ionic_class,(None,None))
        if b.pi_pool_slot is not None and b.mode!='eta':raise FormalError('Pi pool slot is only defined for haptic binding')
        if q is None:raise Unsupported('Unknown ligand electronic class')
        allocation=b.ionic_allocation if b.ionic_allocation is not None else q
        donor_e=b.donated_electrons if b.donated_electrons is not None else (2*len(ds) if b.mode=='dative' else ne)
        cov=sum(edges.get(tuple(sorted((b.metal_atom_id,i))),0) for i in ds)
        if b.mode=='dative':
            if cov or b.hapticity!=1:raise FormalError('Dative bond duplicated or hapticity confused with denticity')
            if donor_e!=2*len(ds):raise FormalError('Localized dative donors contribute one pair each')
            for i in ds:reserved[i]+=2
        elif b.mode=='eta':
            if cov or b.hapticity!=len(ds) or b.hapticity not in {2,3,5}:raise Unsupported('Eta requires distinct haptic donors, not ordinary M-C edges')
            if any(atoms[i].element!='C' for i in ds):raise Unsupported('Hetero-haptic profile not registered')
            graph={i:set() for i in ds}
            for x,y in edges:
                if x in graph and y in graph:graph[x].add(y);graph[y].add(x)
            reached={ds[0]}
            while True:
                nxt=reached|{j for i in reached for j in graph[i]}
                if nxt==reached:break
                reached=nxt
            if len(reached)!=len(ds):raise FormalError('Disconnected haptic region')
            if (b.hapticity,l.ionic_class)==(2,'LX2') and l.ligand_type=='Schrock_alkylidene':
                sig=[x for x in view.metal_bindings if x.ligand_id==l.ligand_id and x.mode=='sigma_covalent']
                if (atoms[b.metal_atom_id].element not in {'Mo','W'} or len(sig)!=1 or len(sig[0].donor_atom_ids)!=1 or
                    edges.get(tuple(sorted((b.metal_atom_id,sig[0].donor_atom_ids[0]))))!=2 or b.pi_pool_slot!=0 or
                    donor_e!=2 or allocation!=0 or edges.get(tuple(sorted(ds)))!=2):
                    raise Unsupported('LX2 eta2 scope is one Mo/W alkylidene plus its actual intramolecular C=C')
            elif (b.hapticity,l.ionic_class)==(2,'L2') and l.ligand_type=='Ru_alkylidene' and neutral_ru_carbene(view,l):
                if b.pi_pool_slot!=0 or donor_e!=2 or allocation!=0:raise FormalError('Ru intramolecular pi binding needs one explicit two-electron pool')
            elif (b.hapticity,l.ionic_class)==(2,'L2'):
                # Finite dba precursor representation: one complete neutral
                # dienone ligand, two separately owned alkene pi pools. A
                # ligand class name alone cannot authorize extra donation.
                ds_all=set(l.donor_atom_ids)
                if len(ds_all)!=4 or l.ionic_charge!=0 or b.pi_pool_slot!=0 or donor_e!=2 or atoms[b.metal_atom_id].element!='Pd':
                    raise Unsupported('L2 eta2 scope is the explicit Pd/dba two-pi-pool precursor')
                pairs={edge for edge,n in edges.items() if n==2 and set(edge)<=ds_all}
                if len(pairs)!=2 or set().union(*map(set,pairs))!=ds_all:raise FormalError('dba requires two distinct disjoint C=C donor pools')
                if not any(atoms[k].element=='C' and any(atoms[o].element=='O' and edges.get(tuple(sorted((k,o))))==2 for o in l.atom_ids)
                    and all(any(edges.get(tuple(sorted((k,a))))==1 for a in pair) for pair in pairs) for k in l.atom_ids):
                    raise Unsupported('L2 pi donor topology must be C=C-C(=O)-C=C')
            elif (b.hapticity,l.ionic_class)==(2,'LX'):
                sig=[x for x in view.metal_bindings if x.ligand_id==l.ligand_id and x.metal_atom_id==b.metal_atom_id and x.mode=='sigma_covalent']
                if l.remote_ionic_charge:
                    remote=[i for i in l.atom_ids if i not in l.donor_atom_ids and atoms[i].formal_charge]
                    if (l.remote_ionic_charge!=-1 or len(remote)!=1 or atoms[remote[0]].element!='O'
                        or atoms[remote[0]].formal_charge!=-1 or sum(o for edge,o in edges.items() if remote[0] in edge)!=1):
                        raise Unsupported('Localized allyl remote charge requires a single terminal alkoxide')
                if l.ligand_type!='allyl' or len(sig)!=1 or len(sig[0].donor_atom_ids)!=1 or l.ionic_charge!=-1+l.remote_ionic_charge or b.pi_pool_slot!=0:
                    raise Unsupported('LX eta2 restricted to localized sigma+pi allyl convention')
                a=sig[0].donor_atom_ids[0]
                if a in ds or edges.get(tuple(sorted(ds)))!=2 or not any(edges.get(tuple(sorted((a,j))))==1 for j in ds):
                    raise FormalError('Localized allyl requires metal-C-C=C topology')
            elif (b.hapticity,l.ionic_class) not in {(2,'L'),(3,'LX'),(5,'L2X')}:raise Unsupported('Unsupported haptic ionic class')
            if b.pi_pool_slot is not None:
                if b.hapticity!=2 or donor_e!=2 or allocation!=0:raise FormalError('Pi pool slot requires neutral two-electron eta2 donation')
                edge=tuple(sorted(ds));capacity=edges.get(edge,0)-1
                if b.pi_pool_slot>=capacity:raise FormalError('Requested pi pool is absent in the input bond')
                token=(*edge,b.pi_pool_slot)
                if token in pi_reservations:raise FormalError('Repeated haptic pi-pair donation')
                pi_reservations.add(token)
        elif b.mode=='sigma_covalent':
            if b.hapticity!=1 or not all(edges.get(tuple(sorted((b.metal_atom_id,i))),0) for i in ds):raise FormalError('Each sigma donor requires an existing M-ligand edge; denticity is not eta')
            if donor_e!=2*cov and not neutral_ru_carbene(view,l):raise FormalError('Ionic covalent donor count differs from bond order')
            if neutral_ru_carbene(view,l) and donor_e!=2:raise FormalError('Neutral carbene sigma depiction represents one L pair, not a second electron inventory')
            covered.update(tuple(sorted((b.metal_atom_id,i))) for i in ds)
        alloc[b.metal_atom_id]+=allocation;donations[b.metal_atom_id]+=donor_e
    for l in ligands.values():
        bs=[b for b in view.metal_bindings if b.ligand_id==l.ligand_id]
        opaque_remote=l.ionic_class=='L' and l.ionic_charge is None and bool(set(l.atom_ids)&symbolic)
        if len({b.metal_atom_id for b in bs})>1 and not opaque_remote and sum(b.ionic_allocation for b in bs)!=l.ionic_charge:
            raise FormalError('Bridge ionic allocation must sum to ligand charge')
        eta=[b for b in bs if b.mode=='eta' and b.hapticity==2]
        if len({b.metal_atom_id for b in eta})>1 and any(b.pi_pool_slot is None for b in eta):
            raise Unsupported('Bridged eta2 donation needs distinct explicit pi-pool allocations')
    for edge in edges:
        if set(edge)&mids and not set(edge)<=mids and edge not in covered:raise FormalError('Metal covalent edge lacks binding')
    from ..chemistry import bond_valences
    adduct_acceptors=set()
    for ad in view.lewis_adducts:
        d,a=ad.donor_atom_id,ad.acceptor_atom_id
        if d not in atoms or a not in atoms or a in adduct_acceptors or edges.get(tuple(sorted((d,a))),0):raise FormalError('Invalid/repeated Lewis adduct')
        if atoms[a].element!='B' or bond_valences(view)[a]!=3 or (atoms[a].formal_charge,atoms[a].nonbonding_electrons,atoms[a].radical_electrons)!=(0,0,0):raise Unsupported('Finite main-group dative acceptor is neutral tricoordinate B')
        if atoms[d].element not in {'N','O','P','S','Se'}:raise Unsupported('B adduct donor element outside finite profile')
        adduct_acceptors.add(a);reserved[d]+=2
    for i,n in reserved.items():
        a=atoms[i]
        if a.nonbonding_electrons is None or a.radical_electrons is None:raise Unsupported('Donor pair occupation unknown')
        if n>a.nonbonding_electrons-a.radical_electrons:raise FormalError('Repeated donor-pair consumption')
    rows=[]
    local = {r.atom_id:r for r in active_rules(view.view_id) if isinstance(r,ValenceRule)}
    for i,m in centers.items():
        a=atoms[i];bs=[b for b in view.metal_bindings if b.metal_atom_id==i]
        if a.element=='Ca' and rules is LEGACY:
            if any(b.mode!='dative' or len(b.donor_atom_ids)!=1 or any(atoms[j].element not in {'N','O'} for j in b.donor_atom_ids) for b in bs):
                raise Unsupported('Ca(II) coordination is limited to donor-owned O/N pairs; no eta, bridging or covalent model')
            if any(sum(x.ligand_id==b.ligand_id for x in view.metal_bindings)>1 for b in bs):
                raise Unsupported('Ca bridging/chelation requires a separate scoped profile')
        if len(set(m.site_ids))!=len(m.site_ids):raise FormalError('Duplicate metal site')
        if m.vacant_site_ids is not None and (not set(m.vacant_site_ids)<=set(m.site_ids) or any((i,s) in sites for s in m.vacant_site_ids)):raise FormalError('Vacant site conflicts with occupied site')
        own={b.binding_id for b in bs}
        if any(p.first_binding_id not in own or p.second_binding_id not in own or p.first_binding_id==p.second_binding_id for p in m.relative_positions):raise FormalError('Invalid geometry binding reference')
        if (m.spin_state=='closed_shell' and a.radical_electrons not in {0,None}) or (m.spin_state=='open_shell' and a.radical_electrons==0):raise FormalError('Local metal spin assertion conflicts with occupation')
        # Covalent partition: each ligand-to-metal ownership shift increases OS;
        # ionic dative ligands retain formal charge and do not shift metal charge.
        shift=0
        for b in bs:
            l=ligands[b.ligand_id]
            if b.mode=='sigma_covalent' and not neutral_ru_carbene(view,l):
                shift+=sum(edges.get(tuple(sorted((i,j))),0) for j in b.donor_atom_ids)
        os=a.formal_charge+shift
        cn=sum(len(b.donor_atom_ids) if b.mode!='eta' else 1 for b in bs)
        maximum=12 if a.element in {'La','Ce','Sm','Ba'} else 8
        rule=local.get(i)
        if rule is not None and rule.occupied_sites is not None:
            require_equal(rule,'occupied_sites',rule.occupied_sites,cn)
        elif cn>maximum:raise Unsupported('Coordination sites exceed finite v0.9 capacity')
        if rule is not None and rule.oxidation_state is not None:
            require_equal(rule,'oxidation_state',rule.oxidation_state,os)
        elif os not in OXIDATION.get(a.element,set()):raise Unsupported(f'{a.element}({os}) finite oxidation profile not registered')
        if m.oxidation_state is not None and m.oxidation_state!=os:raise FormalError('State-local metal oxidation declaration disagrees with ownership')
        d, residual = rules.metal_counts(a.element, os)
        ve=None if residual is None else residual+donations[i]+sum(o for edge,o in edges.items() if i in edge and set(edge)<=mids)
        for key,value in [('d_electrons',d),('valence_electrons',ve)]:
            if getattr(m,key) is not None:
                if value is None:raise Unsupported('Declared metal count cannot be independently derived: '+key)
                if getattr(m,key)!=value:raise FormalError('Metal count declaration mismatch: '+key)
        rows.append({'metal_atom_id':i,'element':a.element,'oxidation_state':os,'atomic_formal_charge':a.formal_charge,
            'd_electrons':d,'valence_electrons':ve,'coordination_number':sum(len(b.donor_atom_ids) for b in bs),'occupied_sites':cn,
            'eighteen_electron_rule':'advisory_only','spin':'unknown' if a.radical_electrons is None else 'local_occupation_declared'})
    return reserved,rows
