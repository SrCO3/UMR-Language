import json
from collections import Counter
from ..chemistry import FormalError,Unsupported,atoms_by_id,identity,electron_count,bond_orders
from ..v07.execution import execute_event as old
from ..v08.execution import execute_event as old08
from ..v07.transaction import transact as transaction,touched
from ..v06.stereo import apply_event
from .ledger import validate_view
from .elements import charge, current_elements
from .coordination import checks

def transact(source,actions,state_id,view_id,*,preserve_stereo=False,after=None,_released_binding_ids=frozenset(),_declared_pair_donors=frozenset()):
    source=current_elements().view_type.model_validate(source.model_dump(mode='json'))
    if touched(actions)&{i for b in source.metal_bindings if b.mode=='eta' and b.binding_id not in _released_binding_ids for i in b.donor_atom_ids}:
        raise Unsupported('Haptic electron redistribution requires a dedicated coordinated profile')
    def state_local(v):
        if after:after(v)
        if not preserve_stereo and source.electronic_state.character=='open_shell_singlet' and not any(a.radical_electrons for a in v.atoms):
            # Recombination of an already declared singlet is not ISC.
            v.electronic_state.character='closed_shell';v.electronic_state.coupling='closed_shell'
            v.electronic_state.label='closed-shell recombination product; preceding singlet declaration retained in evidence'
        from .ast import Ligand,MetalBinding
        v.ligands=[Ligand.model_validate(l.model_dump(mode='json')) for l in v.ligands]
        v.metal_bindings=[MetalBinding.model_validate(b.model_dump(mode='json')) for b in v.metal_bindings]
        if not preserve_stereo:
            # Free components are no longer coordinated ligands. Remove only
            # stale classification metadata; retain every atom and bond.
            bound_ids={b.ligand_id for b in v.metal_bindings}
            v.ligands=[l for l in v.ligands if l.ligand_id in bound_ids or not set(l.atom_ids)&touched(actions)]
            # Refresh an ordinary ligand fragment only when the executed graph
            # changed its membership. Haptic/bridged partitions remain explicit.
            from .metal_events import repartition
            mids={a.atom_id for a in v.atoms if a.element in current_elements().metals};bound={b.ligand_id for b in v.metal_bindings};stale=False
            for l in v.ligands:
                if l.ligand_id not in bound:continue
                ids=set(l.atom_ids);links={i:set() for i in ids}
                for b in v.bonds:
                    if {b.a,b.b}<=ids:links[b.a].add(b.b);links[b.b].add(b.a)
                    if bool(b.a in ids)!=bool(b.b in ids) and not ({b.a,b.b}-ids)<=mids:stale=True
                seen={min(ids)} if ids else set()
                while seen:
                    new=seen|{j for i in seen for j in links[i]}
                    if new==seen:break
                    seen=new
                if seen!=ids and not l.internal_metal_ids:stale=True
                if not l.internal_metal_ids and l.ionic_charge is not None:
                    atoms=atoms_by_id(v);cov=sum(b.order for b in v.bonds if bool(b.a in ids)!=bool(b.b in ids) and ({b.a,b.b}-ids)<=mids)
                    if sum(atoms[i].formal_charge for i in ids)-cov!=l.ionic_charge and l.ligand_type not in {'Fischer_carbene','Ru_alkylidene'}:stale=True
            if stale:repartition(v)
            for d in v.element_declarations:
                if d.atom_id in touched(actions):d.oxidation_state=None;d.coordination_number=None
    def reservations(v):
        from ..v07.ledger import symbols
        reserved,counts=checks(v,_symbolic_ids=symbols(source))
        for b in source.metal_bindings:
            if b.binding_id in _released_binding_ids and b.mode=='dative':
                for i in b.donor_atom_ids:reserved[i]-=2
        return reserved,counts
    return transaction(source,actions,state_id,view_id,validator=validate_view,charge_fn=charge,
        coordination_fn=reservations,preserve_stereo=preserve_stereo,after=state_local,
        _unknown_pair_receivers={a.atom_id for a in source.atoms if a.element in current_elements().metals},
        _declared_pair_donors=_declared_pair_donors)

def execute_event(source,event,*,waivers=frozenset(),_dispatcher=None):
    from .face_branches import KEY, apply
    target,report=(_dispatcher or _execute_event)(source,event,waivers=waivers)
    from .stereo_path import reconcile
    report['stereo_retention_checks']=reconcile(source,target,event)
    if KEY in event.conditions:
        if waivers or report.get('execution_mode')!='electron_actions':
            raise FormalError('Local geometric hypotheses require a strict electronic execution')
        apply(source,target,event,report)
    if 'peracid_syn_centers_v1' in event.conditions:
        from .syn_epoxide import check,KEY as SYN
        report['syn_epoxide_geometry']=check(source,target,event.conditions[SYN])
    return target,report


def _execute_event(source,event,*,waivers=frozenset()):
    if source.view_id!=event.reactant_view_id or source.state_id==event.product_state_id:raise FormalError('Event state/view references invalid')
    from .dienophile import PROFILES as dienophiles, execute as execute_dienophile
    if event.template_id in dienophiles:
        if waivers:raise Unsupported('Heterodienophile role constraints cannot be waived')
        return execute_dienophile(source,event)
    from .sigmatropic import PROFILES as aza_profiles,VERSION as aza_version,execute as execute_aza
    # Historical missing-profile probes stay unsupported until explicitly migrated.
    if event.template_id in aza_profiles and event.template_version==aza_version:
        if waivers:raise Unsupported('Aza topology/occupation checks cannot be waived')
        return execute_aza(source,event)
    active=touched(event.actions or [])|set(event.general_roles.values())|{i for t in event.electron_transfers for i in [t.donor_atom_id,t.acceptor_atom_id]}
    for boundary in source.material_boundaries:
        if active&set(boundary.atom_ids) and (boundary.direction!='inside' or boundary.domain in {'bulk_metal','surface','external_electron_reservoir'}):
            raise Unsupported('External/bulk/surface participant needs a dedicated exchange model; not an isolated molecular atom')
    if event.general_event == 'AUTHORED_ELECTRON_ACTIONS':
        from .authored import execute
        if waivers:raise Unsupported('Authored transactions require the native checks')
        return execute(source,event,transact)
    if event.general_event == 'PILOT_VINYL_TIN_I2_RETENTIVE_V1':
        from .vinyl_tin import execute
        if waivers: raise Unsupported('Vinyl tin pilot requires all structural checks')
        return execute(source,event)
    if event.general_event=='BROMONIUM_OPEN_BROMIDE':
        from .bromonium import opening
        return opening(source,event,waivers)
    from .silyl import SPECS as silyl_specs,execute as execute_silyl
    if event.general_event in silyl_specs:
        if waivers:raise Unsupported('Silyl profile occupancy/topology checks cannot be waived')
        return execute_silyl(source,event)
    from .halonium import PROFILE as halonium_profile,execute as execute_halonium
    if event.general_event==halonium_profile:
        if waivers:raise Unsupported('Halonium occupation and topology cannot be waived')
        return execute_halonium(source,event)
    from .aryne import SPECS as aryne_specs,execute as execute_aryne
    from .pcet import SPECS as pcet_specs,execute as execute_pcet
    if event.general_event in pcet_specs:
        if waivers:raise Unsupported('PCET slot/ledger checks cannot be waived')
        return execute_pcet(source,event)
    from .cluster_events import SPECS as cluster_specs,execute as execute_cluster
    if event.general_event in cluster_specs:
        if waivers:raise Unsupported('Cluster pool allocation cannot be waived')
        return execute_cluster(source,event)
    from .spin_events import SPECS as spin_specs,execute as execute_spin
    if event.general_event in spin_specs:
        if waivers:raise Unsupported('Electronic state transition constraints cannot be waived')
        return execute_spin(source,event)
    if event.element_template in {'BORON_L_ASSOCIATION','BORON_L_DISSOCIATION'}:
        if waivers:raise Unsupported('B coordination ledger/profile cannot be waived')
        from .boron import execute as execute_boron
        return execute_boron(source,event)
    if event.reactive_template in aryne_specs:
        if waivers:raise Unsupported('Aryne ledger/profile checks cannot be waived')
        return execute_aryne(source,event)
    from . import polar_metal
    from .metathesis import SPECS as om_specs,execute as execute_om
    if event.general_event in om_specs:
        if waivers:raise Unsupported('Metathesis occupancy and ligand scope cannot be waived')
        return execute_om(source,event)
    from .cross_coupling import SPECS as coupling_specs,execute as execute_coupling
    if event.general_event in coupling_specs:
        if waivers:raise Unsupported('Coupling ledger and element constraints cannot be waived')
        return execute_coupling(source,event)
    from .catalytic_profiles import SPECS as catalytic_specs,execute as execute_catalytic
    if event.general_event in catalytic_specs:
        if waivers:raise Unsupported('Catalytic profile topology and electron ledger cannot be waived')
        return execute_catalytic(source,event)
    if event.general_event in polar_metal.SPECS:
        if waivers:raise Unsupported('Finite transaction ledger/profile cannot be waived')
        return polar_metal.execute(source,event,event.general_event,event.general_roles)
    if event.element_template=='METAL_X_HETEROLYSIS':
        if waivers:raise Unsupported('Finite transaction ledger/profile cannot be waived')
        roles={'m':event.bindings.get('metal'),'x':event.bindings.get('leaving_group')}
        return polar_metal.execute(source,event,'SIGMA_X_HETEROLYSIS',roles)
    if not any([event.general_event,event.metal_template,event.rule_template,event.element_template,event.template_id,event.reactive_template]) and polar_metal.candidates(source,event):
        if waivers:raise Unsupported('Finite metal transaction checks cannot be waived')
        return polar_metal.execute(source,event)
    if event.general_event=='METAL_SET':
        from .redox import execute
        return execute(source,event)
    if event.general_event=='DISSOCIATIVE_METAL_SET':
        from .redox import execute_dissociative
        return execute_dissociative(source,event)
    if event.general_event in {'MONONUCLEAR_OA_2E','MONONUCLEAR_RE_2E'}:
        expected='OXIDATIVE_ADDITION_2E' if event.general_event=='MONONUCLEAR_OA_2E' else 'REDUCTIVE_ELIMINATION_2E'
        if event.metal_template!=expected:raise FormalError('Explicit general/metal template pair required')
    if event.electron_transfers:raise FormalError('Ownership transfer requires its explicit METAL_SET profile')
    if event.rule_template in {'LEWIS_ASSOCIATION','LEWIS_DISSOCIATION'}:
        validate_view(source);target=source.model_copy(deep=True)
        m=event.bindings.get('metal');donor=event.bindings.get('donor')
        if m not in {c.atom_id for c in source.metal_centers} or donor not in atoms_by_id(source):raise FormalError('Coordination roles missing')
        if event.actions:raise FormalError('Donor-owned coordination must not consume its pair as a covalent action')
        matches=[b for b in target.metal_bindings if b.metal_atom_id==m and donor in b.donor_atom_ids and b.mode=='dative']
        if event.rule_template=='LEWIS_DISSOCIATION':
            if len(matches)!=1:raise FormalError('Expected one existing donor binding')
            b=matches[0];b.donor_atom_ids.remove(donor)
            if not b.donor_atom_ids:target.metal_bindings.remove(b)
        else:
            if matches:raise FormalError('Donor already coordinated')
            ligs=[l for l in target.ligands if donor in l.donor_atom_ids]
            if not ligs:
                from .ast import Ligand
                ats=atoms_by_id(target);mids={i for i,a in ats.items() if a.element in current_elements().metals}
                ids={donor};stack=[donor]
                while stack:
                    i=stack.pop()
                    for b in target.bonds:
                        if i not in {b.a,b.b}:continue
                        j=b.b if b.a==i else b.a
                        if j not in mids|ids:ids.add(j);stack.append(j)
                charges=[ats[i].formal_charge for i in ids]
                q=None if any(x is None for x in charges) else sum(charges)
                if q not in {None,0} or ats[donor].formal_charge!=0:raise Unsupported('Implicit ligand registration limited to explicit neutral donor/local opaque boundary')
                lig=Ligand(ligand_id='incoming.'+event.event_id,atom_ids=sorted(ids),donor_atom_ids=[donor],ionic_class='L',ionic_charge=q,
                   assignment_basis='complete source component; opaque charge unknown if present')
                target.ligands.append(lig);ligs=[lig]
            if len(ligs)!=1:raise FormalError('Association needs one unambiguous complete ligand definition')
            from .ast import MetalBinding
            target.metal_bindings.append(MetalBinding(binding_id=event.event_id+'.binding',metal_atom_id=m,ligand_id=ligs[0].ligand_id,donor_atom_ids=[donor],mode='dative'))
        target.state_id=event.product_state_id;target.view_id=event.product_view_id
        for c in target.metal_centers:
            if c.atom_id==m:c.valence_electrons=None;c.relative_positions=[];c.vacant_site_ids=None
        validate_view(target)
        return target,{'event_id':event.event_id,'event_type':event.event_type,'atomic_commit':True,'formal_status':'passed','execution_mode':'electron_actions','compiled_actions':[],
            'coordination_change':{'donor':donor,'metal':m,'event':event.rule_template},'electron_ownership_changed':False,'chemical_feasibility':'unknown'}
    if not event.metal_template and not event.general_event:
        # Preserve exact legacy semantics unless an explicit v0.9 role profile is requested.
        target,report=old08(source,event,waivers=waivers,_hooks={'validate_view':validate_view,'transact':transact})
        report.setdefault('execution_mode',event.execution_level)
        if source.electronic_state.character=='open_shell_singlet' and not any(a.radical_electrons for a in target.atoms):
            target.electronic_state.character='closed_shell';target.electronic_state.coupling='closed_shell'
        return target,report
    validate_view(source)
    if source.view_id!=event.reactant_view_id or source.state_id==event.product_state_id:raise FormalError('Event state/view references invalid')
    if event.metal_template:
        from ..v03.execution import execute_event as metal
        def reserved(v,actions):
            rs,_=checks(v);ats=atoms_by_id(v)
            spent=Counter(a.a for a in actions if a.op=='PAIR')
            for i,n in spent.items():
                r=ats[i].radical_electrons
                if r is None:
                    count=next((c for c in checks(v)[1] if c['metal_atom_id']==i),{})
                    if event.general_event!='MONONUCLEAR_OA_2E' or count.get('d_electrons')!=10:
                        raise Unsupported('Two-electron donor pair occupation unknown; SET is a separate profile')
                    # A declared d10 ionic counting channel has a filled formal
                    # d shell; this is a template assumption, not a total-spin assignment.
                    r=0
                if 2*n+rs[i]>ats[i].nonbonding_electrons-r:raise FormalError('Repeated consumption of donor pair')
        hooks={'validate_view':validate_view,'charge':charge,'metal_counts':lambda v:checks(v)[1],
            'binding_checks':lambda v:({m.atom_id:m for m in v.metal_centers},{l.ligand_id:l for l in v.ligands},checks(v)[0]),'check_reserved':reserved}
        if event.general_event in {'MONONUCLEAR_OA_2E','MONONUCLEAR_RE_2E'}:
            hooks['metal_scope']={('OXIDATIVE_ADDITION_2E','metal'):{('Pd',0),('Ni',0),('Cu',1)},
                ('REDUCTIVE_ELIMINATION_2E','metal'):{('Pd',2),('Ni',2),('Cu',3)}}
        transaction_source=source
        unknown_d10=set()
        if event.general_event=='MONONUCLEAR_OA_2E':
            donor=event.bindings.get('metal');counts={c['metal_atom_id']:c for c in checks(source)[1]}
            if donor in counts and counts[donor]['d_electrons']==10 and atoms_by_id(source)[donor].radical_electrons is None:
                transaction_source=source.model_copy(deep=True);atoms_by_id(transaction_source)[donor].radical_electrons=0;unknown_d10.add(donor)
        target,report=metal(transaction_source,event,_hooks=hooks)
        for i in unknown_d10:atoms_by_id(target)[i].radical_electrons=None
        if unknown_d10:report['ledger']=validate_view(target)
        if event.general_event=='MONONUCLEAR_OA_2E':report['pair_occupation_assumption']='filled formal d10 donor in declared two-electron channel; global spin unchanged'
    else:
        from .events import compile
        actions,fragment=compile(source,event)
        from .metal_events import released_bindings
        released=released_bindings(source,event)
        key=lambda a:json.dumps(a.model_dump(mode='json'),sort_keys=True)
        if event.actions is not None and Counter(map(key,event.actions))!=Counter(map(key,actions)):raise FormalError('Actions disagree with registered profile')
        def post(v):
            from .textbook_events import postprocess
            postprocess(v,event)
            from .metal_events import SPECS as METAL_SPECS,repartition
            if event.general_event in METAL_SPECS:
                v.metal_bindings=[b for b in v.metal_bindings if b.binding_id not in released]
                repartition(v)
                if event.general_event=='RH_ACYL_DECARBONYLATION':
                    from .ast import Ligand,MetalBinding
                    r=event.general_roles;lid=event.event_id+'.retained_CO'
                    v.ligands.append(Ligand(ligand_id=lid,atom_ids=[r['carbonyl'],r['oxygen']],donor_atom_ids=[r['carbonyl']],ionic_class='L',ionic_charge=0,ligand_type='CO',assignment_basis='source-declared acyl deinsertion leaves C-bound CO'))
                    v.metal_bindings.append(MetalBinding(binding_id=lid+'.binding',ligand_id=lid,metal_atom_id=r['metal'],donor_atom_ids=[r['carbonyl']],mode='dative',donated_electrons=2))
                if event.general_event=='METALLACYCLE_CLEAVAGE':
                    for l in v.ligands:
                        if l.donor_atom_ids==[event.general_roles['third']]:l.ligand_type='Ru_alkylidene' if atoms_by_id(v)[event.general_roles['metal']].element=='Ru' else 'Schrock_alkylidene'
            if event.electronic_outcome:
                if source.electronic_state.coupling=='high_spin' and event.electronic_outcome.coupling=='closed_shell':raise Unsupported('requires_spin_transition')
                if not event.electronic_outcome.evidence:raise FormalError('Electronic outcome requires evidence')
                v.electronic_state=event.electronic_outcome.model_copy(deep=True)
            for d in v.element_declarations:
                if d.atom_id in touched(actions):d.oxidation_state=None;d.coordination_number=None
            for m in v.metal_centers:
                if m.atom_id in touched(actions):m.oxidation_state=None;m.d_electrons=None;m.valence_electrons=None
            if event.general_event=='DIAZO_N2_EXTRUSION':
                from ..v04.ast import ReactiveCenter
                v.reactive_centers.append(ReactiveCenter(center_id=event.event_id+'.carbene',atom_id=event.general_roles['carbon'],species='carbene',profile='carbene_paired',environment='free',existence_basis='formal_model',reaction_context=event.explanation or event.mechanism_assumption,evidence=event.evidence))
            if event.general_event in {'CARBENE_1_2_H_MIGRATION','ALKYLIDENECARBENE_1_2_H_MIGRATION'}:v.reactive_centers=[c for c in v.reactive_centers if c.atom_id!=event.general_roles['center']]
        target=transact(source,actions,event.product_state_id,event.product_view_id,after=post,_released_binding_ids=released)
        event.actions=actions
        if fragment:
            ids=set(fragment)
            if any(bool(b.a in ids)!=bool(b.b in ids) for b in target.bonds):raise FormalError('Released fragment remains bonded to inventory')
        report={'event_id':event.event_id,'event_type':event.event_type,'general_event':event.general_event,'released_fragments':[fragment] if fragment else [],
                'compiled_actions':[a.model_dump(mode='json') for a in actions],'stereochemistry':apply_event(source,target,event)}
    from ..v07.ledger import project_symbols
    s,t=project_symbols(source),project_symbols(target)
    if identity(source)!=identity(target) or electron_count(s)!=electron_count(t) or sum(a.formal_charge for a in s.atoms)!=sum(a.formal_charge for a in t.atoms):raise FormalError('Whole transaction conservation failure')
    from .ast import RULE_VERSION
    report.update(rule_version=RULE_VERSION,execution_mode='electron_actions',atomic_commit=True,formal_status='passed',chemical_feasibility='unknown')
    return target,report
