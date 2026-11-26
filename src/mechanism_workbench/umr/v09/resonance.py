import sys
from ..v07.resonance import candidates as old_candidates
from ..v07.transaction import touched
from .elements import current_elements
from ..v04.resonance import resonance_scan as scan,action_key
from ..v04.ledger import electronic_signature
from ..chemistry import FormalError,Unsupported,identity,hydrogen_connections,bond_orders
from .ledger import validate_view,signature,compare_views
from .execution import transact

def product_spin(origin,event,target):
    from .spin_events import SPECS,product_spin as declared_transition
    if event.general_event in SPECS:return declared_transition(origin,event,target)
    from ..v05.execution import product_spin as legacy
    from ..v07.ledger import symbols,project_symbols
    if symbols(target):
        mult,report=legacy(origin,event,project_symbols(target))
        report.update(scope='known fragment with unchanged opaque boundaries',global_spin_verified=False)
        return mult,report
    return legacy(origin,event,target)
def protected_region(v):
    mids={a.atom_id for a in v.atoms if a.element in current_elements().metals}
    return mids|{i for b in v.bonds if {b.a,b.b}&mids for i in (b.a,b.b)}|{i for b in v.metal_bindings for i in b.donor_atom_ids}
def candidates(v):
    yield from mcbd_candidates(v)
    shadow=v.model_copy(deep=True);mids={a.atom_id for a in v.atoms if a.element in current_elements().metals}
    shadow.atoms=[a for a in shadow.atoms if a.atom_id not in mids]
    shadow.bonds=[b for b in shadow.bonds if not {b.a,b.b}&mids]
    shadow.metal_centers=[];shadow.metal_bindings=[];shadow.ligands=[]
    shadow.element_declarations=[d for d in shadow.element_declarations if d.atom_id not in mids]
    if not shadow.atoms:return
    for rule,mode,actions in old_candidates(shadow):
        if not touched(actions)&protected_region(v):yield rule,mode,actions
def mcbd_candidates(v):
    from ..ast import Shift
    from ..chemistry import atoms_by_id,order
    from .coordination import checks
    counts={c['metal_atom_id']:c for c in checks(v)[1]}
    for l in v.ligands:
        if l.ligand_type!='metallacyclobutadiene' or l.ionic_class!='X3':continue
        for m,row in counts.items():
            if row['element'] not in {'W','Mo'} or row['oxidation_state']!=6:continue
            for a in l.donor_atom_ids:
                for c in l.donor_atom_ids:
                    if a==c or order(v,m,a)!=2 or order(v,m,c)!=1:continue
                    for b in l.atom_ids:
                        if order(v,a,b)==1 and order(v,b,c)==2 and all(atoms_by_id(v)[i].element=='C' for i in [a,b,c]):
                            yield 'mcbd_pi_redistribution_r10','POLAR',[Shift(a=a,b=m,c=c),Shift(a=c,b=b,c=a)]

def apply(v,actions,target_id):
    from .metal_events import repartition
    specific=any(action_key(actions)==action_key(a) for _,_,a in mcbd_candidates(v))
    return transact(v,actions,v.state_id,target_id,preserve_stereo=True,after=repartition if specific else None)
def validate_resonance(s,t,r):
    if (r.source_view_id,r.target_view_id,r.state_id)!=(s.view_id,t.view_id,s.state_id):raise FormalError('Resonance references differ from supplied views')
    if electronic_signature(s)!=electronic_signature(t):raise Unsupported('Electronic-state change is not resonance')
    from ..stereochemistry import check_resonance_stereo
    check_resonance_stereo(s,t)
    if s.state_id!=t.state_id or identity(s)!=identity(t) or hydrogen_connections(s)!=hydrogen_connections(t) or set(bond_orders(s))!=set(bond_orders(t)):raise FormalError('Resonance identity/sigma/H changed')
    if not any(rule==r.rule and mode==r.resonance_mode and action_key(actions)==action_key(r.actions) for rule,mode,actions in candidates(s)):raise Unsupported('unsupported_for_resonance')
    compare_views(apply(s,r.actions,t.view_id),t)
def resonance_scan(view,*,selection_policy=None,scope_atom_ids=None,**kw):
    from types import SimpleNamespace
    from dataclasses import asdict
    from .resonance_policy import ResonancePolicy,decision,VERSION
    policy=selection_policy or ResonancePolicy();suppressed=[];outside_scope=0
    scope=None if scope_atom_ids is None else set(scope_atom_ids)
    if scope is not None and not scope <= {a.atom_id for a in view.atoms}:
        raise FormalError('Resonance scope references missing atoms')
    def selected(current):
        nonlocal outside_scope
        for rule,mode,actions in candidates(current):
            if scope is not None and not touched(actions)<=scope:
                outside_scope+=1;continue
            try:target=apply(current,actions,'selection-preview')
            except (FormalError,Unsupported):continue
            check=decision(current,target,policy)
            if check is not None and not check['selected']:
                suppressed.append({'rule':rule,**check});continue
            yield rule,mode,actions
    # max_views includes the input view. Previously max_views=1 still ran
    # potentially hundreds of expensive preview transactions before noticing
    # it could not retain a single additional contributor (notably PPh3 rings).
    # Validate the state as usual and report an explicit exhausted budget;
    # never misreport this as an exhaustive "no resonance" result.
    exhausted_view_budget=kw.get('max_views')==1
    backend=SimpleNamespace(validate_view=validate_view,protected_region=protected_region,
        candidates=(lambda _:iter(())) if exhausted_view_budget else selected,apply=apply,signature=signature,validate_resonance=validate_resonance)
    result=scan(view,_backend=backend,**kw);result[2].rule_version='umr-resonance/0.9-r2'
    if exhausted_view_budget and result[2].output.get('reason') is None:
        result[2].status='truncated';result[2].rule_version='umr-resonance/0.9-r3'
        result[2].output.update(truncated=True,exhaustive=False,expanded_views=0,
            truncation_reason='max_views_includes_input',candidate_transactions=0)
    result[2].input['selection_policy']={'version':VERSION,**asdict(policy)}
    result[2].output.update(selection_suppressed_count=len(suppressed),selection_suppressed=suppressed,
        selection_is_formal_prohibition=False)
    if scope is not None:
        result[2].input['explicit_search_scope_atom_ids']=sorted(scope)
        result[2].output['outside_search_scope_candidates']=outside_scope
        result[2].uncertainties.append('Explicit search scope omits other regions; omitted species and their electrons remain in the full validated state.')
    result[2].uncertainties.append('Neutral C/C charge separation is selected by declared electronegativity evidence; no site values inferred from partial charges or substituent constants.')
    for link in result[1]:link.rule_version='umr-resonance/0.9'
    return result

def prepare_source(source,event):
    from ..v06.execution import prepare_source as legacy
    from ..chemistry import order
    if event.general_event=='ALKYLIDYNE_CYCLOREVERSION':
        m,c=event.general_roles['metal'],event.general_roles['first']
        if order(source,m,c)==1:return source,None
        from .ast import ResonanceTransform
        for rule,mode,actions in mcbd_candidates(source):
            target=apply(source,actions,source.view_id+'.mcbd.'+event.event_id)
            if order(target,m,c)!=1:continue
            link=ResonanceTransform(transform_id='mcbd.'+event.event_id,state_id=source.state_id,source_view_id=source.view_id,
                target_view_id=target.view_id,actions=actions,rule=rule,resonance_mode=mode,rule_version='umr-resonance/0.9')
            validate_resonance(source,target,link);return target,link
        return source,None
    if event.reactive_template=='ARYNE_CONCERTED_ELIMINATION':
        trial=event.model_copy(deep=True);trial.reactive_template='ARYNE_DEPROTONATION'
        return legacy(source,trial)
    if event.general_event!='GASSMAN_2_3':return legacy(source,event)
    d,e=event.general_roles.get('d'),event.general_roles.get('e')
    if order(source,d,e)==2:return source,None
    from .ast import ResonanceTransform
    for rule,mode,actions in candidates(source):
        if not rule.startswith('kekule_'):continue
        target=apply(source,actions,source.view_id+'.kekule.'+event.event_id)
        if order(target,d,e)!=2:continue
        link=ResonanceTransform(transform_id='kekule.'+event.event_id,state_id=source.state_id,source_view_id=source.view_id,
          target_view_id=target.view_id,actions=actions,rule=rule,resonance_mode=mode,rule_version='umr-resonance/0.9')
        validate_resonance(source,target,link);return target,link
    return source,None
