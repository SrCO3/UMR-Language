from ..ast import Pair, Split
from ..chemistry import FormalError, Unsupported, atoms_by_id, order, bond_valences
from ..v04.templates import compile_template
from ..v04.execution import execute_event as reactive_execute
from ..v05.execution import execute_event as execute05, make_state, product_spin
from .ledger import validate_view, validate_ledger
from ..v05.elements import charge
from .ast import ResonanceTransform
from . import stereo


def prepare_source(source,event):
    if event.reactive_template not in {'ARYNE_DEPROTONATION','ARYNE_ELIMINATION'}:return source,None
    roles=event.bindings
    if not all(k in roles for k in ('alpha_carbon','beta_carbon','ring')):raise FormalError('Aryne template lacks ring roles')
    if order(source,roles['alpha_carbon'],roles['beta_carbon'])==2:return source,None
    from .resonance import kekule_actions, apply, validate_resonance
    actions=kekule_actions(source,roles['ring'])
    target=apply(source,actions,source.view_id+'.kekule.'+event.event_id)
    if order(target,roles['alpha_carbon'],roles['beta_carbon'])!=2:raise FormalError('Aryne sites must be adjacent')
    link=ResonanceTransform(transform_id='kekule.'+event.event_id,state_id=source.state_id,
        source_view_id=source.view_id,target_view_id=target.view_id,actions=actions,
        rule='kekule_six_carbon',resonance_mode='POLAR',rule_version='umr-resonance/0.6')
    validate_resonance(source,target,link)
    return target,link


def check_applicability(source,event):
    atoms=atoms_by_id(source);vals=bond_valences(source)
    if event.event_type=='SN2':
        c=event.bindings.get('electrophile')
        if c not in atoms or atoms[c].element!='C' or vals[c]!=4 or any(b.order!=1 for b in source.bonds if c in {b.a,b.b}):
            raise FormalError('Ordinary SN2 requires the directly substituted sp3 carbon; aryl/vinyl/alkynyl center is outside scope')
        if sum(atoms[b.b if b.a==c else b.a].element=='C' for b in source.bonds if c in {b.a,b.b})>2:
            raise Unsupported('Ordinary SN2 template excludes tertiary carbon; empirical steric parameters do not override scope')
    if event.event_type=='POLAR' and not event.reactive_template:
        pairs=[a for a in event.actions or [] if isinstance(a,Pair)]
        splits=[a for a in event.actions or [] if isinstance(a,Split)]
        if len(pairs)==len(splits)==1 and len(event.actions)==2:
            p,s=pairs[0],splits[0];c=p.b
            if c==s.a and order(source,c,s.b)==1 and atoms[c].element=='C' and any(b.order>1 for b in source.bonds if c in {b.a,b.b}):
                if event.foundation_template!='SNAR_CONCERTED':
                    raise Unsupported('POLAR net substitution at sp2/sp center requires a dedicated applicability template; no generic SN2 bypass')
                from ...physical.foundation import check_input
                check_input('SNAR_CONCERTED',source,{'ipso':c,'leaving':s.b,'nucleophile':p.a})
                if not event.evidence or not event.mechanism_assumption:
                    raise Unsupported('Concerted SNAr requires source/context and an explicit mechanism assumption')
    missing=[k for k in ('solvent','temperature_kelvin') if not event.conditions.get(k)]
    return {'status':'unknown_mechanism_scope' if event.event_type=='POLAR' and not (event.foundation_template or event.reactive_template or event.element_template) else 'structural_scope_passed','missing_context':missing,'phase':event.phase,
            'physical_evidence':'unknown','chemical_plausibility':'unknown','reference_agreement':'not_evaluated'}


def compile_reactive(source,event):
    if any(a.radical_electrons is None for a in source.atoms):
        raise Unsupported('Reactive-intermediate execution requires resolved occupation; unknown spectator spin is preserved, not inferred')
    result=compile_template(source,event,_aryne_halides=frozenset({'Cl','Br'}),_allow_spectator_metals=True)
    if event.reactive_template in {'ARYNE_DEPROTONATION','ARYNE_ELIMINATION'}:
        result[-1].update(rule_version='umr-reactive/0.6',scope='isolated six-carbon aryl chloride/bromide; explicit adjacent H/NH2- and formal stepwise elimination model',
            source='https://doi.org/10.1021/ja01109a523',
            additional_context='https://archive.nptel.ac.in/content/storage2/courses/104103071/module7/lec19/3.html',
            selectivity='unknown; both aryne endpoints may be attacked; intermediate lifetime not established')
    return result


def execute_event(source,event):
    applicability=check_applicability(source,event)
    if event.reactive_template:
        target,report=reactive_execute(source,event,_compiler=compile_reactive,_validator=validate_view,
            _ledger=validate_ledger,_charge=charge)
    else:target,report=execute05(source,event,_validator=validate_view)
    report['stereochemistry']=stereo.apply_event(source,target,event)
    validate_view(target)
    report.update(applicability=applicability,representation='supported',electron_ledger='passed',
        chemical_evidence='unknown',reference_agreement='not_evaluated',rule_pack='umr/0.6')
    return target,report
