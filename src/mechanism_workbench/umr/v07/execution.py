import json
from collections import Counter
from ..chemistry import FormalError,Unsupported,atoms_by_id,order,bond_valences
from ..engine import check_event_rule
from ..v02.radical import compile_radical
from ..v02.pericyclic import compile_pericyclic
from ..v04.execution import execute_event as execute_reactive
from ..v06.execution import compile_reactive,prepare_source
from ..v06 import stereo
from .ledger import validate_view,validate_ledger
from .transaction import transact,touched
from .elements import METAL_GROUP,charge
from .substitution import FAMILIES,compile_stage,aromatic_ring
from .templates import SPECS,compile_template

WAIVABLE={'substitution.unsaturated_scope','snar.ring_scope','snar.leaving_scope',
          'template.unregistered','template.legacy_scope','geometry.template_scope'}


class RuleUnsupported(Unsupported):
    def __init__(self,rule,message):
        super().__init__(message);self.report={'rule_id':rule,'rule_version':'umr-rules/0.7-r1'}


def execute_event(source,event,*,waivers=frozenset(),_hooks=None):
    hooks=_hooks or {}
    validate_view=hooks.get('validate_view',globals()['validate_view'])
    transact=hooks.get('transact',globals()['transact'])
    validate_view(source)
    if event.reactant_view_id!=source.view_id or event.product_state_id==source.state_id:raise FormalError('Invalid event state/view references')
    if event.event_type=='ISC':raise RuleUnsupported('spin.transition','requires_spin_transition: no spin dynamics executor')
    violations=[]
    def require(ok,rule,message):
        if not ok:
            violations.append({'rule_id':rule,'original_status':'unsupported','reason':message})
            if rule not in waivers:raise RuleUnsupported(rule,message)
    if event.reactive_template:
        if waivers:raise Unsupported('Reactive-center validation is not waived by the generic adapter')
        return execute_reactive(source,event,_compiler=compile_reactive,_validator=validate_view,_ledger=validate_ledger,_charge=charge)
    atoms=atoms_by_id(source);actions=event.actions;after=None;report={}
    if event.rule_template in SPECS:
        try:actions,after,report=compile_template(source,event)
        except Unsupported as exc:
            require(False,'template.legacy_scope',str(exc))
    elif event.rule_template and event.rule_template.split('.')[0] in FAMILIES:
        try:actions,report=compile_stage(source,event)
        except Unsupported as exc:
            code='snar.ring_scope' if 'ring' in str(exc) else 'snar.leaving_scope' if 'Leaving group' in str(exc) else 'template.legacy_scope'
            require(False,code,str(exc))
    elif event.rule_template:
        require(False,'template.unregistered','Unregistered rule template: '+event.rule_template)
    elif event.template_id and event.event_type in {'RADICAL','PERICYCLIC'}:
        compiler=compile_radical if event.event_type=='RADICAL' else compile_pericyclic
        try:actions,report=compiler(source,event)
        except Unsupported as exc:
            code='template.unregistered' if 'Unregistered' in str(exc) else 'template.legacy_scope'
            require(False,code,str(exc))
    elif event.element_template:
        # Reuse old finite element recipes, with the new state validator.
        from ..v05.execution import execute_event as execute05
        if set(touched(actions or []))&{a.atom_id for a in source.atoms if a.element in {'Fe','Cu','Sn'}}:
            raise RuleUnsupported('template.unregistered','Use a registered v0.7 element recipe')
        target,report=execute05(source,event,_validator=validate_view)
        return target,report
    else:
        if event.event_type=='SN2':
            c=event.bindings.get('electrophile')
            if c not in atoms or atoms[c].element!='C' or bond_valences(source)[c]!=4 or any(b.order!=1 for b in source.bonds if c in {b.a,b.b}):
                raise FormalError('Ordinary SN2 requires a saturated sp3 carbon')
            if sum(atoms[b.b if b.a==c else b.a].element=='C' for b in source.bonds if c in {b.a,b.b})>2:raise Unsupported('Tertiary carbon outside ordinary SN2 profile')
        pairs=[a for a in actions or [] if a.op=='PAIR'];splits=[a for a in actions or [] if a.op=='SPLIT']
        if event.event_type=='POLAR' and len(actions or [])==2 and len(pairs)==len(splits)==1:
            p,s=pairs[0],splits[0];c=p.b
            # Existing bond order increase is NOT nucleophilic replacement (Swern).
            if c==s.a and order(source,p.a,c)==0 and atoms[c].element=='C' and order(source,c,s.b)==1 and any(b.order>1 for b in source.bonds if c in {b.a,b.b}):
                if event.foundation_template=='SNAR_CONCERTED':
                    require(aromatic_ring(source,c) is not None,'snar.ring_scope','Concerted SNAr needs a supported aromatic ring')
                    require(atoms[s.b].element in {'F','Cl','Br','I'},'snar.leaving_scope','SNAr leaving group outside default scope')
                    if not event.evidence or event.mechanism_assumption!='concerted':raise FormalError('Explicit concerted hypothesis and evidence required')
                else:require(False,'substitution.unsaturated_scope','Unclassified sp2/sp substitution needs a dedicated recipe or human exception')
        if event.event_type in {'RADICAL','PERICYCLIC','SET','COORDINATION'}:
            require(False,'template.unregistered','Specialized event needs a registered template')
        else:report=check_event_rule(source,event)
        if touched(actions or [])&{a.atom_id for a in source.atoms if a.element in METAL_GROUP}:
            require(False,'template.unregistered','Metal changes require an explicit element recipe')
    if actions is None:raise FormalError('No typed electron actions supplied')
    if event.actions is not None:
        key=lambda a:json.dumps(a.model_dump(mode='json'),sort_keys=True)
        if Counter(map(key,event.actions))!=Counter(map(key,actions)):raise FormalError('Actions disagree with registered recipe')
    event.actions=actions
    target=transact(source,actions,event.product_state_id,event.product_view_id,after=after)
    if source.electronic_state.coupling=='high_spin' and sum(a.radical_electrons or 0 for a in source.atoms)!=sum(a.radical_electrons or 0 for a in target.atoms):
        raise Unsupported('requires_spin_transition: cannot silently consume high-spin centers')
    report['stereochemistry']=stereo.apply_event(source,target,event)
    report.update(event_id=event.event_id,event_type=event.event_type,rule_template=event.rule_template,
        atomic_commit=True,formal_status='passed',ledger=validate_view(target),execution_mode='electron_actions',
        compiled_actions=[a.model_dump(mode='json') for a in actions],waived_checks=violations,
        rule_version='umr-rules/0.7-r1',chemical_feasibility='unknown',human_approved=False)
    return target,report
