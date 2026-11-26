from collections import Counter
import json
from ..chemistry import FormalError, Unsupported, atoms_by_id
from ..v07.execution import execute_event as legacy
from .ledger import validate_view
from .transaction import transact
from ..v06.stereo import apply_event
from .pericyclic import PROFILES, compile_pericyclic
from . import net


def execute_event(source,event,*,waivers=frozenset(),_hooks=None):
    hooks=_hooks or {}
    validate_view=hooks.get('validate_view',globals()['validate_view'])
    transact=hooks.get('transact',globals()['transact'])
    validate_view(source)
    if event.reactant_view_id!=source.view_id or event.product_state_id==source.state_id:
        raise FormalError('Invalid event state/view references')
    if event.element_template in {'LA_L_ASSOCIATION','LA_L_DISSOCIATION','CEIV_ONE_ELECTRON_OXIDATION','SMII_ONE_ELECTRON_REDUCTION'}:
        from .lanthanides import execute
        target,report=execute(source,event)
    elif event.execution_level=='net_transformation':
        if waivers:raise Unsupported('Net templates cannot waive their own scope')
        target,report=net.execute(source,event)
    elif event.template_id in PROFILES:
        if event.event_type!='PERICYCLIC':raise FormalError('Family must be declared as a pericyclic event')
        actions,report=compile_pericyclic(source,event)
        key=lambda a:json.dumps(a.model_dump(mode='json'),sort_keys=True)
        if event.actions is not None and Counter(map(key,event.actions))!=Counter(map(key,actions)):
            raise FormalError('Declared actions differ from family compilation')
        event.actions=actions
        target=transact(source,actions,event.product_state_id,event.product_view_id)
        # Preserve unknown constraints as side evidence, not invented atom tags.
        trial=event.model_copy(deep=True);trial.stereo_effects=[]
        report.update(event_id=event.event_id,event_type=event.event_type,formal_status='passed',
                      atomic_commit=True,compiled_actions=[a.model_dump(mode='json') for a in actions],
                      execution_mode='electron_actions',stereochemistry=apply_event(source,target,trial),
                      unresolved_stereo_effects=[s.model_dump(mode='json') for s in event.stereo_effects])
    else:
        from .elements import NEW_METALS
        if any(a.element in NEW_METALS for a in source.atoms):
            # Existing organic polar rules remain usable with remote explicit ions.
            from ..engine import check_event_rule
            from ..chemistry import bond_valences
            from ..v07.transaction import touched
            active=touched(event.actions or [])
            if active & {a.atom_id for a in source.atoms if a.element in NEW_METALS} or event.element_template or event.template_id or event.rule_template:
                raise Unsupported('New-metal reactive profile must be registered explicitly')
            if event.event_type not in {'POLAR','PROTON_TRANSFER','HETEROLYSIS','CARBONYL_ADDITION','SN2','E1_IONIZATION','E1_ELIMINATION','E1CB_DEPROTONATION','E1CB_ELIMINATION','E2'}:
                raise Unsupported('Specialized event with new metal has no compatible adapter')
            pairs=[a for a in event.actions or [] if a.op=='PAIR']; splits=[a for a in event.actions or [] if a.op=='SPLIT']
            if event.event_type=='POLAR' and len(event.actions or [])==2 and len(pairs)==len(splits)==1:
                from ..chemistry import order
                p,s=pairs[0],splits[0]
                if p.b==s.a and order(source,p.a,p.b)==0 and atoms_by_id(source)[p.b].element=='C' and order(source,s.a,s.b)==1 and any(b.order>1 for b in source.bonds if p.b in {b.a,b.b}):
                    raise Unsupported('Unclassified sp2/sp substitution remains restricted even with a spectator metal')
            if event.event_type=='SN2':
                c=event.bindings.get('electrophile');atoms=atoms_by_id(source)
                if c not in atoms or atoms[c].element!='C' or bond_valences(source)[c]!=4 or any(b.order!=1 for b in source.bonds if c in {b.a,b.b}):
                    raise FormalError('Ordinary SN2 requires sp3 carbon')
                if sum(atoms[b.b if b.a==c else b.a].element=='C' for b in source.bonds if c in {b.a,b.b})>2:raise Unsupported('Tertiary ordinary SN2 not registered')
            report=check_event_rule(source,event)
            if not event.actions:raise FormalError('Electron event needs nonempty actions')
            target=transact(source,event.actions,event.product_state_id,event.product_view_id)
            report.update(event_id=event.event_id,event_type=event.event_type,formal_status='passed',execution_mode='electron_actions',
                          atomic_commit=True,compiled_actions=[a.model_dump(mode='json') for a in event.actions],stereochemistry=apply_event(source,target,event))
        else:
            target,report=legacy(source,event,waivers=waivers,_hooks=hooks)
    for claim in event.source_stereo:
        if not set(claim.atom_ids)<=atoms_by_id(target).keys():raise FormalError('Source stereo atom reference absent')
    report.update(rule_version='umr-rules/0.8-r1',source_stereo=[s.model_dump(mode='json') for s in event.source_stereo],
                  source_stereo_applied_to_execution=False,chemical_feasibility='unknown')
    return target,report
