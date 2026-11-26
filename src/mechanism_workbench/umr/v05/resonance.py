"""Reuse the finite BFS scanner; exclude metal/expanded-valence transformations."""
import sys
from collections import Counter
from ..chemistry import FormalError, Unsupported, identity, hydrogen_connections, bond_orders
from ..engine import apply_actions
from ..v02.ast import Fishhook
from ..v02.electrons import apply_fishhooks
from ..v04.resonance import (candidates as candidates04,protected_region as protected04,
                             action_key,resonance_scan as scan04)
from ..v04.ledger import electronic_signature
from .elements import MAIN, METAL_GROUP, charge
from .ledger import signature, compare_views, validate_view
from .coordination import checks

def protected_region(view):
    return protected04(view)|{a.atom_id for a in view.atoms if a.element in MAIN|set(METAL_GROUP)}

def candidates(view):
    protected=protected_region(view)
    for rule,mode,actions in candidates04(view):
        touched=set()
        for a in actions:
            if a.op=='FISHHOOK':
                for endpoint in (a.source,a.destination):
                    touched.update([endpoint.atom_id] if hasattr(endpoint,'atom_id') else [endpoint.a,endpoint.b])
            else:touched.update(getattr(a,k) for k in ('a','b','c') if hasattr(a,k))
        if not touched&protected:yield rule,mode,actions

def apply(source,actions,target_id,*,_validator=validate_view):
    reserved,_=checks(source)
    spent=Counter(a.a for a in actions if a.op=='PAIR')
    for atom_id,count in spent.items():
        atom=next(x for x in source.atoms if x.atom_id==atom_id)
        if atom.nonbonding_electrons-(atom.radical_electrons or 0)-reserved[atom_id]<2*count:
            raise FormalError('Reserved coordination pair cannot be spent by resonance')
    executor=apply_fishhooks if actions and all(isinstance(a,Fishhook) for a in actions) else apply_actions
    return executor(source,actions,source.state_id,target_id,preserve_stereo=True,_validator=_validator,_charge=charge)

def validate_resonance(source,target,transform,*,_validator=validate_view,_apply=apply):
    validate_view=_validator
    if (source.state_id,target.state_id)!=(transform.state_id,transform.state_id) or (source.view_id,target.view_id)!=(transform.source_view_id,transform.target_view_id):
        raise FormalError('Resonance references or state_id differ')
    if electronic_signature(source)!=electronic_signature(target):raise Unsupported('unsupported_for_resonance: electronic state differs')
    if identity(source)!=identity(target) or hydrogen_connections(source)!=hydrogen_connections(target) or set(bond_orders(source))!=set(bond_orders(target)):
        raise FormalError('Resonance must preserve atoms, hydrogen and sigma skeleton')
    validate_view(source);validate_view(target)
    if not any(rule==transform.rule and mode==transform.resonance_mode and action_key(actions)==action_key(transform.actions) for rule,mode,actions in candidates(source)):
        raise Unsupported('unsupported_for_resonance: element/metal conversion is outside registered organic rules')
    compare_views(_apply(source,transform.actions,target.view_id),target)
    from ..stereochemistry import check_resonance_stereo
    check_resonance_stereo(source,target)

def resonance_scan(view,**kwargs):
    views,links,a=scan04(view,_backend=sys.modules[__name__],**kwargs)
    a.rule_version='umr-resonance/0.5'
    a.uncertainties.append('Metal/expanded-valence representation rules excluded; organic ligand fragments remain eligible.')
    a.output['supported_scope']='v0.4 organic polar/single-radical rules outside protected atom centers; partial coverage'
    for link in links:link.rule_version='umr-resonance/0.5'
    return views,links,a
