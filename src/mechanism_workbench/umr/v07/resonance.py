import sys
from ..chemistry import FormalError,Unsupported,identity,hydrogen_connections,bond_orders,order
from ..ast import Shift
from ..v04.resonance import resonance_scan as scan04,action_key
from ..v04.ledger import electronic_signature
from ..v06.resonance import candidates as candidates06
from ..v05.ledger import signature,compare_views
from .ledger import validate_view,project_symbols,symbols
from .elements import METAL_GROUP,MAIN
from .transaction import transact,touched
from .substitution import aromatic_ring


def protected_region(view):
    from ..v04.resonance import protected_region as special
    return special(project_symbols(view))|{a.atom_id for a in view.atoms if a.element in MAIN|set(METAL_GROUP)}|symbols(view)


def candidates(view):
    shadow=project_symbols(view)
    for a in shadow.atoms:
        if a.radical_electrons is None:a.radical_electrons=0
    protected=protected_region(view)
    for rule,mode,actions in candidates06(shadow):
        if not touched(actions)&protected:yield rule,mode,actions
    seen=set()
    for a in view.atoms:
        if a.element!='C':continue
        ring=aromatic_ring(view,a.atom_id)
        if not ring or tuple(sorted(ring)) in seen or set(ring)&protected:continue
        seen.add(tuple(sorted(ring)))
        yield 'kekule_six_member_CN','POLAR',[Shift(a=ring[i],b=ring[(i+1)%6],c=ring[(i+2)%6]) for i in (0,2,4)]


def apply(source,actions,target_id):
    return transact(source,actions,source.state_id,target_id,preserve_stereo=True)


def validate_resonance(source,target,transform):
    if source.state_id!=target.state_id or source.state_id!=transform.state_id:raise FormalError('Resonance cannot change chemical state')
    if (source.view_id,target.view_id)!=(transform.source_view_id,transform.target_view_id):raise FormalError('Resonance view references disagree')
    if electronic_signature(source)!=electronic_signature(target):raise Unsupported('Electronic-state change is not resonance')
    if identity(source)!=identity(target) or hydrogen_connections(source)!=hydrogen_connections(target) or set(bond_orders(source))!=set(bond_orders(target)):
        raise FormalError('Resonance must preserve atom identities, hydrogen and sigma skeleton')
    if not any(rule==transform.rule and mode==transform.resonance_mode and action_key(actions)==action_key(transform.actions) for rule,mode,actions in candidates(source)):
        raise Unsupported('unsupported_for_resonance: unregistered transformation')
    compare_views(apply(source,transform.actions,target.view_id),target)
    from ..stereochemistry import check_resonance_stereo
    check_resonance_stereo(source,target)


def resonance_scan(view,**kwargs):
    found,links,a=scan04(view,_backend=sys.modules[__name__],**kwargs)
    a.rule_version='umr-resonance/0.7'
    for link in links:link.rule_version=a.rule_version
    return found,links,a
