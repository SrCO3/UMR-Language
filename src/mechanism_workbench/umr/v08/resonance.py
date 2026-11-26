"""Shared finite BFS, preserving independent ligand resonance and stereo constraints."""
import sys
from ..v07.resonance import candidates as old_candidates, protected_region as old_protected
from ..v04.resonance import resonance_scan as scan, action_key
from ..v04.ledger import electronic_signature
from ..chemistry import FormalError, Unsupported, identity, hydrogen_connections, bond_orders
from ..v07.transaction import touched
from ..stereochemistry import check_resonance_stereo
from .elements import NEW_METALS
from .ledger import validate_view, signature, compare_views
from .transaction import transact


def organic_shadow(view):
    shadow=view.model_copy(deep=True);ids={a.atom_id for a in shadow.atoms if a.element in NEW_METALS}
    # Covalent new-metal bonds are not capped for resonance; reserve their neighbors.
    shadow.atoms=[a for a in shadow.atoms if a.atom_id not in ids]
    shadow.bonds=[b for b in shadow.bonds if not {b.a,b.b}&ids]
    shadow.metal_centers=[m for m in shadow.metal_centers if m.atom_id not in ids]
    removed={b.ligand_id for b in shadow.metal_bindings if b.metal_atom_id in ids}
    shadow.metal_bindings=[b for b in shadow.metal_bindings if b.metal_atom_id not in ids]
    shadow.ligands=[l for l in shadow.ligands if l.ligand_id not in removed]
    shadow.element_declarations=[d for d in shadow.element_declarations if d.atom_id not in ids]
    return shadow,ids


def protected_region(view):
    shadow,ids=organic_shadow(view)
    neighbors={i for b in view.bonds if {b.a,b.b}&ids for i in (b.a,b.b)}
    return (old_protected(shadow) if shadow.atoms else set())|ids|neighbors


def candidates(view):
    shadow,_=organic_shadow(view);protected=protected_region(view)
    if not shadow.atoms:return
    for rule,mode,actions in old_candidates(shadow):
        if not touched(actions)&protected:yield rule,mode,actions


def apply(source,actions,target_id):
    return transact(source,actions,source.state_id,target_id,preserve_stereo=True)


def validate_resonance(source,target,transform):
    if (source.state_id,target.state_id)!=(transform.state_id,transform.state_id) or (source.view_id,target.view_id)!=(transform.source_view_id,transform.target_view_id):
        raise FormalError('Resonance state/view references differ')
    if electronic_signature(source)!=electronic_signature(target):raise Unsupported('Electronic-state changes are not resonance')
    if identity(source)!=identity(target) or hydrogen_connections(source)!=hydrogen_connections(target) or set(bond_orders(source))!=set(bond_orders(target)):
        raise FormalError('Resonance must preserve identities, H and sigma skeleton')
    validate_view(source);validate_view(target)
    if not any(rule==transform.rule and mode==transform.resonance_mode and action_key(actions)==action_key(transform.actions) for rule,mode,actions in candidates(source)):
        raise Unsupported('unsupported_for_resonance: no registered compatible rule')
    compare_views(apply(source,transform.actions,target.view_id),target)
    check_resonance_stereo(source,target)
    for bond in source.bonds:
        if bond.stereo.get('status')=='specified':
            new=next((b for b in target.bonds if {b.a,b.b}=={bond.a,bond.b}),None)
            if new is None or new.order!=bond.order or new.stereo!=bond.stereo:
                raise Unsupported('Resonance involving a specified alkene needs a geometry-preserving representation rule')


def resonance_scan(view,**kwargs):
    found,links,a=scan(view,_backend=sys.modules[__name__],**kwargs)
    a.rule_version='umr-resonance/0.8'
    for link in links:link.rule_version=a.rule_version
    return found,links,a
