"""Registered alternating six-carbon Kekule transformation, plus existing scans."""
import sys
from ..ast import Shift
from ..chemistry import FormalError, Unsupported, order, atoms_by_id, identity, hydrogen_connections, bond_orders
from ..v04.ledger import check_isolated_ring, electronic_signature
from ..v04.resonance import resonance_scan as shared_scan, action_key
from ..v05.resonance import candidates as old_candidates, apply as apply05, protected_region, validate_resonance as old_validate
from .ledger import validate_view, signature, compare_views

def apply(source,actions,target_id):return apply05(source,actions,target_id,_validator=validate_view)


def kekule_actions(view,ring):
    atoms=atoms_by_id(view)
    if len(ring)!=6 or len(set(ring))!=6 or any(i not in atoms or atoms[i].element!='C' for i in ring):
        raise Unsupported('Kekule rule requires six distinct carbon atoms')
    check_isolated_ring(view,ring)
    if any(atoms[i].radical_electrons or any(i in r.ring_atom_ids for r in view.aryne_regions) for i in ring):
        raise Unsupported('Registered Kekule scope excludes aryne and radical centers')
    orders=[order(view,i,ring[(j+1)%6]) for j,i in enumerate(ring)]
    if orders not in ([1,2,1,2,1,2],[2,1,2,1,2,1]):raise FormalError('Not an alternating six-carbon ring')
    # Three pi-pairs relocate as a single representation transaction.
    return [Shift(a=ring[j],b=ring[(j+1)%6],c=ring[(j+2)%6]) for j in range(6) if orders[j]==2]


def rings(view):
    atoms=atoms_by_id(view);adj={i:[] for i,a in atoms.items() if a.element=='C'}
    for b in view.bonds:
        if b.a in adj and b.b in adj:
            adj[b.a].append(b.b);adj[b.b].append(b.a)
    seen=set()
    for start in sorted(adj):
        todo=[[start]]
        while todo:
            path=todo.pop()
            if len(path)==6:
                if start in adj[path[-1]]:
                    key=tuple(sorted(path))
                    if key not in seen:seen.add(key);yield path
                continue
            todo.extend(path+[j] for j in adj[path[-1]] if j>start and j not in path)


def candidates(view):
    yield from old_candidates(view)
    for ring in rings(view):
        try:actions=kekule_actions(view,ring)
        except (FormalError,Unsupported):continue
        yield 'kekule_six_carbon','POLAR',actions
        yield 'kekule_six_carbon','POLAR',kekule_actions(view,[ring[0]]+list(reversed(ring[1:])))


def validate_resonance(source,target,transform):
    if transform.rule!='kekule_six_carbon':return old_validate(source,target,transform,_validator=validate_view,_apply=apply)
    if (source.state_id,target.state_id)!=(transform.state_id,transform.state_id) or (source.view_id,target.view_id)!=(transform.source_view_id,transform.target_view_id):
        raise FormalError('Kekule view/state references disagree')
    if electronic_signature(source)!=electronic_signature(target):raise Unsupported('Kekule transformation cannot change electronic state')
    if identity(source)!=identity(target) or hydrogen_connections(source)!=hydrogen_connections(target) or set(bond_orders(source))!=set(bond_orders(target)):
        raise FormalError('Kekule transformation must preserve atom/H/sigma skeleton')
    if not any(rule==transform.rule and mode==transform.resonance_mode and action_key(actions)==action_key(transform.actions) for rule,mode,actions in candidates(source)):
        raise Unsupported('No matching registered Kekule transform')
    compare_views(apply(source,transform.actions,target.view_id),target)
    from ..stereochemistry import check_resonance_stereo
    check_resonance_stereo(source,target)


def resonance_scan(view,**kwargs):
    found,links,analysis=shared_scan(view,_backend=sys.modules[__name__],**kwargs)
    analysis.rule_version='umr-resonance/0.6'
    analysis.uncertainties.append('Registered isolated six-carbon Kekule views included; no aryne triple/diradical equivalence inferred.')
    for link in links:link.rule_version='umr-resonance/0.6'
    return found,links,analysis
