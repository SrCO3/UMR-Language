"""Existing finite local rules, with protected special electronic-structure regions."""
import hashlib
import json
from collections import Counter

from ..chemistry import FormalError, Unsupported, atoms_by_id, bond_orders, hydrogen_connections, identity, order
from ..engine import apply_actions
from ..enumeration import ResonanceFrontier
from ..resonance import _candidates
from ..v02.ast import Fishhook
from ..v02.electrons import apply_fishhooks, atom_to_bond, bond_to_atom, bond_to_bond
from ..v02.resonance import conjugated_region
from .ast import AnalysisRecord, ResonanceTransform
from .ledger import compare_views, electronic_signature, signature, validate_view


def protected_region(view):
    ids={c.atom_id for c in view.reactive_centers}|{i for r in view.aryne_regions for i in r.ring_atom_ids}
    radicals={a.atom_id for a in view.atoms if a.radical_electrons}
    if sum(a.radical_electrons or 0 for a in view.atoms)>1: ids |= radicals
    return conjugated_region(view,ids) if ids else set()


def candidates(view):
    region=conjugated_region(view)-protected_region(view)
    for rule,actions in _candidates(view,region): yield rule,"POLAR",actions
    atoms=atoms_by_id(view)
    for r in sorted(region):
        if atoms[r].radical_electrons!=1: continue
        for a in sorted(region-{r}):
            if order(view,r,a)!=1: continue
            for b in sorted(region-{r,a}):
                if order(view,a,b)==2:
                    yield "radical_allyl_pi","RADICAL",[atom_to_bond(r,r,a,"pi1"),bond_to_bond(a,b,"pi1",0,r,a,"pi1"),bond_to_atom(a,b,"pi1",1,b)]


def action_key(actions):
    return Counter(json.dumps(a.model_dump(mode="json"),sort_keys=True) for a in actions)


def apply(source,actions,target_id):
    if all(isinstance(a,Fishhook) for a in actions):
        return apply_fishhooks(source,actions,source.state_id,target_id,preserve_stereo=True,_validator=validate_view)
    if any(isinstance(a,Fishhook) for a in actions): raise Unsupported("Mixed resonance mode is not registered")
    return apply_actions(source,actions,source.state_id,target_id,preserve_stereo=True,_validator=validate_view)


def validate_resonance(source,target,transform):
    if source.state_id!=target.state_id or source.state_id!=transform.state_id:
        raise FormalError("Resonance must preserve state_id; spin transitions are distinct states")
    if (source.view_id,target.view_id)!=(transform.source_view_id,transform.target_view_id):
        raise FormalError("Resonance references disagree")
    if electronic_signature(source)!=electronic_signature(target):
        raise Unsupported("unsupported_for_resonance: electronic-state/profile/aryne-representation change is not resonance")
    if source.metal_centers or target.metal_centers:
        validate_view(source);validate_view(target)
        from ..v03.resonance import validate_resonance as validate03
        return validate03(source,target,transform)
    if identity(source)!=identity(target) or hydrogen_connections(source)!=hydrogen_connections(target) or set(bond_orders(source))!=set(bond_orders(target)):
        raise FormalError("Resonance must preserve atoms, hydrogens and sigma skeleton")
    validate_view(source);validate_view(target)
    key=action_key(transform.actions)
    if not any(rule==transform.rule and mode==transform.resonance_mode and action_key(actions)==key for rule,mode,actions in candidates(source)):
        raise Unsupported("unsupported_for_resonance: no registered electronic-state-compatible transformation")
    compare_views(apply(source,transform.actions,target.view_id),target)
    from ..stereochemistry import check_resonance_stereo
    check_resonance_stereo(source,target)


def resonance_scan(view,*,state=None,seeds=None,max_views=32,max_expansions=128,analysis_id="scan_initial",strategy="bfs",heuristic_pruning=False,_backend=None):
    validate_view=globals()['validate_view'] if _backend is None else _backend.validate_view
    protected_region=globals()['protected_region'] if _backend is None else _backend.protected_region
    candidates=globals()['candidates'] if _backend is None else _backend.candidates
    apply=globals()['apply'] if _backend is None else _backend.apply
    signature=globals()['signature'] if _backend is None else _backend.signature
    validate_resonance=globals()['validate_resonance'] if _backend is None else _backend.validate_resonance
    if max_views<1 or max_expansions<1: raise ValueError("Positive resonance enumeration budgets required")
    analysis=AnalysisRecord(analysis_id=analysis_id,operation="RESONANCE_SCAN",state_id=view.state_id,view_id=view.view_id,
        rule_version="umr-resonance/0.4",input={"max_views":max_views,"max_expansions":max_expansions,"reaction_center":sorted(seeds or [])},
        uncertainties=["Finite registered Lewis rules, no contributor weights or wavefunction equivalence.",
                       "Special-center conjugation, multi-radical coupling and aryne triple/diradical interchange are unimplemented."])
    try: validate_view(view,state)
    except Unsupported as exc:
        analysis.status="unsupported_for_resonance";analysis.output={"reason":str(exc),"views_found":1,"truncated":False,"exhaustive":False}
        return [view],[],analysis
    if view.metal_centers and _backend is None:
        from ..v03.resonance import resonance_scan as scan03
        return scan03(view,state=state,seeds=seeds,max_views=max_views,max_expansions=max_expansions,analysis_id=analysis_id,strategy=strategy,heuristic_pruning=heuristic_pruning)
    protected=protected_region(view)
    todo=ResonanceFrontier(strategy);todo.push(view,0)
    found,links,seen=[view],[],{signature(view)}
    expanded=duplicates=pruned=0;truncated=False
    while todo:
        if expanded>=max_expansions: truncated=True;break
        depth,current=todo.pop();expanded+=1
        for rule,mode,actions in candidates(current):
            try:
                target=apply(current,actions,"temporary");sig=signature(target)
                if sig in seen: duplicates+=1;continue
                if heuristic_pruning and sum(abs(a.formal_charge) for a in target.atoms)>sum(abs(a.formal_charge) for a in view.atoms)+2:
                    pruned+=1;continue
                if len(found)>=max_views: truncated=True;break
                target.view_id=f"{view.state_id}.r{hashlib.sha256(repr(sig).encode()).hexdigest()[:12]}"
                link=ResonanceTransform(transform_id=f"res.{analysis_id}.{len(links)+1}",state_id=view.state_id,
                    source_view_id=current.view_id,target_view_id=target.view_id,actions=actions,rule=rule,resonance_mode=mode,rule_version="umr-resonance/0.4")
                validate_resonance(current,target,link)
            except (FormalError,Unsupported): continue
            seen.add(sig);found.append(target);links.append(link);todo.push(target,depth+1)
        if truncated: break
    analysis.status="truncated" if truncated else "unsupported_for_resonance" if protected else "complete" if links else "not_found"
    analysis.output={"views_found":len(found),"expanded_views":expanded,"duplicate_or_cycle_edges_skipped":duplicates,
        "truncated":truncated,"exhaustive":False,"strategy":strategy,"heuristic_pruned":pruned,
        "unsupported_region_atom_ids":sorted(protected),"supported_scope":"registered polar / single-radical rules outside protected conjugated regions",
        "reaction_execution_status":"independent"}
    return found,links,analysis
