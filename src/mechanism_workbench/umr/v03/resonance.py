"""Organic ligand resonance, holding metal bindings/geometry/electronic state fixed."""
import hashlib
from heapq import heappop, heappush

from ..chemistry import (FormalError, Unsupported, atoms_by_id, bond_orders, electron_count,
                         hydrogen_connections, identity)
from ..resonance import _candidates
from ..v02.ast import Fishhook, ResonanceTransform
from ..v02.electrons import atom_to_bond, bond_to_atom, bond_to_bond
from ..v02.resonance import conjugated_region
from .ast import AnalysisRecord
from .execution import apply_organic
from .ledger import compare_views, metal_signature, signature, validate_view

RULE_VERSION="umr-metal-ligand-resonance/0.3.1"


def validate_resonance(source,target,transform):
    if (source.state_id,target.state_id)!=(transform.state_id,transform.state_id): raise FormalError("Resonance must preserve state_id")
    if (source.view_id,target.view_id)!=(transform.source_view_id,transform.target_view_id): raise FormalError("Resonance view references disagree")
    if metal_signature(source)!=metal_signature(target):
        raise FormalError("Coordination, ligand identity, oxidation state, spin or geometry change is not organic resonance; representation conversions need a separate adapter")
    if identity(source)!=identity(target) or hydrogen_connections(source)!=hydrogen_connections(target):
        raise FormalError("Resonance cannot transfer H or change atoms/isotopes")
    if set(bond_orders(source))!=set(bond_orders(target)): raise FormalError("Resonance cannot change sigma skeleton")
    radicals=lambda v:sum(a.radical_electrons for a in v.atoms)
    if radicals(source)!=radicals(target) or radicals(source)>1: raise Unsupported("Only zero/one organic-radical resonance is supported")
    fish=any(isinstance(a,Fishhook) for a in transform.actions)
    if transform.resonance_mode!=("RADICAL" if fish else "POLAR"): raise FormalError("Resonance mode disagrees with electron actions")
    for a in transform.actions:
        if a.op in {"SPLIT","SHIFT"} and bond_orders(source).get(tuple(sorted((a.a,a.b))),0)<2:
            raise FormalError("Organic resonance cannot consume a sigma pair")
    validate_view(target)
    actual=apply_organic(source,transform.actions,source.state_id,target.view_id,preserve_stereo=True)
    compare_views(actual,target)
    if electron_count(source)!=electron_count(target) or sum(a.formal_charge for a in source.atoms)!=sum(a.formal_charge for a in target.atoms):
        raise FormalError("Resonance must preserve electron state and total charge")
    from ..stereochemistry import check_resonance_stereo
    check_resonance_stereo(source,target)
    orders=bond_orders(target)
    if any(b.stereo.get("status")=="specified" and orders[tuple(sorted((b.a,b.b)))]!=b.order for b in source.bonds):
        raise Unsupported("Resonance across specified double-bond stereo needs a stereochemical adapter")


def resonance_scan(view,*,state=None,seeds=None,max_views=32,max_expansions=128,analysis_id="scan_initial",strategy="charge_priority",heuristic_pruning=None):
    if max_views<1 or max_expansions<1: raise ValueError("Positive resonance enumeration bounds required")
    mids={m.atom_id for m in view.metal_centers}
    region=set(conjugated_region(view,seeds))-mids
    analysis=AnalysisRecord(analysis_id=analysis_id,operation="RESONANCE_SCAN",state_id=view.state_id,
        view_id=view.view_id,rule_version=RULE_VERSION,input={"organic_region":sorted(region),
            "excluded_metal_atom_ids":sorted(mids),"max_views":max_views,"max_expansions":max_expansions},
        uncertainties=["Finite local organic rules; not exhaustive. Metal-ligand delocalization is not enumerated.",
                       "Orbital alignment, contributor weights and ambiguous oxidation-state assignments remain unknown."])
    try:
        validate_view(view,state)
        if sum(a.radical_electrons for a in view.atoms)>1: raise Unsupported("Multi-radical resonance not implemented")
    except Unsupported as exc:
        analysis.status="unsupported_for_resonance";analysis.output={"reason":str(exc),"views_found":1,"exhaustive":False}
        return [view],[],analysis
    found,links,seen=[view],[],{signature(view)}
    # Search priority only; no physical contributor weights are inferred.
    from ..enumeration import ResonanceFrontier
    todo=ResonanceFrontier(strategy);todo.push(view,0)
    heuristic_pruning=strategy!="bfs" if heuristic_pruning is None else heuristic_pruning
    expanded=duplicates=heuristic_pruned=0;truncated=False
    while todo:
        if expanded>=max_expansions: truncated=True;break
        depth,current=todo.pop();expanded+=1
        candidates=[(rule,"POLAR",a) for rule,a in _candidates(current,region)]
        atoms=atoms_by_id(current);orders=bond_orders(current)
        for r in sorted(region):
            if atoms[r].radical_electrons!=1: continue
            for a in sorted(region-{r}):
                if orders.get(tuple(sorted((r,a))))!=1: continue
                for b in sorted(region-{r,a}):
                    if orders.get(tuple(sorted((a,b))))==2:
                        candidates.append(("radical_allyl_pi","RADICAL",[atom_to_bond(r,r,a,"pi1"),
                            bond_to_bond(a,b,"pi1",0,r,a,"pi1"),bond_to_atom(a,b,"pi1",1,b)]))
        for rule,mode,actions in candidates:
            try:
                candidate=apply_organic(current,actions,view.state_id,"temporary",preserve_stereo=True)
                key=signature(candidate)
                if heuristic_pruning and sum(abs(a.formal_charge) for a in candidate.atoms)>sum(abs(a.formal_charge) for a in view.atoms)+2:
                    heuristic_pruned+=1;continue
                if key in seen: duplicates+=1;continue
                if len(found)>=max_views: truncated=True;break
                candidate.view_id=f"{view.state_id}.r{hashlib.sha256(repr(key).encode()).hexdigest()[:12]}"
                link=ResonanceTransform(transform_id=f"res.{analysis_id}.{len(links)+1}",state_id=view.state_id,
                    source_view_id=current.view_id,target_view_id=candidate.view_id,actions=actions,rule=rule,
                    resonance_mode=mode,rule_version=RULE_VERSION)
                validate_resonance(current,candidate,link)
            except (FormalError,Unsupported): continue
            seen.add(key);found.append(candidate);links.append(link)
            todo.push(candidate,depth+1)
        if truncated: break
    analysis.status="truncated" if truncated else "complete" if links else "not_found"
    analysis.output={"views_found":len(found),"expanded_views":expanded,"duplicate_or_cycle_edges_skipped":duplicates,
        "truncated":truncated,"scope":["organic_paired_local_rules","organic_allyl_radical"],
        "metal_ligand_resonance":"unsupported","exhaustive":False,
        "enumeration_priority":"lower_absolute_formal_charge_first"}
    analysis.output.update(strategy=strategy,heuristic_pruning=heuristic_pruning,heuristic_pruned=heuristic_pruned,
        enumeration_priority="depth_then_absolute_charge" if strategy=="bfs" else "lower_absolute_formal_charge_first")
    return found,links,analysis
