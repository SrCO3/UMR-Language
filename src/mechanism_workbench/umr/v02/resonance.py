import hashlib
from collections import deque

from ..chemistry import (FormalError, Unsupported, atoms_by_id, bond_orders, compare_views,
                         electron_count, graph_signature, hydrogen_connections, identity)
from ..engine import apply_actions
from ..resonance import RULE_VERSION as POLAR_RULE_VERSION, _candidates, resonance_scan as scan01
from .ast import AnalysisRecord, Fishhook, ResonanceTransform
from .electrons import apply_fishhooks, atom_to_bond, bond_to_atom, bond_to_bond
from .ledger import radical_count, validate_view

RULE_VERSION = "umr-resonance/0.2"


def validate_resonance(source,target,transform):
    if source.state_id!=target.state_id or source.state_id!=transform.state_id:
        raise FormalError("Resonance must preserve state_id")
    if (source.view_id,target.view_id)!=(transform.source_view_id,transform.target_view_id):
        raise FormalError("Resonance view references disagree")
    if identity(source)!=identity(target) or hydrogen_connections(source)!=hydrogen_connections(target):
        raise FormalError("Resonance cannot change atoms, isotopes or hydrogen connections (tautomerization is a reaction)")
    if set(bond_orders(source))!=set(bond_orders(target)):
        raise FormalError("Resonance cannot change the sigma skeleton")
    if radical_count(source)!=radical_count(target):
        raise FormalError("Implemented radical resonance preserves local unpaired-electron count")
    if radical_count(source)>1:
        raise Unsupported("unsupported_for_resonance: multi-radical resonance not implemented")
    fish = any(isinstance(a,Fishhook) for a in transform.actions)
    if fish:
        if transform.resonance_mode!="RADICAL" or radical_count(source)!=1 or not all(isinstance(a,Fishhook) for a in transform.actions):
            raise FormalError("Single-electron resonance requires resonance_mode=RADICAL and one declared radical")
        actual = apply_fishhooks(source,transform.actions,source.state_id,target.view_id,preserve_stereo=True)
    else:
        if transform.resonance_mode!="POLAR":
            raise FormalError("Paired resonance requires POLAR mode")
        for a in transform.actions:
            if a.op in {"SPLIT","SHIFT"} and bond_orders(source).get(tuple(sorted((a.a,a.b))),0)<2:
                raise FormalError("Paired resonance cannot consume a sigma pair")
        actual = apply_actions(source,transform.actions,source.state_id,target.view_id,preserve_stereo=True,_validator=validate_view)
    validate_view(target)
    if electron_count(source)!=electron_count(target) or sum(a.formal_charge for a in source.atoms)!=sum(a.formal_charge for a in target.atoms):
        raise FormalError("Resonance charge/electron conservation failed")
    compare_views(actual,target)
    from ..stereochemistry import check_resonance_stereo
    check_resonance_stereo(source,target)


def conjugated_region(view,seeds=None):
    atoms = atoms_by_id(view)
    neighbors = {i:set() for i in atoms}; pi=set()
    for b in view.bonds:
        neighbors[b.a].add(b.b); neighbors[b.b].add(b.a)
        if b.order>=2: pi.update((b.a,b.b))
    capable = {i for i,a in atoms.items() if a.element not in {"H","Li","Na","K","*"}
               and (i in pi or a.radical_electrons or (a.nonbonding_electrons or 0)>=2 or a.formal_charge==1)}
    if seeds is None:
        return capable
    seen = (set(seeds)|{j for i in seeds for j in neighbors.get(i,set())}) & capable
    todo=list(seen)
    while todo:
        for j in neighbors[todo.pop()] & capable - seen:
            seen.add(j); todo.append(j)
    return seen


def resonance_scan(view,*,state=None,seeds=None,max_views=32,max_expansions=128,analysis_id="scan_initial",strategy="charge_priority",heuristic_pruning=None):
    if max_views<1 or max_expansions<1:
        raise ValueError("Resonance limits must be positive")
    analysis = AnalysisRecord(analysis_id=analysis_id,operation="RESONANCE_SCAN",state_id=view.state_id,
        view_id=view.view_id,rule_version=RULE_VERSION,
        input={"conjugated_region":sorted(conjugated_region(view,seeds)),"max_views":max_views,"max_expansions":max_expansions},
        uncertainties=["Finite graph rules; not exhaustive. Orbital geometry and contributor weights remain unknown."])
    try:
        validate_view(view,state)
        if radical_count(view)>1:
            raise Unsupported("Multi-radical resonance is not implemented")
    except Unsupported as exc:
        analysis.status="unsupported_for_resonance"
        analysis.output={"reason":str(exc),"views_found":1,"exhaustive":False}
        return [view],[],analysis
    if radical_count(view)==0:
        found,links,old = scan01(view,seeds=seeds,max_views=max_views,max_expansions=max_expansions,analysis_id=analysis_id,strategy=strategy,heuristic_pruning=heuristic_pruning)
        analysis.status,analysis.output = old.status,old.output
        analysis.output["polar_rule_version"]=POLAR_RULE_VERSION
        return found,[ResonanceTransform(**r.model_dump(),resonance_mode="POLAR") for r in links],analysis
    region=conjugated_region(view,seeds)
    from ..enumeration import ResonanceFrontier
    todo=ResonanceFrontier(strategy);todo.push(view,0)
    found,links,seen=[view],[],{graph_signature(view)}
    expanded=duplicates=heuristic_pruned=0; truncated=False
    heuristic_pruning=strategy!="bfs" if heuristic_pruning is None else heuristic_pruning
    while todo:
        if expanded>=max_expansions:
            truncated=True; break
        depth,current=todo.pop(); expanded+=1
        candidates=[]
        for r in sorted(region):
            if atoms_by_id(current)[r].radical_electrons!=1: continue
            for a in sorted(region-{r}):
                if bond_orders(current).get(tuple(sorted((r,a))),0)!=1: continue
                for b in sorted(region-{r,a}):
                    if bond_orders(current).get(tuple(sorted((a,b))),0)==2:
                        candidates.append(("radical_allyl_pi","RADICAL",[
                            atom_to_bond(r,r,a,"pi1"),bond_to_bond(a,b,"pi1",0,r,a,"pi1"),bond_to_atom(a,b,"pi1",1,b)]))
        candidates.extend((rule,"POLAR",actions) for rule,actions in _candidates(current,region))
        for rule,mode,actions in candidates:
            try:
                candidate = apply_fishhooks(current,actions,view.state_id,"temporary",preserve_stereo=True) if mode=="RADICAL" else apply_actions(current,actions,view.state_id,"temporary",preserve_stereo=True,_validator=validate_view)
                signature=graph_signature(candidate)
                if heuristic_pruning and sum(abs(a.formal_charge) for a in candidate.atoms)>sum(abs(a.formal_charge) for a in view.atoms)+2:
                    heuristic_pruned+=1; continue
                if signature in seen:
                    duplicates+=1; continue
                if len(found)>=max_views:
                    truncated=True; break
                candidate.view_id=f"{view.state_id}.r{hashlib.sha256(repr(signature).encode()).hexdigest()[:12]}"
                transform=ResonanceTransform(transform_id=f"res.{analysis_id}.{len(links)+1}",state_id=view.state_id,
                    source_view_id=current.view_id,target_view_id=candidate.view_id,actions=actions,rule=rule,resonance_mode=mode)
                validate_resonance(current,candidate,transform)
            except (FormalError,Unsupported):
                continue
            seen.add(signature); found.append(candidate); links.append(transform); todo.push(candidate,depth+1)
        if truncated: break
    analysis.status="truncated" if truncated else "complete" if links else "not_found"
    analysis.output={"views_found":len(found),"expanded_views":expanded,"duplicate_or_cycle_edges_skipped":duplicates,
                     "truncated":truncated,"scope":["radical_allyl_pi","benzyl_Kekule_via_allyl","paired_local_rules"],"exhaustive":False}
    analysis.output.update(strategy=strategy,heuristic_pruning=heuristic_pruning,heuristic_pruned=heuristic_pruned)
    return found,links,analysis
