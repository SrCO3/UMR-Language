from __future__ import annotations

import hashlib
from collections import deque
from heapq import heappop, heappush

from .ast import AnalysisRecord, LewisView, Pair, ResonanceTransform, Shift, Split
from .chemistry import (FormalError, Unsupported, atoms_by_id, bond_orders, bond_valences,
                        compare_views, electron_count, graph_signature, hydrogen_connections,
                        identity, validate_view)
from .engine import apply_actions

RULE_VERSION = "umr-resonance/0.1.1"


def validate_resonance(source: LewisView, target: LewisView, transform: ResonanceTransform) -> None:
    if source.state_id != target.state_id or source.state_id != transform.state_id:
        raise FormalError("Resonance views must share one ChemicalState/state_id")
    if (transform.source_view_id, transform.target_view_id) != (source.view_id, target.view_id):
        raise FormalError("Resonance transform view references are inconsistent")
    if identity(source) != identity(target):
        raise FormalError("Resonance cannot change atom identity, isotopes or attached hydrogens")
    if hydrogen_connections(source) != hydrogen_connections(target):
        raise FormalError("Hydrogen movement/tautomerization is a reaction, not resonance")
    if set(bond_orders(source)) != set(bond_orders(target)):
        raise FormalError("Resonance cannot change the sigma connectivity skeleton")
    target_bonds = bond_orders(target)
    for bond in source.bonds:
        if bond.stereo.get("status") == "specified" and target_bonds[tuple(sorted((bond.a, bond.b)))] != bond.order:
            raise Unsupported("Resonance across specified double-bond stereochemistry requires a future stereochemical adapter")
    validate_view(source)
    validate_view(target)
    if (electron_count(source) != electron_count(target) or
            sum(a.formal_charge for a in source.atoms) != sum(a.formal_charge for a in target.atoms)):
        raise FormalError("Resonance must preserve total charge and electron state")
    for action in transform.actions:
        if isinstance(action, (Split, Shift)) and bond_orders(source).get(tuple(sorted((action.a, action.b))), 0) < 2:
            raise FormalError("Resonance may redistribute pi pairs, not consume a sigma pair")
    actual = apply_actions(source, transform.actions, source.state_id, target.view_id, preserve_stereo=True)
    compare_views(actual, target)
    from .stereochemistry import check_resonance_stereo
    check_resonance_stereo(source,target)


def conjugated_region(view: LewisView, seeds: set[int] | None) -> list[int]:
    atoms, valences = atoms_by_id(view), bond_valences(view)
    adjacent = {i: set() for i in atoms}
    pi = set()
    for bond in view.bonds:
        adjacent[bond.a].add(bond.b)
        adjacent[bond.b].add(bond.a)
        if bond.order >= 2:
            pi.update((bond.a, bond.b))
    capable = {i for i, a in atoms.items() if a.element not in {"H", "Li", "Na", "K", "*"}
               and (i in pi or (a.nonbonding_electrons or 0) >= 2 or
                    (a.element in {"C", "B"} and valences[i] == 3 and a.nonbonding_electrons == 0))}
    if seeds is None:
        return sorted(capable)
    initial = set(seeds) | {j for i in seeds for j in adjacent.get(i, set())}
    pending = deque(initial & capable)
    result = set(pending)
    while pending:
        for neighbor in adjacent[pending.popleft()] & capable:
            if neighbor not in result:
                result.add(neighbor)
                pending.append(neighbor)
    return sorted(result)


def _candidates(view: LewisView, region: set[int]):
    atoms, valences = atoms_by_id(view), bond_valences(view)
    adjacency = {i: {} for i in atoms}
    for bond in view.bonds:
        adjacency[bond.a][bond.b] = bond.order
        adjacency[bond.b][bond.a] = bond.order
    def empty_p(i):
        return atoms[i].element in {"C", "B"} and valences[i] == 3 and atoms[i].nonbonding_electrons == 0
    for a in sorted(region):
        for b, ab in sorted(adjacency[a].items()):
            if b not in region:
                continue
            if ab == 1 and atoms[a].nonbonding_electrons >= 2:
                if empty_p(b):
                    yield "lone_pair_empty_p", [Pair(a=a, b=b)]
                for c, bc in sorted(adjacency[b].items()):
                    if c in region and c != a and bc >= 2:
                        yield "lone_pair_pi", [Pair(a=a, b=b), Split(a=b, b=c, target=c)]
            if ab >= 2:
                if atoms[b].element in {"N", "O", "S"} and atoms[a].element == "C":
                    yield "polarized_pi", [Split(a=a, b=b, target=b)]
                for c, bc in sorted(adjacency[b].items()):
                    if c not in region or c == a or bc != 1:
                        continue
                    if empty_p(c):
                        yield "pi_empty_p", [Shift(a=a, b=b, c=c)]
                    for d, cd in sorted(adjacency[c].items()):
                        if d in region and d not in {a, b} and cd >= 2:
                            yield "pi_pi", [Shift(a=a, b=b, c=c), Split(a=c, b=d, target=d)]


def resonance_scan(view: LewisView, *, seeds: set[int] | None = None, max_views: int = 32,
                   max_expansions: int = 128, analysis_id: str = "scan_initial",
                   strategy: str = "charge_priority", heuristic_pruning: bool | None = None) -> tuple[list[LewisView], list[ResonanceTransform], AnalysisRecord]:
    if max_views < 1 or max_expansions < 1:
        raise ValueError("Resonance enumeration bounds must be positive")
    region = conjugated_region(view, seeds)
    analysis = AnalysisRecord(analysis_id=analysis_id, operation="RESONANCE_SCAN", state_id=view.state_id,
                              view_id=view.view_id, rule_version=RULE_VERSION,
                              input={"reaction_center": sorted(seeds) if seeds is not None else [],
                                     "conjugated_region": region, "max_views": max_views, "max_expansions": max_expansions},
                              uncertainties=["Finite local Lewis rules, not an exhaustive quantum resonance expansion.",
                                             "Orbital alignment and contributor weights are not established from 2D connectivity."])
    try:
        validate_view(view)
    except Unsupported as exc:
        analysis.status, analysis.output = "unsupported", {"reason": str(exc), "views_found": 1}
        return [view], [], analysis
    # Explore less charge-separated views first so limited scans reach ordinary
    # cation contributors before expanding the combinatorial charge-separated set.
    # This is an enumeration heuristic, not an energy or contributor-weight model.
    from .enumeration import ResonanceFrontier
    pending = ResonanceFrontier(strategy)
    pending.push(view, 0)
    heuristic_pruning = strategy != "bfs" if heuristic_pruning is None else heuristic_pruning
    heuristic_pruned = 0
    found = [view]
    visited = {graph_signature(view)}
    transforms, duplicates, expanded, truncated = [], 0, 0, False
    charge_limit = sum(abs(a.formal_charge) for a in view.atoms) + 2
    while pending:
        if expanded >= max_expansions:
            truncated = True
            break
        depth, current = pending.pop()
        expanded += 1
        for rule, actions in _candidates(current, set(region)):
            try:
                candidate = apply_actions(current, actions, view.state_id, "temporary", preserve_stereo=True)
                if heuristic_pruning and sum(abs(a.formal_charge) for a in candidate.atoms) > charge_limit:
                    heuristic_pruned += 1
                    continue
                signature = graph_signature(candidate)
                if signature in visited:
                    duplicates += 1
                    continue
                if len(found) >= max_views:
                    truncated = True
                    break
                short_hash = hashlib.sha256(repr(signature).encode()).hexdigest()[:12]
                candidate.view_id = f"{view.state_id}.r{short_hash}"
                transform = ResonanceTransform(
                    transform_id=f"res.{analysis_id}.{len(transforms)+1}", state_id=view.state_id,
                    source_view_id=current.view_id, target_view_id=candidate.view_id,
                    actions=actions, rule=rule, rule_version=RULE_VERSION,
                )
                validate_resonance(current, candidate, transform)
            except (FormalError, Unsupported):
                continue
            visited.add(signature)
            found.append(candidate)
            transforms.append(transform)
            pending.push(candidate, depth + 1)
        if truncated:
            break
    analysis.status = "truncated" if truncated else "complete" if transforms else "not_found"
    analysis.output = {"views_found": len(found), "duplicate_or_cycle_edges_skipped": duplicates,
                       "expanded_views": expanded, "truncated": truncated,
                       "scope": ["lone_pair_pi", "pi_empty_p", "lone_pair_empty_p", "polarized_pi", "pi_pi"],
                       "charge_separation_budget": charge_limit, "exhaustive": False}
    analysis.output["enumeration_priority"] = "lower_absolute_formal_charge_first"
    analysis.output.update(strategy=strategy, heuristic_pruning=heuristic_pruning,
                           heuristic_pruned=heuristic_pruned,
                           enumeration_priority="depth_then_absolute_charge" if strategy == "bfs" else "lower_absolute_formal_charge_first")
    return found, transforms, analysis
