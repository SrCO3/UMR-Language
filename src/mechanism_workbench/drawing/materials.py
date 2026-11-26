"""Projection of a frozen event onto its depicted components and full material pool."""
from collections import Counter

from ..umr.stereochemistry import stereo_assertion


def raw(obj):
    return obj.model_dump(mode="json") if hasattr(obj, "model_dump") else obj


def contacts(view):
    for b in getattr(view, "metal_bindings", []):
        yield (b.mode, b.metal_atom_id, tuple(sorted(b.donor_atom_ids)), b.hapticity)


def components(view):
    neighbors = {a.atom_id: set() for a in view.atoms}
    edges = [(b.a, b.b) for b in view.bonds]
    edges.extend((m, i) for _, m, donors, _ in contacts(view) for i in donors)
    for a, b in edges:
        neighbors[a].add(b)
        neighbors[b].add(a)
    rows = []
    while neighbors:
        todo = [next(iter(neighbors))]
        group = set()
        while todo:
            i = todo.pop()
            if i in neighbors:
                group.add(i)
                todo.extend(neighbors.pop(i) - group)
        rows.append(group)
    return sorted(rows, key=min)


def action_sites(actions):
    touched, pairs = set(), set()
    for action in map(raw, actions):
        touched.update(action[k] for k in ("a","b","c","target") if type(action.get(k)) is int)
        if action.get("op") == "PAIR":
            pairs.add((action["a"], action["b"]))
        elif action.get("op") == "SHIFT":
            pairs.add((action["b"], action["c"]))
        elif action.get("op") == "FISHHOOK":
            ends = [[end[k] for k in ("atom_id","a","b") if k in end]
                    for end in (action["source"], action["destination"])]
            for ids in ends:
                touched.update(ids)
            pairs.update((a,b) for a in ends[0] for b in ends[1] if a != b)
    return touched, pairs


def differences(source, target):
    ids, pairs = set(), set()
    if target is None:
        return ids, pairs
    def atoms(v):
        return {a.atom_id: (a.element,a.isotope,a.hydrogens,a.formal_charge,a.nonbonding_electrons,
                            a.radical_electrons,stereo_assertion(a.stereo)) for a in v.atoms}
    def bonds(v):
        return {tuple(sorted((b.a,b.b))): (b.order,stereo_assertion(b.stereo)) for b in v.bonds}
    old, new = atoms(source), atoms(target)
    ids.update(i for i in old.keys() | new.keys() if old.get(i) != new.get(i))
    old, new = bonds(source), bonds(target)
    for edge in old.keys() | new.keys():
        if old.get(edge) != new.get(edge):
            ids.update(edge)
            if edge not in old:
                pairs.add(edge)
    for _, metal, donors, _ in set(contacts(source)) ^ set(contacts(target)):
        ids.update((metal,*donors))
        pairs.update((metal,i) for i in donors)
    old = {m.atom_id: raw(m) for m in getattr(source,"metal_centers",[])}
    new = {m.atom_id: raw(m) for m in getattr(target,"metal_centers",[])}
    ids.update(i for i in old.keys() | new.keys() if old.get(i) != new.get(i))
    return ids, pairs


def phase_inventory(document, event_id):
    phases = getattr(document, "phases", [])
    current = next((n for n,p in enumerate(phases) if event_id in p.event_ids), None)
    return {i for n,p in enumerate(phases) if current is not None and n > current for i in p.activate_atom_ids}


def display_plan(source, target, actions, event, released, inactive, config, state_only):
    active, pairs = action_sites(actions)
    changed, new_pairs = differences(source, target)
    active |= changed
    pairs |= new_pairs
    event = raw(event) or {}
    for role_key in ("bindings","general_roles","metal_roles"):
        for value in (event.get(role_key) or {}).values():
            active.update(i for i in (value if isinstance(value,list) else [value]) if type(i) is int)
    for tx in event.get("electron_transfers",[]):
        active.update((tx["donor_atom_id"],tx["acceptor_atom_id"]))
        pairs.add((tx["donor_atom_id"],tx["acceptor_atom_id"]))
    for v in (source,target):
        if v is None:
            continue
        ligand_roles=set((event.get("ligand_roles") or {}).values())
        binding_roles=set((event.get("binding_roles") or {}).values())
        for ligand in getattr(v,"ligands",[]):
            if ligand.ligand_id in ligand_roles:active.update(ligand.donor_atom_ids)
        for b in getattr(v,"metal_bindings",[]):
            if b.binding_id in binding_roles:
                active.update((b.metal_atom_id,*b.donor_atom_ids))
    released = {frozenset(r) for r in [*released,*event.get("released_fragments",[])] if r}
    leaving = (event.get("bindings") or {}).get("leaving_group")
    detached_from=set()
    if type(leaving) is int and event.get("event_type") in {"SN2","E2","E1_IONIZATION","E1CB_ELIMINATION"}:
        output_edges={frozenset((b.a,b.b)) for b in target.bonds} if target else set()
        detached_from={a["a"] for a in map(raw,actions) if a.get("op")=="SPLIT" and
                       a.get("target")==leaving and frozenset((a["a"],leaving)) not in output_edges}
    inactive = set(inactive)
    plan = {}
    for side,view in (("input",source),("output",target)):
        if view is None:
            continue
        rows = []
        atoms = {a.atom_id:a for a in view.atoms}
        for group in components(view):
            if group <= inactive:
                reason = "not_yet_added"
            elif side == "output" and (frozenset(group) in released or
                    (type(leaving) is int and leaving in group and detached_from and group.isdisjoint(detached_from))):
                reason = "released_byproduct"
            elif group & active:
                reason = "participant"
            elif state_only or not active:
                reason = "state_snapshot" if state_only else "participation_unknown"
            else:
                reason = "unchanged_inventory"
            display = config.material_display == "all" or reason in {"participant","state_snapshot","participation_unknown"}
            counts = Counter(atoms[i].element for i in group)
            counts["H"] += sum(atoms[i].hydrogens for i in group)
            names = sorted(counts, key=lambda s: (0 if s=="C" else 1 if s=="H" else 2,s))
            formula = "".join(s + (str(counts[s]) if counts[s]>1 else "") for s in names if counts[s])
            q = [atoms[i].formal_charge for i in group]
            rows.append(dict(atom_ids=sorted(group),formula=formula,charge=None if None in q else sum(q),
                             reason=reason,displayed=display))
        plan[side] = rows
    return plan, active, sorted(pairs)


def subset(view, ids):
    """An ephemeral rendering view. The executor only sees the original view."""
    if ids == {a.atom_id for a in view.atoms}:
        return view
    fields = dict(atoms=[a for a in view.atoms if a.atom_id in ids],
                  bonds=[b for b in view.bonds if b.a in ids and b.b in ids])
    for name in ("metal_centers","reactive_centers"):
        if hasattr(view,name):
            fields[name] = [m for m in getattr(view,name) if m.atom_id in ids]
    if hasattr(view,"ligands"):
        fields["ligands"] = [l.model_copy(update={"atom_ids":[i for i in l.atom_ids if i in ids]})
                              for l in view.ligands if set(l.donor_atom_ids) <= ids]
        fields["metal_bindings"] = [b for b in view.metal_bindings
                                    if b.metal_atom_id in ids and set(b.donor_atom_ids) <= ids]
    return view.model_copy(update=fields)


def pool_details(folder,panel,prefix=""):
    """Offline HTML: open the inventory or look up any folded atom."""
    import json
    from html import escape
    if not panel.get("inventory"):
        return ""
    file=folder/panel["inventory"]
    data=json.loads(file.read_text(encoding="utf-8"))
    names={"participant":"参与","unchanged_inventory":"旁观／留存","released_byproduct":"生成后入库",
           "not_yet_added":"后续阶段加入","state_snapshot":"状态快照","participation_unknown":"参与范围待定"}
    body=[]
    for key, side in data["sides"].items():
        body.append("<b>"+("反应前" if key=="input" else "反应后")+"</b><ul>")
        for c in side["components"]:
            charge=c["charge"]
            q="（电荷未知）" if charge is None else (f" ({charge:+d})" if charge else "")
            ids=", ".join(map(str,c["atom_ids"]))
            body.append(f'<li>{escape(c["formula"]+q)} · {names[c["reason"]]} · 原子 {ids}</li>')
        for a in side["abbreviations"]:
            anchors=", ".join(map(str,a.get("anchor_atom_ids",[a["anchor_atom_id"]])))
            body.append(f'<li>{escape(a["label"])} · 锚点 {anchors} · 原子 {escape(str(a["atom_ids"]))} · {escape(a["source"])}</li>')
        body.append("</ul>")
    return ('<details class="material-pool"><summary>物料区与缩写对照</summary>'+''.join(body)+
            '<a href="'+escape(prefix+file.name)+'">完整映射结构与电子账本 JSON</a></details>')
