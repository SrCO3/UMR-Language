def components(view):
    neighbors = {a.atom_id: set() for a in view.atoms}
    for b in view.bonds:
        neighbors[b.a].add(b.b); neighbors[b.b].add(b.a)
    result = {}
    for start in sorted(neighbors):
        if start in result:
            continue
        todo = [start]
        while todo:
            i = todo.pop()
            if i in result:
                continue
            result[i] = start
            todo.extend(neighbors[i] - result.keys())
    return result
