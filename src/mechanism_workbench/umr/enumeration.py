"""Shared finite resonance frontier. Priority never crosses BFS depth."""
from heapq import heappop, heappush


class ResonanceFrontier:
    def __init__(self, strategy="bfs"):
        if strategy not in {"bfs", "charge_priority"}:
            raise ValueError("Unknown resonance enumeration strategy")
        self.strategy, self.pending, self.serial = strategy, [], 0

    def push(self, view, depth):
        charge = sum(abs(a.formal_charge or 0) for a in view.atoms)
        key = (depth, charge) if self.strategy == "bfs" else (charge, depth)
        heappush(self.pending, (*key, self.serial, depth, view))
        self.serial += 1

    def pop(self):
        item = heappop(self.pending)
        return item[-2:]

    def __bool__(self):
        return bool(self.pending)
