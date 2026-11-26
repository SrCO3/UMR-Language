"""Event-local electron slots and complete ordinary two-electron bond pools."""

from collections import Counter, defaultdict

from ..ast import Bond
from ..chemistry import (FormalError, atoms_by_id, bond_orders, bond_valences, derived_charge,
                         electron_count, identity)
from .ast import AtomDestination, BondSource, Fishhook, RadicalSource
from .ledger import validate_view

PAIRS = {"sigma": 1, "pi1": 2, "pi2": 3}


def edge(a, b):
    if a == b:
        raise FormalError("Electron source/pool cannot be a self bond")
    return tuple(sorted((a, b)))


def fish_key(action):
    src, dst = action.source, action.destination
    s = ('lone_pair',src.atom_id,src.pair_index,src.slot) if src.kind=='lone_pair' else ("radical", src.atom_id, src.slot) if isinstance(src, RadicalSource) else ("bond", *edge(src.a, src.b), src.pair, src.slot)
    d = ('lone_pair',dst.atom_id) if dst.kind=='lone_pair' else ("atom", dst.atom_id) if isinstance(dst, AtomDestination) else ("bond", *edge(dst.a, dst.b), dst.pair)
    return s, d


def apply_fishhooks(source, actions, state_id, view_id, *, preserve_stereo=False,
                    _validator=validate_view, _charge=derived_charge, _allow_lone_pairs=False):
    _validator(source)
    atoms, orders = atoms_by_id(source), bond_orders(source)
    used, bond_slots, pools = set(), defaultdict(set), Counter()
    radical_debits, received, touched = Counter(), Counter(), set()
    lone_slots=defaultdict(set); paired_received=Counter()
    for action in actions:
        if not isinstance(action, Fishhook):
            raise FormalError("Single-electron executor accepts only FISHHOOK actions")
        src, dst = action.source, action.destination
        if not _allow_lone_pairs and (src.kind=='lone_pair' or dst.kind=='lone_pair'):
            raise FormalError('Explicit nonbonding pair slots require the v0.7 adapter')
        key, _ = fish_key(action)
        if key in used:
            raise FormalError(f"Duplicate electron slot consumption: {key}")
        used.add(key)
        if src.kind=='lone_pair':
            if src.atom_id not in atoms or src.pair_index >= (atoms[src.atom_id].nonbonding_electrons-atoms[src.atom_id].radical_electrons)//2:
                raise FormalError('Lone-pair source slot is not occupied')
            lone_slots[(src.atom_id,src.pair_index)].add(src.slot);touched.add(src.atom_id)
        elif isinstance(src, RadicalSource):
            if src.atom_id not in atoms or src.slot >= atoms[src.atom_id].radical_electrons:
                raise FormalError("Radical source slot is not occupied in the input")
            radical_debits[src.atom_id] += 1
            touched.add(src.atom_id)
        else:
            pair = edge(src.a, src.b)
            if orders.get(pair, 0) < PAIRS[src.pair]:
                raise FormalError("Requested sigma/pi electron pair does not exist in the input")
            bond_slots[(*pair, src.pair)].add(src.slot)
            touched.update(pair)
        if dst.kind=='lone_pair':
            if dst.atom_id not in atoms:raise FormalError('Lone-pair destination missing')
            paired_received[dst.atom_id]+=1;touched.add(dst.atom_id)
        elif isinstance(dst, AtomDestination):
            if dst.atom_id not in atoms:
                raise FormalError("Destination atom is missing")
            received[dst.atom_id] += 1
            touched.add(dst.atom_id)
        else:
            pair = edge(dst.a, dst.b)
            if any(i not in atoms for i in pair):
                raise FormalError("Destination bond references missing atoms")
            pools[(*pair, dst.pair)] += 1
            touched.update(pair)
    if any(slots != {0, 1} for slots in bond_slots.values()):
        raise FormalError("Both slots of each consumed bond pair require destinations")
    if any(count != 2 for count in pools.values()):
        raise FormalError("A new ordinary bond pool must receive exactly two electrons")
    if any(count%2 for count in paired_received.values()):raise FormalError('A nonbonding pair destination requires complete pairs')
    consumed, formed = defaultdict(set), defaultdict(set)
    for a, b, name in bond_slots:
        consumed[(a,b)].add(PAIRS[name])
    for a, b, name in pools:
        formed[(a,b)].add(PAIRS[name])
    result = source.model_copy(deep=True)
    result.state_id, result.view_id = state_id, view_id
    result.bonds = []
    old = {edge(b.a,b.b): b for b in source.bonds}
    for pair in sorted(set(orders) | set(formed)):
        retained = set(range(1, orders.get(pair,0)+1)) - consumed[pair]
        if retained & formed[pair]:
            raise FormalError("A destination pool tries to create an already occupied bond pair")
        final = retained | formed[pair]
        if final != set(range(1,len(final)+1)):
            raise FormalError("Cannot retain pi electrons without a sigma bond or skip a bond-pair level")
        if final:
            bond = old[pair].model_copy(deep=True) if pair in old else Bond(a=pair[0],b=pair[1],order=len(final))
            bond.order = len(final)
            if not preserve_stereo and any(i in touched for i in pair):
                bond.stereo = {"status":"unknown","reason":"radical event stereochemistry not assigned"}
            result.bonds.append(bond)
    for atom in result.atoms:
        change = received[atom.atom_id]-radical_debits[atom.atom_id]
        spent=sum(len(slots) for (i,_),slots in lone_slots.items() if i==atom.atom_id)
        unpaired=sum(len(slots)==1 for (i,_),slots in lone_slots.items() if i==atom.atom_id)
        atom.nonbonding_electrons += change-spent+paired_received[atom.atom_id]
        atom.radical_electrons += change+unpaired
        if not preserve_stereo and atom.atom_id in touched:
            atom.stereo = {"status":"unknown","reason":"radical event stereochemistry not assigned"}
    for atom in result.atoms:
        atom.formal_charge = _charge(atom, bond_valences(result)[atom.atom_id])
    from ..stereochemistry import reassess_centers
    reassess_centers(source,result,touched,resonance=preserve_stereo)
    _validator(result)
    if identity(source) != identity(result) or electron_count(source) != electron_count(result):
        raise FormalError("Atom/isotope/H or electron conservation failed")
    if sum(a.formal_charge for a in source.atoms) != sum(a.formal_charge for a in result.atoms):
        raise FormalError("Total charge conservation failed")
    return result


def atom_to_bond(atom, a, b, pair="sigma"):
    from .ast import BondDestination
    return Fishhook(source=RadicalSource(atom_id=atom), destination=BondDestination(a=a,b=b,pair=pair))


def bond_to_atom(a,b,pair,slot,destination):
    return Fishhook(source=BondSource(a=a,b=b,pair=pair,slot=slot), destination=AtomDestination(atom_id=destination))


def bond_to_bond(a,b,pair,slot,c,d,target_pair="sigma"):
    from .ast import BondDestination
    return Fishhook(source=BondSource(a=a,b=b,pair=pair,slot=slot), destination=BondDestination(a=c,b=d,pair=target_pair))
