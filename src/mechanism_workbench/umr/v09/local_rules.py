"""Author-declared Lewis occupations and structural tests on a named Lewis view."""
from contextlib import contextmanager
from contextvars import ContextVar
from typing import Annotated

from pydantic import Field, StringConstraints, model_validator

from ..ast import AtomID, Identifier, Node
from ..chemistry import FormalError, atoms_by_id, bond_orders, bond_valences


Text = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]
Count = Annotated[int, Field(ge=0)]


class ViewRule(Node):
    rule_id: Identifier
    view_id: Identifier
    description: Text
    source: Text | None = None


class ReagentDeclaration(Node):
    reagent_id: Identifier
    view_id: Identifier
    description: Text
    atom_ids: list[AtomID] = Field(min_length=1)
    role: Text
    phase_id: Identifier | None = None
    source: Text | None = None

    @model_validator(mode='after')
    def distinct_atoms(self):
        if len(self.atom_ids) != len(set(self.atom_ids)):
            raise ValueError('Duplicate atom in reagent declaration')
        return self


class ValenceRule(ViewRule):
    atom_id: AtomID
    element: Text
    bond_order_sum: Count
    nonbonding_electrons: Count
    radical_electrons: Count | None
    formal_charge: int
    oxidation_state: int | None = None
    occupied_sites: Count | None = None


class AtomTest(Node):
    atom_id: AtomID
    element: Text | None = None
    formal_charge: int | None = None
    bond_order_sum: Count | None = None
    nonbonding_electrons: Count | None = None
    radical_electrons: Count | None = None

    @model_validator(mode='after')
    def has_test(self):
        if not self.model_dump(exclude_none=True, exclude={'atom_id'}):
            raise ValueError('An atom test needs at least one property')
        return self


class BondTest(Node):
    a: AtomID
    b: AtomID
    order: Annotated[int, Field(ge=0, le=3)]

    @model_validator(mode='after')
    def distinct_atoms(self):
        if self.a == self.b:
            raise ValueError('A bond test needs two distinct atoms')
        return self


class IntermediateRule(ViewRule):
    atoms: list[AtomTest] = Field(default_factory=list)
    bonds: list[BondTest] = Field(default_factory=list)
    total_charge: int | None = None

    @model_validator(mode='after')
    def has_tests(self):
        if not (self.atoms or self.bonds or self.total_charge is not None):
            raise ValueError('An intermediate rule needs a structural or charge test')
        if len({a.atom_id for a in self.atoms}) != len(self.atoms):
            raise ValueError('Combine tests on the same atom into one entry')
        if len({tuple(sorted((b.a,b.b))) for b in self.bonds}) != len(self.bonds):
            raise ValueError('Duplicate intermediate bond test')
        return self


def rules_by_view(document):
    grouped = {}
    for rule in (*getattr(document, 'valence_rules', ()), *getattr(document, 'intermediate_rules', ()),
                 *getattr(document, 'reagents', ())):
        grouped.setdefault(rule.view_id, []).append(rule)
    return grouped


# Transactions and resonance exploration call the same view validator. ContextVar
# keeps independent replays isolated, including when a quench runs recursively.
_active_rules = ContextVar('umr_view_rules', default={})


@contextmanager
def use_rules(document):
    from .elements import LEGACY, use_elements
    elements = LEGACY
    if getattr(document, 'schema_version', None) == 'umr/1.0':
        from ..v10.elements import RULES as elements
    token = _active_rules.set(rules_by_view(document))
    try:
        with use_elements(elements):
            yield
    finally:
        _active_rules.reset(token)


def active_rules(view_id):
    return [item for item in _active_rules.get().get(view_id, ()) if isinstance(item, ViewRule)]


def active_reagents(view_id):
    return [item for item in _active_rules.get().get(view_id, ()) if isinstance(item, ReagentDeclaration)]


def check_reagents(view, reagents, symbolic_ids):
    known = set(atoms_by_id(view)).difference(symbolic_ids)
    contacts = [{b.a, b.b} for b in view.bonds]
    ligands = {ligand.ligand_id: ligand.atom_ids for ligand in view.ligands}
    for binding in view.metal_bindings:
        contacts.append({binding.metal_atom_id, *ligands.get(binding.ligand_id, binding.donor_atom_ids)})
    for reagent in reagents:
        members = set(reagent.atom_ids)
        if not members <= known:
            raise FormalError(f'{reagent.reagent_id}: missing or symbolic reagent atoms {sorted(members - known)}')
        if any(members & group and group - members for group in contacts):
            raise FormalError(f'{reagent.reagent_id}: reagent declaration splits a bonded or coordinated component')


def require_equal(rule, field, expected, actual):
    if expected != actual:
        raise FormalError(f'{rule.rule_id} ({rule.view_id}): {field} expected {expected!r}, observed {actual!r}')


def check_local_rules(view, rules, *, symbolic_ids=()):
    """Return verified local occupations. Coordination counts are checked separately."""
    atoms, sums, orders = atoms_by_id(view), bond_valences(view), bond_orders(view)
    symbolic_ids = set(symbolic_ids)
    def atom(rule, atom_id):
        if atom_id not in atoms or atom_id in symbolic_ids:
            raise FormalError(f'{rule.rule_id} ({view.view_id}): missing or symbolic atom #{atom_id}')
        return atoms[atom_id]

    def check_atom(rule, item, fields):
        a = atom(rule, item.atom_id)
        for field, expected in fields.items():
            actual = sums[a.atom_id] if field == 'bond_order_sum' else getattr(a, field)
            require_equal(rule, f'#{a.atom_id} {field}', expected, actual)

    profiles = {}
    for rule in rules:
        if isinstance(rule, ValenceRule):
            fields = ('element', 'bond_order_sum', 'nonbonding_electrons', 'radical_electrons', 'formal_charge')
            check_atom(rule, rule, {k: getattr(rule, k) for k in fields})
            profiles[rule.atom_id] = rule
        else:
            for test in rule.atoms:
                check_atom(rule, test, test.model_dump(exclude={'atom_id'}, exclude_none=True))
            for bond in rule.bonds:
                atom(rule, bond.a)
                atom(rule, bond.b)
                order = orders.get(tuple(sorted((bond.a, bond.b))), 0)
                require_equal(rule, f'bond #{bond.a}–#{bond.b}', bond.order, order)
            if rule.total_charge is not None:
                if symbolic_ids:
                    raise FormalError(f'{rule.rule_id}: total charge requires a fully specified inventory')
                require_equal(rule, 'total_charge', rule.total_charge, sum(a.formal_charge for a in atoms.values()))
    return profiles
