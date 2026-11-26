"""Versioned aza profiles; reuse ordered-topology compilation and native ledger.

These are formal six-electron models, not eligibility/kinetic predictions.  Local
occupation and topology are checked on both sides; spectator chemistry is kept.
"""
from collections import Counter
import json
from ..chemistry import FormalError, Unsupported, atoms_by_id, bond_valences
from ..v08.pericyclic import compile_profile

VERSION = '0.9-r6'
RULE_VERSION = 'umr-sigmatropic/0.9-r6'
PROFILES = {
    'AZA_CLAISEN_3_3': ('33', ('CCNCCC',), 'neutral N-allyl enamine to imine; free, closed-shell 3-aza-1,5-hexadiene'),
    'AZA_COPE_3_3': ('33', ('CNCCCC',), 'cationic 2-aza-Cope, iminium N(+); no neutral/diaza or Mannich inference'),
    'WITTIG_N_2_3': ('23', ('CNCCC',), 'alpha-carbanion attached to neutral tertiary N; amide anion retained; not ammonium ylide'),
}
SOURCES = {
    'AZA_CLAISEN_3_3': {'book_pdf_page': 74, 'book_printed_page': 20,
        'primary_url': 'https://pubs.acs.org/doi/10.1021/jo00019a021',
        'limit': 'Primary paper is charge-accelerated scope; neutral thermal scope comes from textbook, not that experiment.'},
    'AZA_COPE_3_3': {'book_pdf_page': 76, 'book_printed_page': 22,
        'primary_url': 'https://pubs.acs.org/doi/10.1021/jo00221a012',
        'limit': 'Finite iminium rearrangement, not certification of a full alkaloid synthesis.'},
    'WITTIG_N_2_3': {'book_pdf_page': 80, 'book_printed_page': 26,
        'primary_url': 'https://pubs.acs.org/doi/10.1021/jo0056343',
        'limit': 'Specific substrate activation and metal aggregation need case evidence; no general chirality assignment.'},
}


def occupation(view, event, *, product=False):
    """Check the exact local Lewis model; N includes all unpaired electrons."""
    spec = PROFILES[event.template_id]
    ids = event.bindings.get('ordered_atoms')
    if (set(event.bindings) != {'ordered_atoms'} or not isinstance(ids, list)
            or any(type(i) is not int for i in ids) or len(ids) != len(spec[1][0])
            or len(ids) != len(set(ids))):
        raise FormalError('Exact distinct integer ordered_atoms roles required')
    atoms = atoms_by_id(view)
    if not set(ids) <= atoms.keys():
        raise FormalError('Missing aza profile atom reference')
    if ''.join(atoms[i].element for i in ids) not in spec[1]:
        raise Unsupported('Element ordering outside registered aza profile')
    if any(atoms[i].radical_electrons != 0 for i in ids):
        raise Unsupported('Aza profiles require known closed-shell local occupations')
    if any(c.atom_id in ids for c in view.reactive_centers):
        raise Unsupported('Special reactive center needs its own sigmatropic profile')
    if any(set(b.donor_atom_ids) & set(ids) for b in view.metal_bindings):
        raise Unsupported('Coordinated aza participants require a separately evidenced profile')
    if any({b.donor_atom_id,b.acceptor_atom_id} & set(ids) for b in view.lewis_adducts):
        raise Unsupported('Lewis-acid-bound aza participants require a separate profile')
    if any(set(b.atom_ids) & set(ids) and (b.direction!='inside' or b.domain in
            {'bulk_metal','surface','external_electron_reservoir'}) for b in view.material_boundaries):
        raise Unsupported('External/bulk/surface aza participant needs a dedicated exchange model')
    # Test against the forward starting/ending occupation, including reverse use.
    end = (event.direction == 'reverse') != product
    if spec[0] == '33':
        edges = [(0,1,1),(1,2,2),(3,4,2),(4,5,1),(0,5,1)] if end else [(0,1,2),(1,2,1),(2,3,1),(3,4,1),(4,5,2)]
    else:
        edges = [(0,1,1),(2,3,2),(3,4,1),(0,4,1)] if end else [(0,1,1),(1,2,1),(2,3,1),(3,4,2)]
    required = {tuple(sorted((ids[a],ids[b]))):n for a,b,n in edges}
    actual = {tuple(sorted((b.a,b.b))):b.order for b in view.bonds if b.a in ids and b.b in ids}
    if actual != required:
        raise FormalError('Aza ordered fragment bond topology mismatch (including extra participant bonds)')
    expected = {i: (0, 0, 4) for i in ids}  # charge, nonbonding electrons, bond-order sum
    if event.template_id == 'AZA_CLAISEN_3_3':
        expected[ids[2]] = (0, 2, 3)
    elif event.template_id == 'AZA_COPE_3_3':
        expected[ids[1]] = (1, 0, 4)
    else:
        expected[ids[0]] = (0, 0, 4) if end else (-1, 2, 3)
        expected[ids[1]] = (-1, 4, 2) if end else (0, 2, 3)
    valences = bond_valences(view)
    for i, (q, nonbonding, valence) in expected.items():
        a = atoms[i]
        if (a.formal_charge, a.nonbonding_electrons, valences[i]) != (q, nonbonding, valence):
            raise Unsupported(f'Aza occupation mismatch at atom {i}: requires q={q}, N={nonbonding}, b={valence}')


def compile_pericyclic(source, event):
    if event.template_id not in PROFILES:
        raise Unsupported('Unregistered aza profile')
    if event.template_version != VERSION:
        raise Unsupported('Aza profiles require template_version ' + VERSION)
    if event.event_type != 'PERICYCLIC' or event.execution_level != 'electron_actions':
        raise FormalError('Aza profile requires a native PERICYCLIC electron-actions event')
    if any([event.general_event, event.rule_template, event.element_template,
            event.metal_template, event.reactive_template, event.electron_transfers]):
        raise FormalError('Conflicting selectors on aza event')
    if not event.evidence or event.mechanism_assumption != 'concerted':
        raise FormalError('Aza event requires evidence and an explicit concerted model declaration')
    if event.electronic_outcome is not None:
        raise Unsupported('Sigmatropic profile cannot declare an electronic-state transition')
    if event.rotation_mode != 'unknown':
        raise Unsupported('Electrocyclic rotation labels do not specify aza sigmatropic geometry')
    for c in event.stereo_constraints:
        if c.get('kind') != 'facial':
            raise Unsupported('Aza profile does not assign absolute/relative stereochemical constraints')
    occupation(source, event)
    actions, report = compile_profile(source, event, profiles=PROFILES, template_version=VERSION,
        rule_version=RULE_VERSION, neutral_linkers={'WITTIG_N_2_3'})
    report.update(source=SOURCES[event.template_id], participating_electrons=6,
        product_selection='not_predicted', electronic_state_transition=False)
    return actions, report


def execute(source, event):
    from .execution import transact
    from .ledger import validate_view
    from .ast import RULE_VERSION as rules
    from ..v06.stereo import apply_event
    validate_view(source)
    actions, report = compile_pericyclic(source, event)
    key = lambda a: json.dumps(a.model_dump(mode='json'), sort_keys=True)
    if event.actions is not None and Counter(map(key, event.actions)) != Counter(map(key, actions)):
        raise FormalError('Declared actions differ from aza profile compilation')
    # Native transaction checks identity, H, charge, electrons and duplicate use.
    target = transact(source, actions, event.product_state_id, event.product_view_id)
    occupation(target, event, product=True)
    event.actions = actions
    trial = event.model_copy(deep=True)
    trial.stereo_effects = []
    report.update(event_id=event.event_id, event_type=event.event_type,
        profile_rule_version=RULE_VERSION, rule_version=rules,
        atomic_commit=True, formal_status='passed', execution_mode='electron_actions',
        compiled_actions=[a.model_dump(mode='json') for a in actions],
        stereochemistry=apply_event(source, target, trial),
        unresolved_stereo_effects=[s.model_dump(mode='json') for s in event.stereo_effects],
        chemical_feasibility='unknown')
    return target, report
