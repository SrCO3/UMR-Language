"""Extend only registered element profiles; keep the existing main-group ledger."""
from ..v07.ledger import validate_view as old_validate, signature, compare_views
from ..v05.ledger import validate_ledger as ledger05
from ..v04.ledger import validate_view as validate04
from ..v07.ledger import project_symbols, symbols
from ..chemistry import Unsupported
from .elements import NEW_METALS, METAL_GROUP, MAIN, profile, coordination, oxidation


def validate_ledger(view):
    return validate_view(view)


def validate_view(view,state=None):
    if not any(a.element in NEW_METALS for a in view.atoms):return old_validate(view,state)
    if symbols(view):raise Unsupported('New-metal plus opaque-symbol composite profile not yet registered')
    report=ledger05(view,_profile=profile,_coordination=coordination,_main=MAIN,_metals=METAL_GROUP,_oxidation=oxidation)
    shadow=view.model_copy(deep=True)
    unknown=[a for a in shadow.atoms if a.radical_electrons is None]
    if unknown:
        if any(a.element not in METAL_GROUP for a in unknown):raise Unsupported('Unknown main-group local occupation')
        if view.electronic_state.character!='unknown' or view.electronic_state.coupling!='unknown' or state and state.multiplicity is not None:
            raise Unsupported('Unknown metal coupling cannot establish a global spin declaration')
        for a in unknown:a.radical_electrons=0
    checked={**report,'unpaired_electrons':sum(a.radical_electrons for a in shadow.atoms)}
    validate04(shadow,state,_ledger=lambda v:checked)
    report.update(rule_version='umr-elements/0.8',spin_status='unknown_not_inferred' if unknown else 'declared_local_occupation')
    return report
