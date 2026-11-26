import json
from collections import Counter

from ..chemistry import FormalError, Unsupported, order
from ..engine import apply_actions
from ..v02.ast import Fishhook
from ..v02.electrons import apply_fishhooks
from ..v02.execution import product_spin as product_spin02
from ..v03.execution import execute_event as execute03
from .ast import AryneRegion, ChemicalState, ElectronicState, ReactiveCenter
from .ledger import PROFILES, occupation_profiles, validate_ledger, validate_view
from .templates import compile_template


def execute_event(source,event,*,_compiler=compile_template,_validator=validate_view,_ledger=validate_ledger,_charge=None):
    validate_view=_validator
    validate_view(source)
    if event.reactant_view_id!=source.view_id or event.product_state_id==source.state_id:
        raise FormalError("Reaction must use its input view and create a new state")
    if event.event_type=="ISC":
        raise Unsupported("requires_spin_transition: dedicated ISC request retained; spin dynamics execution unsupported")
    if not event.reactive_template:
        if source.reactive_centers or source.aryne_regions:
            raise Unsupported("Reactive intermediate events require a registered reactive_template")
        target,report=execute03(source,event)
    else:
        actions,created,aryne,report=_compiler(source,event)
        key=lambda a:json.dumps(a.model_dump(mode="json"),sort_keys=True)
        if event.actions is not None and Counter(map(key,event.actions))!=Counter(map(key,actions)):
            raise FormalError("Explicit action group disagrees with registered template (including electron slots)")
        event.actions=actions
        if all(isinstance(a,Fishhook) for a in actions):
            target=apply_fishhooks(source,actions,event.product_state_id,event.product_view_id,_validator=_ledger,**({'_charge':_charge} if _charge else {}))
        elif any(isinstance(a,Fishhook) for a in actions):
            raise Unsupported("Mixed pair/fishhook reactive templates are not registered")
        else:
            target=apply_actions(source,actions,event.product_state_id,event.product_view_id,_validator=_ledger,**({'_charge':_charge} if _charge else {}))
        actual=occupation_profiles(target)
        expected={c.atom_id:c.profile for c in source.reactive_centers if actual.get(c.atom_id)==c.profile}
        expected.update(created)
        if actual!=expected: raise FormalError("Product reactive occupations differ from template constraints")
        target.reactive_centers=[ReactiveCenter(center_id=f"center.{i}",atom_id=i,species=PROFILES[p][0],profile=p,evidence=event.evidence,
            existence_basis="formal_model",reaction_context=report["scope"]) for i,p in sorted(actual.items())]
        target.aryne_regions=[r.model_copy(deep=True) for r in source.aryne_regions if order(target,*r.aryne_atom_ids)==3]
        if aryne: target.aryne_regions.append(AryneRegion(**aryne,evidence=event.evidence,existence_basis="formal_model",reaction_context=report["scope"]))
    r=sum(a.radical_electrons for a in target.atoms)
    if source.electronic_state.coupling=="high_spin" and r!=sum(a.radical_electrons for a in source.atoms):
        raise Unsupported("requires_spin_transition: cannot remove high-spin singly occupied centers by direct closure")
    # Record local occupation character, never infer total multiplicity from r.
    if source.electronic_state.character!="unknown":
        target.electronic_state=ElectronicState(character="open_shell" if r else "closed_shell",
            coupling=source.electronic_state.coupling if r else "closed_shell",
            label=source.electronic_state.label if r else "closed-shell formal occupation",
            evidence=source.electronic_state.evidence)
    report.update(event_id=event.event_id,event_type=event.event_type,reactive_template=event.reactive_template,
        formal_status="passed",atomic_commit=True,ledger=validate_view(target),human_approved=False,
        compiled_actions=[a.model_dump(mode="json") for a in event.actions or []],
        product_species=species_inventory(target))
    return target,report


def species_inventory(view):
    from ...research.search import components
    groups=components(view)
    return [{"atom_ids":sorted(i for i,g in groups.items() if g==group),
             "formal_charge":sum(a.formal_charge for a in view.atoms if groups[a.atom_id]==group)}
            for group in sorted(set(groups.values()))]


def product_spin(origin,event,target):
    name=event.reactive_template or ""
    if name.startswith("TRIPLET") and origin.multiplicity!=3:
        raise Unsupported("Triplet template requires declared total triplet for this one-active-center scope")
    if name=="CARBYNE_DOUBLET_HAT" and origin.multiplicity!=2:
        raise Unsupported("Carbyne HAT requires explicit total doublet assumption")
    # Legacy product_spin checks parity/range and prohibits multiplicity changes.
    return product_spin02(origin,event,target)


def make_state(origin,event,target,multiplicity,report):
    return ChemicalState(state_id=target.state_id,primary_view_id=target.view_id,total_charge=origin.total_charge,
        multiplicity=multiplicity,electronic_mode="RADICAL" if any(a.radical_electrons for a in target.atoms) else "CLOSED_SHELL",
        electronic_state=target.electronic_state,conditions={**origin.conditions,**event.conditions})
