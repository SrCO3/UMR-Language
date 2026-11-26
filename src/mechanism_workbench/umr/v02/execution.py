from collections import Counter

from ..ast import ReactionEvent as Event01
from ..chemistry import FormalError, Unsupported, atoms_by_id, bond_orders, identity
from ..engine import action_key, apply_actions, check_event_rule
from .ast import BondSource, Fishhook
from .electrons import apply_fishhooks, fish_key
from .ledger import radical_count, validate_spin, validate_view
from .pericyclic import compile_pericyclic
from .radical import compile_radical


def execute_event(source,event):
    if event.reactant_view_id!=source.view_id or event.product_state_id==source.state_id:
        raise FormalError("Event input identity or new-state identity is invalid")
    if event.event_type in {"SET","ISC","COORDINATION"}:
        raise Unsupported(f"{event.event_type} execution remains unsupported")
    if event.event_type in {"RADICAL","PERICYCLIC"}:
        actions,report = (compile_radical if event.event_type=="RADICAL" else compile_pericyclic)(source,event)
        if event.actions is not None:
            if not event.actions:
                raise FormalError("An explicit empty DO group is not a template omission")
            def semantic_key(a):
                source_key,destination_key=fish_key(a)
                return (source_key[:-1] if isinstance(a.source,BondSource) else source_key,destination_key)
            key = semantic_key if event.event_type=="RADICAL" else action_key
            if any(isinstance(a,Fishhook)!=(event.event_type=="RADICAL") for a in event.actions):
                raise FormalError("Single/pair electron action mode disagrees with event classification")
            if Counter(key(a) for a in actions)!=Counter(key(a) for a in event.actions):
                raise FormalError("Explicit DO arrows disagree with the declared template")
        else:
            event.actions = actions
        target = apply_fishhooks(source,event.actions,event.product_state_id,event.product_view_id) if event.event_type=="RADICAL" else apply_actions(source,event.actions,event.product_state_id,event.product_view_id,_validator=validate_view)
    else:
        if event.template_id or any(isinstance(a,Fishhook) for a in event.actions or []):
            raise FormalError("Explicit RADICAL/PERICYCLIC classification is required for those templates/actions")
        if not event.actions or any(type(v) is not int for v in event.bindings.values()):
            raise FormalError("Polar events require explicit paired actions and scalar bindings")
        base_fields = Event01.model_fields
        data = {k:v for k,v in event.model_dump(mode="json").items() if k in base_fields}
        base = Event01.model_validate(data)
        report = check_event_rule(source,base)
        target = apply_actions(source,base.actions,event.product_state_id,event.product_view_id,_validator=validate_view)
        report.update({"formal_execution":"passed","template_scope":"matched","orbital_symmetry":"not_applicable",
                       "geometry_status":"unverified","stereochemistry":"unknown",
                       "mechanism_evidence":event.mechanism_evidence,"kinetics_selectivity":"unknown"})
    region = set(event.region_atom_ids)
    if region:
        before,after = bond_orders(source),bond_orders(target)
        if any(any(i not in region for i in pair) for pair in set(before)|set(after) if before.get(pair,0)!=after.get(pair,0)):
            raise FormalError("Template changed a bond outside its declared region")
        source_atoms,target_atoms = atoms_by_id(source),atoms_by_id(target)
        if any(source_atoms[i].model_dump()!=target_atoms[i].model_dump() for i in source_atoms if i not in region):
            raise FormalError("Template changed an atom outside its declared region")
    report.update({"event_id":event.event_id,"event_type":event.event_type,"template_id":event.template_id,
                   "template_version":event.template_version,"atomic_commit":True,"formal_status":"passed",
                   "ledger":validate_view(target),"unpaired_before":radical_count(source),
                   "unpaired_after":radical_count(target),"human_approved":False})
    return target,report


def product_spin(source_state,event,target):
    initial = source_state.multiplicity
    if event.spin_correlation != "unknown":
        correlated = 1 if event.spin_correlation=="singlet" else 3
        if initial is not None and initial!=correlated:
            raise Unsupported("requires_spin_transition: spin correlation contradicts supplied total multiplicity")
        initial = correlated
    final = event.product_multiplicity if event.product_multiplicity is not None else initial
    if initial is not None and final is not None and initial!=final:
        raise Unsupported("requires_spin_transition: ISC is not implemented")
    validate_spin(target,final)
    return final,{"input_multiplicity":source_state.multiplicity,"effective_correlation":initial,
                  "output_multiplicity":final,"requires_spin_transition":False,
                  "status":"unknown" if final is None else "compatible_local_occupancy_not_spin_dynamics"}
