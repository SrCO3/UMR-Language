"""Finite elemental-step templates. No free-form charge/oxidation-state patches."""
from collections import Counter

from ..ast import Pair, Shift, Split
from ..chemistry import FormalError, Unsupported, atoms_by_id, electron_count, identity, order
from ..engine import action_key, apply_actions
from ..v02.ast import Fishhook
from ..v02.electrons import apply_fishhooks
from .ast import MetalBinding
from .ledger import binding_checks, charge, metal_counts, validate_view

TEMPLATES = {"L_COORDINATION", "L_DISSOCIATION", "OXIDATIVE_ADDITION_2E", "REDUCTIVE_ELIMINATION_2E", "TRANSMETALATION_PD_ZN"}
DEFERRED = {"LIGAND_SUBSTITUTION", "CO_INSERTION", "ALKENE_INSERTION", "BETA_H_ELIMINATION",
            "PROTONOLYSIS", "METAL_SET", "SPIN_CROSSING", "MULTIMETAL_BRIDGE"}


def check_reserved(source, actions):
    _,_,reserved=binding_checks(source)
    spent=Counter(a.a for a in actions if isinstance(a,Pair))
    atoms=atoms_by_id(source)
    for i,n in spent.items():
        if 2*n+reserved[i]>atoms[i].nonbonding_electrons-atoms[i].radical_electrons:
            raise FormalError("Action double-spends a dative-reserved donor pair; dissociation must be explicit")
    eta_pairs={tuple(sorted(b.donor_atom_ids)) for b in source.metal_bindings if b.mode=="eta"}
    for a in actions:
        if isinstance(a,(Shift,Split)) and tuple(sorted((a.a,a.b))) in eta_pairs:
            raise Unsupported("Redistribution of a coordinated eta2 pi pair needs a metal-specific template")
        if isinstance(a,Fishhook):
            src=a.source
            if src.kind=="bond" and tuple(sorted((src.a,src.b))) in eta_pairs:
                raise Unsupported("Coordinated eta2 pi electron slots are reserved")


def apply_organic(source, actions, state_id, view_id, *, preserve_stereo=False):
    validate_view(source)
    mids={m.atom_id for m in source.metal_centers}
    for a in actions:
        data=a.model_dump()
        ids={v for k,v in data.items() if k in {"a","b","c"}}
        if isinstance(a,Fishhook):
            for end in (a.source,a.destination):
                ids|={v for k,v in end.model_dump().items() if k in {"a","b","atom_id"}}
        if ids&mids: raise Unsupported("Generic organic electron actions cannot touch a metal; use a scoped metal template")
    check_reserved(source,actions)
    if any(isinstance(a,Fishhook) for a in actions):
        return apply_fishhooks(source,actions,state_id,view_id,preserve_stereo=preserve_stereo,
                              _validator=validate_view,_charge=charge)
    return apply_actions(source,actions,state_id,view_id,preserve_stereo=preserve_stereo,
                         _validator=validate_view,_charge=charge)


def execute_event(source,event,*,_hooks=None):
    # Version adapters inject scoped validators without process-global mutation.
    hooks=_hooks or {}
    validate_view=hooks.get('validate_view',globals()['validate_view'])
    binding_checks=hooks.get('binding_checks',globals()['binding_checks'])
    metal_counts=hooks.get('metal_counts',globals()['metal_counts'])
    check_reserved=hooks.get('check_reserved',globals()['check_reserved'])
    charge=hooks.get('charge',globals()['charge'])
    validate_view(source)
    if source.view_id!=event.reactant_view_id or source.state_id==event.product_state_id:
        raise FormalError("Event source reference or distinct product state is invalid")
    if not event.metal_template:
        if event.metal_participation: raise Unsupported("Metal participation requires an explicit supported elemental template")
        if not source.metal_centers:
            from ..v02.execution import execute_event as execute02
            return execute02(source,event)
        if not event.actions: raise FormalError("An organic ReactionEvent requires a nonempty electron action group")
        if event.event_type in {"RADICAL","PERICYCLIC","SET","ISC","COORDINATION"}:
            raise Unsupported("Metal-containing radical/pericyclic reaction execution requires a future scoped adapter; organic resonance is supported")
        from ..engine import check_event_rule
        report=check_event_rule(source,event)
        result=apply_organic(source,event.actions or [],event.product_state_id,event.product_view_id)
        report.update(event_id=event.event_id,event_type=event.event_type,atomic_commit=True,formal_status="passed",ledger=validate_view(result))
        return result,report
    if not event.metal_participation: raise FormalError("Metal template requires metal_participation=true")
    template=event.metal_template
    if template not in TEMPLATES: raise Unsupported(f"Unsupported metal elemental template: {template}")
    expected_type="COORDINATION" if template.startswith("L_") else "POLAR"
    if event.event_type!=expected_type:
        raise Unsupported("This template is limited to declared closed-shell coordination or two-electron POLAR paths; radical/stepwise paths are separate events")
    if not event.mechanism_assumption or event.mechanism_evidence=="unknown":
        raise FormalError("Metal template requires an explicit mechanism assumption and evidence status; net change alone does not establish a path")
    if any(a.radical_electrons for a in source.atoms): raise Unsupported("Closed-shell metal templates cannot establish a radical pathway")
    centers,ligands,_=binding_checks(source)
    atoms=atoms_by_id(source)
    counts={m["metal_atom_id"]:m for m in metal_counts(source)}
    output=source.model_copy(deep=True)
    output.state_id,output.view_id=event.product_state_id,event.product_view_id
    bindings={b.binding_id:b for b in output.metal_bindings}
    actions=[]
    affected=set()
    geometry={"status":"unknown","reason":"2D connectivity does not establish spatial requirements"}

    def role(name):
        value=event.bindings.get(name)
        if not isinstance(value,int) or value not in atoms: raise FormalError(f"Missing metal template atom role {name}")
        return value

    def metal(name,element,os):
        i=role(name)
        allowed=hooks.get('metal_scope',{}).get((template,name),{(element,os)})
        if i not in centers or (atoms[i].element,counts[i]["oxidation_state"]) not in allowed:
            raise Unsupported(f"{name} requires {element}({os}) under known ionic ligand assignments")
        affected.add(i)
        return i

    def ligand(name):
        key=event.ligand_roles.get(name)
        if key not in ligands: raise FormalError(f"Missing ligand role {name}")
        return ligands[key]

    def binding(name):
        key=event.binding_roles.get(name)
        if key not in bindings: raise FormalError(f"Missing binding role {name}")
        return bindings[key]

    def add_binding(i,lig,mode):
        key=f"{event.event_id}.{lig.ligand_id}"
        site=event.site_roles.get(lig.ligand_id)
        if key in bindings: raise FormalError("Generated binding ID collides")
        bindings[key]=MetalBinding(binding_id=key,metal_atom_id=i,ligand_id=lig.ligand_id,
            donor_atom_ids=lig.donor_atom_ids,mode=mode,hapticity=2 if mode=="eta" else 1,site_id=site)

    if template=="L_COORDINATION":
        m=metal("metal","Pd",0) if counts.get(role("metal"),{}).get("oxidation_state")==0 else metal("metal","Pd",2)
        lig=ligand("incoming")
        if any(b.ligand_id==lig.ligand_id for b in bindings.values()): raise FormalError("Incoming ligand is already bound")
        if lig.ionic_class!="L": raise Unsupported("L coordination accepts only neutral 2e L ligands")
        mode="eta" if len(lig.donor_atom_ids)==2 else "dative"
        add_binding(m,lig,mode)
    elif template=="L_DISSOCIATION":
        old=binding("outgoing")
        m=role("metal")
        if m!=old.metal_atom_id or atoms[m].element!="Pd" or counts[m]["oxidation_state"] not in {0,2}:
            raise Unsupported("L dissociation requires the declared Pd(0/II) binding")
        affected.add(m)
        if old.mode not in {"dative","eta"}: raise Unsupported("Covalent X cleavage is not neutral L dissociation")
        del bindings[old.binding_id]
    elif template=="OXIDATIVE_ADDITION_2E":
        m=metal("metal","Pd",0); r,x=role("organic_atom"),role("leaving_atom")
        if len({m,r,x})!=3 or atoms[r].element!="C" or atoms[x].element not in {"Cl","Br","I"} or order(source,r,x)!=1:
            raise Unsupported("OA currently supports a single C-halogen bond into mononuclear Pd(0)")
        lr,lx=ligand("organic"),ligand("leaving")
        if lr.donor_atom_ids!=[r] or lx.donor_atom_ids!=[x]: raise FormalError("OA fragment donors do not match atom roles")
        if any(b.ligand_id in {lr.ligand_id,lx.ligand_id} for b in bindings.values()): raise FormalError("OA fragments were already metal-bound")
        actions=[Pair(a=m,b=r),Shift(a=r,b=x,c=m)]
        add_binding(m,lr,"sigma_covalent");add_binding(m,lx,"sigma_covalent")
    elif template=="REDUCTIVE_ELIMINATION_2E":
        m=metal("metal","Pd",2); one,two=binding("first"),binding("second")
        if one.binding_id==two.binding_id or any(b.metal_atom_id!=m or b.mode!="sigma_covalent" for b in (one,two)):
            raise FormalError("RE requires two different X sigma ligands on the same metal")
        r,x=one.donor_atom_ids[0],two.donor_atom_ids[0]
        if atoms[r].element!="C" or atoms[x].element not in {"C","Cl","Br","I"} or order(source,r,x):
            raise Unsupported("RE adapter supports C-C/C-halogen formation only")
        for relation in centers[m].relative_positions:
            if {relation.first_binding_id,relation.second_binding_id}=={one.binding_id,two.binding_id}:
                if relation.relation=="trans": raise FormalError("Declared trans ligands require isomerization before this cis RE template")
                if relation.relation=="cis": geometry={"status":"cis_declared","is_geometric_proof":False}
        actions=[Split(a=r,b=m,target=m),Shift(a=m,b=x,c=r)]
        del bindings[one.binding_id];del bindings[two.binding_id]
    else:
        pd=metal("acceptor","Pd",2);zn=metal("donor","Zn",2)
        leaving,organic=binding("acceptor_ligand"),binding("donor_ligand")
        if leaving.metal_atom_id!=pd or organic.metal_atom_id!=zn or any(b.mode!="sigma_covalent" for b in (leaving,organic)):
            raise FormalError("Transmetalation must identify explicit Pd-X and Zn-C bindings")
        x,r=leaving.donor_atom_ids[0],organic.donor_atom_ids[0]
        if atoms[x].element not in {"Cl","Br","I"} or atoms[r].element!="C": raise Unsupported("Only explicit neutral Pd-X / organozinc X-ligand exchange is supported")
        actions=[Shift(a=pd,b=x,c=zn),Shift(a=zn,b=r,c=pd)]
        del bindings[leaving.binding_id];del bindings[organic.binding_id]
        add_binding(zn,ligands[leaving.ligand_id],"sigma_covalent")
        add_binding(pd,ligands[organic.ligand_id],"sigma_covalent")
    if any(isinstance(a,Fishhook) for a in event.actions or []):
        raise Unsupported("Two-electron metal templates cannot be replaced by fishhook arrows")
    if event.actions is not None and Counter(action_key(a) for a in event.actions)!=Counter(action_key(a) for a in actions):
        raise FormalError("Explicit actions differ from the metal template's compiled electron group")
    if actions:
        check_reserved(source,actions)
        # Binding/OS metadata is updated before whole-result validation; no partial
        # graph is a state or an exported intermediate.
        raw=apply_actions(source,actions,event.product_state_id,event.product_view_id,
                          _validator=lambda v:None,_charge=charge)
        output.atoms,output.bonds=raw.atoms,raw.bonds
    output.metal_bindings=list(bindings.values())
    for m in output.metal_centers:
        if m.atom_id not in affected: continue
        if template=="OXIDATIVE_ADDITION_2E": m.oxidation_state=counts[m.atom_id]['oxidation_state']+2
        elif template=="REDUCTIVE_ELIMINATION_2E": m.oxidation_state=counts[m.atom_id]['oxidation_state']-2
        m.d_electrons=None;m.valence_electrons=None
        # Spatial consequences are not silently inferred from graph edits.
        m.geometry="unknown";m.relative_positions=[];m.vacant_site_ids=None
    ledger=validate_view(output)
    if identity(source)!=identity(output) or electron_count(source)!=electron_count(output) or sum(a.formal_charge for a in source.atoms)!=ledger["total_charge"]:
        raise FormalError("Metal event lost identity, H, charge or electrons")
    return output,{"event_id":event.event_id,"event_type":event.event_type,"metal_participation":True,
        "metal_template":template,"rule_version":"umr-metal/0.3","atomic_commit":True,
        "formal_status":"passed","chemical_feasibility":"unknown","geometry":geometry,
        "mechanism_evidence":event.mechanism_evidence,"hydrogen_transfer_class":"none",
        "compiled_actions":[a.model_dump(mode="json") for a in actions],"ledger":ledger,
        "uncertainties":["Declared/template-assumed path is not established by its net graph change.",
                         "No kinetic, orbital, solvent or 3D feasibility guarantee."]}
