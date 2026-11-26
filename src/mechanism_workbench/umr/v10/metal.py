"""Finite molecular channels for titanium tartrate and palladium allylic epoxide chemistry."""
from collections import Counter

from ..ast import Pair, Shift, Split
from ..chemistry import FormalError, atoms_by_id, order
from ..v09.ast import Ligand, MetalBinding
from ..v09.coordination import checks
from ..v09.execution import transact
from ..v09.ledger import validate_view


SPECS = {
    "TI_ALKOXIDE_PROTON_EXCHANGE": ("ti", "leaving", "entering", "hydrogen"),
    "TI_TARTRATE_PEROXIDE_TRANSFER": ("ti", "alcohol_o", "alcohol_c", "proximal_c",
                                     "distal_c", "chain_c", "transferred_o", "alkyl_o", "tert_c"),
    "PD_VINYL_EPOXIDE_OPENING": ("pd", "first", "middle", "last", "epoxide_o", "retained_c"),
    "PD_ALLYL_AZOLE_CAPTURE": ("pd", "first", "middle", "last", "nitrogen", "alkoxide"),
}


def require(test, message):
    if not test:
        raise FormalError(message)


def neighbours(view, index):
    return {b.b if b.a == index else b.a for b in view.bonds if index in (b.a, b.b)}


def organic_component(view, start):
    metals = {m.atom_id for m in view.metal_centers}
    seen, pending = {start}, [start]
    while pending:
        new = neighbours(view, pending.pop()) - seen - metals
        seen.update(new)
        pending.extend(new)
    return sorted(seen)


def single_ligand(view, metal, donor):
    hits = [b for b in view.metal_bindings if b.metal_atom_id == metal
            and b.donor_atom_ids == [donor] and b.mode == "sigma_covalent"]
    require(len(hits) == 1, "Expected a terminal metal–oxygen ligand")
    ligand_id = hits[0].ligand_id
    require(sum(b.ligand_id == ligand_id for b in view.metal_bindings) == 1,
            "The reacting ligand must be terminal and monodentate")
    return ligand_id


def remove_ligands(view, ids):
    view.ligands = [l for l in view.ligands if l.ligand_id not in ids]
    view.metal_bindings = [b for b in view.metal_bindings if b.ligand_id not in ids]


def add_sigma_x(view, metal, donor, key):
    ligand = Ligand(ligand_id=key, atom_ids=organic_component(view, donor), donor_atom_ids=[donor],
                    ionic_class="X", ionic_charge=-1, assignment_basis="terminal covalent X ligand")
    view.ligands.append(ligand)
    view.metal_bindings.append(MetalBinding(binding_id=key+".sigma", ligand_id=key,
        metal_atom_id=metal, donor_atom_ids=[donor], mode="sigma_covalent",
        ionic_allocation=-1, donated_electrons=2))


def titanium(view, event, ids, counts):
    atoms = atoms_by_id(view)
    ti = ids[0]
    require(ti in counts and atoms[ti].element == "Ti" and counts[ti]["oxidation_state"] == 4,
            "Titanium(IV) required at the reacting site")
    # Each tartrate chelates one Ti with two alkoxides and bridges the other.
    det = [l for l in view.ligands if l.ionic_class == "L2X2" and len(l.donor_atom_ids) >= 3
           and all(atoms[i].element == "O" for i in l.donor_atom_ids)
           and {atoms[i].element for i in l.atom_ids} <= {"C", "O", "H"}
           and sum(atoms[i].element == "O" for i in l.atom_ids) == 6
           and sum(atoms[i].element == "C" for i in l.atom_ids) == 8]
    partners = [{b.metal_atom_id for b in view.metal_bindings if b.ligand_id == l.ligand_id} for l in det]
    require(len(partners) == 2 and partners[0] == partners[1] and ti in partners[0]
            and len(partners[0]) == 2 and all(counts[i]["element"] == "Ti"
            and counts[i]["oxidation_state"] == 4 for i in partners[0]),
            "The declared binuclear Ti(IV) tartrate environment is required")
    if event.general_event == "TI_ALKOXIDE_PROTON_EXCHANGE":
        _, old, new, h = ids
        require([atoms[i].element for i in (old, new, h)] == ["O", "O", "H"]
                and order(view, ti, old) == order(view, new, h) == 1
                and not order(view, ti, new) and not order(view, old, h),
                "Exchange requires a terminal Ti–O ligand and an explicit incoming O–H")
        leaving = single_ligand(view, ti, old)
        actions = [Shift(a=ti, b=old, c=h), Shift(a=h, b=new, c=ti)]
        def partition(target):
            remove_ligands(target, {leaving})
            add_sigma_x(target, ti, new, event.event_id+".alkoxide")
        return actions, partition, set(), {}
    _, ox, carbinol, near, far, chain, oxy, alkoxy, tert = ids
    require([atoms[i].element for i in ids[1:]] == ["O","C","C","C","C","O","O","C"],
            "Ti transfer requires an allylic alkoxide and tert-butyl hydroperoxide")
    for left, right, bond_order in ((ti,ox,1),(ox,carbinol,1),(carbinol,near,1),
            (near,far,2),(far,chain,1),(ti,oxy,1),(oxy,alkoxy,1),(alkoxy,tert,1)):
        require(order(view,left,right) == bond_order,
                f"Incorrect oxygen identity or reacting topology at {left}–{right}")
    require(neighbours(view,oxy) == {ti,alkoxy} and neighbours(view,alkoxy) == {oxy,tert}
            and sum(atoms[i].element == "C" for i in neighbours(view,tert)) == 3
            and atoms[near].hydrogens == atoms[far].hydrogens == 1,
            "Only the terminal peroxide oxygen is delivered in this E-disubstituted substrate")
    alkene = next(b for b in view.bonds if {b.a,b.b} == {near,far})
    require(alkene.stereo.get("rdkit_stereo") in {"STEREOE","STEREOTRANS"}
            and set(alkene.stereo.get("reference_atoms",[])) == {carbinol,chain},
            "The two declared carbon substituents must define an E alkene")
    for ligand in det:
        centers = [i for i in ligand.atom_ids if atoms[i].element=="C" and atoms[i].hydrogens==1]
        require(len(centers)==2 and all(atoms[i].stereo.get("status")=="specified" for i in centers),
                "Both tartrate stereocenters must be specified")
    face = event.conditions.get("approach_face")
    require(face in {"toward","away"}, "Specify the approach face in the substrate reference frame")
    old = {single_ligand(view,ti,ox), single_ligand(view,ti,oxy)}
    actions = [Pair(a=oxy,b=near), Shift(a=near,b=far,c=oxy),
               Split(a=oxy,b=alkoxy,target=alkoxy), Split(a=ti,b=oxy,target=oxy), Pair(a=alkoxy,b=ti)]
    def partition(target):
        remove_ligands(target, old)
        add_sigma_x(target,ti,ox,event.event_id+".epoxyalcohol")
        add_sigma_x(target,ti,alkoxy,event.event_id+".tert_butoxide")
        tag = "CHI_TETRAHEDRAL_CCW" if face == "toward" else "CHI_TETRAHEDRAL_CW"
        for center, refs in ((near,[oxy,far,carbinol]),(far,[oxy,near,chain])):
            atoms_by_id(target)[center].stereo = dict(status="specified",rdkit_tag=tag,
                neighbor_order=refs,basis="declared syn approach to the E substrate")
    return actions, partition, set(), dict(transferred_oxygen=oxy, retained_peroxide_oxygen=alkoxy,
            approach_face=face, frame="proximal alkene carbon left, distal right, CH2OH above",
            face_source="declared model; the checker evaluates the stereochemical consequence")


def palladium(view, event, ids, counts):
    atoms = atoms_by_id(view)
    pd, a, b, c, x, y = ids
    require(pd in counts and atoms[pd].element == "Pd" and all(atoms[i].element == "C" for i in (a,b,c))
            and order(view,a,b)==1 and order(view,b,c)==2,
            "Pd and a contiguous allyl fragment are required")
    phosphines = [binding for binding in view.metal_bindings if binding.metal_atom_id == pd
                  and binding.mode == "dative" and all(atoms[i].element=="P" for i in binding.donor_atom_ids)]
    require(len(phosphines)==2, "The active molecular model uses bis(phosphine)palladium")
    if event.general_event == "PD_VINYL_EPOXIDE_OPENING":
        o, retained = x, y
        require(counts[pd]["oxidation_state"]==0 and counts[pd]["d_electrons"]==10
                and atoms[o].element=="O" and atoms[retained].element=="C"
                and all(order(view,i,j)==1 for i,j in ((a,o),(retained,o),(a,retained)))
                and order(view,pd,a)==0,
                "Pd(0) and an intact three-membered allylic epoxide are required")
        actions = [Pair(a=pd,b=a), Split(a=a,b=o,target=o)]
        def partition(target):
            key = event.event_id+".allyl"
            target.ligands.append(Ligand(ligand_id=key, atom_ids=organic_component(target,a),
                donor_atom_ids=[a,b,c], ionic_class="LX", ionic_charge=-2, remote_ionic_charge=-1,
                ligand_type="allyl", assignment_basis="localized sigma + eta2, remote alkoxide"))
            target.metal_bindings.extend([
                MetalBinding(binding_id=key+".sigma",ligand_id=key,metal_atom_id=pd,donor_atom_ids=[a],
                             mode="sigma_covalent",ionic_allocation=-1,donated_electrons=2),
                MetalBinding(binding_id=key+".pi",ligand_id=key,metal_atom_id=pd,donor_atom_ids=[b,c],
                             mode="eta",hapticity=2,ionic_allocation=0,donated_electrons=2,pi_pool_slot=0)])
        return actions, partition, set(), dict(haptic_region=[a,b,c],remote_alkoxide=o)
    n,o = x,y
    require(counts[pd]["oxidation_state"]==2 and atoms[n].element=="N" and atoms[n].formal_charge==0
            and atoms[n].nonbonding_electrons==2 and atoms[n].hydrogens==0
            and atoms[o].element=="O" and atoms[o].formal_charge==-1 and order(view,pd,a)==1,
            "Pd(II) allyl–alkoxide and a neutral imidazole nitrogen are required")
    # Verify the five-membered 1,3-diazole ring and its pyridine-like donor.
    ring = organic_component(view,n)
    heavy = [i for i in ring if atoms[i].element!="H"]
    require(len(heavy)==5 and Counter(atoms[i].element for i in heavy)=={"C":3,"N":2}
            and all(len(neighbours(view,i)&set(heavy))==2 for i in heavy)
            and all(atoms[i].element=="C" for i in neighbours(view,n))
            and any(order(view,n,i)==2 for i in neighbours(view,n)), "Expected free imidazole")
    hits = [l for l in view.ligands if l.ligand_type=="allyl" and set(l.donor_atom_ids)=={a,b,c}
            and o in l.atom_ids and l.remote_ionic_charge==-1]
    require(len(hits)==1, "Missing allyl–alkoxide ligand")
    bound = [k for k in view.metal_bindings if k.ligand_id==hits[0].ligand_id and k.metal_atom_id==pd]
    require(any(k.mode=="eta" and set(k.donor_atom_ids)=={b,c} for k in bound),
            "The adjacent pi bond must belong to the Pd allyl ligand")
    actions = [Pair(a=n,b=c), Shift(a=c,b=b,c=a), Split(a=a,b=pd,target=pd)]
    def partition(target):
        remove_ligands(target,{hits[0].ligand_id})
    return actions, partition, {k.binding_id for k in bound}, dict(attack_atom=c,nucleophile=n)


def execute(view, event):
    name = event.general_event
    require(event.event_type=="POLAR" and event.execution_level=="electron_actions" and event.evidence
            and event.mechanism_assumption and not event.exception_ids, "An evidenced polar channel is required")
    require(not any((event.template_id,event.rule_template,event.metal_template,event.element_template,
                     event.reactive_template,event.net_template,event.electron_transfers,event.electronic_outcome)),
            "One finite channel per event")
    roles = event.general_roles
    require(set(roles)==set(SPECS[name]) and len(set(roles.values()))==len(roles)
            and set(roles.values())<=atoms_by_id(view).keys(), "Exact and distinct mapped roles are required")
    ids = [roles[key] for key in SPECS[name]]
    validate_view(view)
    counts = {row["metal_atom_id"]:row for row in checks(view)[1]}
    compiler = titanium if name.startswith("TI_") else palladium
    actions, after, released, details = compiler(view,event,ids,counts)
    require(event.actions is None or Counter(a.model_dump_json() for a in event.actions)==
            Counter(a.model_dump_json() for a in actions), "Declared arrows disagree with this channel")
    target = transact(view,actions,event.product_state_id,event.product_view_id,after=after,
                      _released_binding_ids=released,
                      _declared_pair_donors={ids[0]} if name=="PD_VINYL_EPOXIDE_OPENING" else set())
    event.actions = actions
    return target, dict(event_id=event.event_id,event_type=event.event_type,general_event=name,
        rule_version="umr-rules/1.0",formal_status="passed",atomic_commit=True,
        execution_mode="electron_actions",compiled_actions=[a.model_dump(mode="json") for a in actions],
        ledger=validate_view(target), chemical_feasibility="unknown", **details)
