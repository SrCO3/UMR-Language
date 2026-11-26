from ..ast import Pair, Shift, Split
from ..chemistry import FormalError, Unsupported, atoms_by_id, order

PERICYCLIC_REGISTRY = {
    "DA_4_2": {"roles":{"diene":4,"dienophile":2},"electrons":6,"profile":"ordinary carbon diene/alkene"},
    "ENE": {"roles":{"ene":4,"enophile":2},"electrons":6,"profile":"carbon ene/alkene, explicit H"},
    "COPE_3_3": {"roles":{"ordered_chain":6},"electrons":6,"profile":"hydrocarbon 1,5-diene"},
    "CLAISEN_3_3": {"roles":{"ordered_chain":6},"electrons":6,"profile":"acyclic aliphatic allyl vinyl ether"},
    "H_SHIFT_1_5": {"roles":{"ordered_chain":5,"hydrogen":1},"electrons":6,"profile":"explicit H pentadienyl"},
    "H_SHIFT_1_3": {"roles":{"ordered_chain":3,"hydrogen":1},"electrons":4,"profile":"explicit H allylic, condition-limited"},
    "H_SHIFT_1_7": {"roles":{"ordered_chain":7,"hydrogen":1},"electrons":8,"profile":"explicit H heptatrienyl"},
    "ELECTROCYCLIC_4PI": {"roles":{"ordered_chain":4},"electrons":4,"profile":"carbon diene closure only"},
    "ELECTROCYCLIC_6PI": {"roles":{"ordered_chain":6},"electrons":6,"profile":"carbon triene closure only"},
    "DIPOLAR_3_2": {"roles":{"dipole":3,"dipolarophile":2},"electrons":6,"profile":"azomethine ylide C(-)-N(+)=C only"},
    "CYCLOADDITION_2_2": {"roles":{"alkene_a":2,"alkene_b":2},"electrons":4,"profile":"ordinary carbon alkenes, photo supra/supra"},
}
PERICYCLIC_DEFERRED = ["RETRO_DA_4_2","RETRO_ENE","ELECTROCYCLIC_4PI_OPEN","ELECTROCYCLIC_6PI_OPEN",
                       "SIGMATROPIC_2_3","CHELETROPIC","HETERO_DA","ALKYNE_DA","OXY_COPE","AROMATIC_CLAISEN",
                       "NITRILE_OXIDE_3_2","NITRONE_3_2","AZIDE_3_2","KETENE_2_2"]


class ConstraintError(FormalError):
    def __init__(self, message, report):
        super().__init__(message)
        self.report = {**report,"formal_execution":"failed"}


def _in_cycle(view,a,b):
    graph = {i.atom_id:set() for i in view.atoms}
    for bond in view.bonds:
        if {bond.a,bond.b} == {a,b}:
            continue
        graph[bond.a].add(bond.b); graph[bond.b].add(bond.a)
    seen, todo = {a}, [a]
    while todo:
        for neighbor in graph[todo.pop()]-seen:
            seen.add(neighbor); todo.append(neighbor)
    return b in seen


def compile_pericyclic(source,event):
    t = event.template_id
    if t not in PERICYCLIC_REGISTRY:
        raise Unsupported(f"Unregistered pericyclic profile: {t}")
    if event.template_version != "0.2" or event.direction != "forward":
        raise Unsupported("Reverse/opening templates and other template versions are not implemented")
    profile = PERICYCLIC_REGISTRY[t]
    if set(event.bindings) != set(profile["roles"]):
        raise FormalError(f"Ordered pericyclic bindings required: {profile['roles']}")
    groups = {}
    for key,length in profile["roles"].items():
        value = event.bindings[key]
        if length == 1:
            if type(value) is not int:
                raise FormalError("Hydrogen role must be one explicit atom_id")
            value = [value]
        if not isinstance(value,list) or len(value) != length:
            raise FormalError(f"Role {key} requires {length} explicitly ordered atoms")
        groups[key] = value
    ids = [i for group in groups.values() for i in group]
    atoms = atoms_by_id(source)
    if len(ids) != len(set(ids)) or any(i not in atoms for i in ids):
        raise FormalError("Missing/repeated pericyclic role atom")
    if len(event.region_atom_ids) != len(set(event.region_atom_ids)) or set(event.region_atom_ids) != set(ids):
        raise FormalError("region_atom_ids must exactly match all ordered template participants")
    if any(atoms[i].radical_electrons for i in ids):
        raise Unsupported("These pericyclic profiles require closed-shell participants")
    if event.mechanism_assumption not in {None,"concerted"}:
        raise FormalError("A declared stepwise mechanism must use explicit stepwise events")
    event.mechanism_assumption = "concerted"
    if event.mechanism_evidence == "unknown":
        event.mechanism_evidence = "assumed_by_template"
    activation = "photochemical" if event.activation == "hν" else event.activation
    report = {"formal_execution":"passed","template_scope":"matched","orbital_symmetry":"unknown",
              "geometry_status":"unverified","stereochemistry":"relative_constraints_only",
              "mechanism_evidence":event.mechanism_evidence,"mechanism_assumption":"concerted",
              "kinetics_selectivity":"unknown","participating_electrons":profile["electrons"],
              "relative_constraints":[],"absolute_configuration":"unknown"}
    def bond(a,b,n):
        if order(source,a,b) != n:
            raise FormalError(f"{t} requires mapped bond {a}-{b} order {n}")
    def element(i,symbol):
        if atoms[i].element != symbol:
            raise Unsupported(f"{t} profile requires {symbol} at mapped atom {i}")
    for i in ids:
        expected = "H" if i in groups.get("hydrogen",[]) or (t=="ENE" and i==groups["ene"][-1]) else "N" if t=="DIPOLAR_3_2" and i==groups["dipole"][1] else "O" if t=="CLAISEN_3_3" and i==groups["ordered_chain"][2] else "C"
        element(i,expected)
        if t != "DIPOLAR_3_2" and atoms[i].formal_charge != 0:
            raise Unsupported("Charged variants require an independently registered profile")
        if t=="DIPOLAR_3_2":
            expected_charge=-1 if i==groups["dipole"][0] else 1 if i==groups["dipole"][1] else 0
            if atoms[i].formal_charge!=expected_charge:
                raise Unsupported("Additional charge states are outside the azomethine ylide profile")
    if t in {"DA_4_2","ENE","COPE_3_3","CLAISEN_3_3","DIPOLAR_3_2"}:
        if activation == "photochemical":
            raise Unsupported("Only the conventional thermal profile is registered for this template")
        expected_face = "supra_supra"
    elif t.startswith("H_SHIFT"):
        length = len(groups["ordered_chain"])
        expected_face = None if activation=="unknown" else "suprafacial" if (length==5 and activation!="photochemical") or (length in {3,7} and activation=="photochemical") else "antarafacial"
    else:
        expected_face = "supra_supra" if t=="CYCLOADDITION_2_2" else None
    if expected_face:
        if event.facial_mode == "unknown":
            event.facial_mode = expected_face
        if event.facial_mode != expected_face:
            if not t.startswith("H_SHIFT"):
                raise Unsupported("This facial configuration is outside the registered profile")
            report["orbital_symmetry"] = "forbidden"
            raise ConstraintError("Facial mode contradicts the supported symmetry configuration",report)
        report["relative_constraints"].append({"kind":"facial","mode":event.facial_mode,"basis":"template"})
    if t=="CYCLOADDITION_2_2" and activation=="thermal":
        report["orbital_symmetry"] = "forbidden"
        raise ConstraintError("Ordinary thermal supra/supra [2+2] is symmetry-forbidden in this profile",report)
    if t.startswith("ELECTROCYCLIC"):
        length = len(groups["ordered_chain"])
        expected_rotation = ("conrotatory" if length==4 else "disrotatory") if activation!="photochemical" else ("disrotatory" if length==4 else "conrotatory")
        if event.rotation_mode == "unknown" and activation!="unknown":
            event.rotation_mode = expected_rotation
        if event.rotation_mode != expected_rotation and activation!="unknown":
            report["orbital_symmetry"] = "forbidden"
            raise ConstraintError("Electrocyclic rotation contradicts the thermal/photo rule",report)
        report["relative_constraints"].append({"kind":"rotation","mode":event.rotation_mode,"basis":"template"})
    if activation != "unknown":
        report["orbital_symmetry"] = "allowed"
    if event.conditions.get("geometry_status") == "incompatible" or (t=="DA_4_2" and event.conditions.get("s_cis") is False):
        report["geometry_status"] = "incompatible"
        raise ConstraintError("Supplied geometry is explicitly incompatible with the template",report)
    if t=="CYCLOADDITION_2_2" and activation=="photochemical":
        event.model_assumption = "textbook_concerted_default"
        report["model_assumption"] = event.model_assumption
        report["execution_note"] = "按教材级默认协同模板执行"
    if t=="DA_4_2":
        d1,d2,d3,d4 = groups["diene"]; e1,e2 = groups["dienophile"]
        for a,b,n in [(d1,d2,2),(d2,d3,1),(d3,d4,2),(e1,e2,2),(d1,e1,0),(d4,e2,0)]: bond(a,b,n)
        actions = [Shift(a=d2,b=d1,c=e1),Shift(a=e1,b=e2,c=d4),Shift(a=d4,b=d3,c=d2)]
        _relative_alkene(source,e1,e2,report)
    elif t=="ENE":
        a1,a2,a3,h = groups["ene"]; b1,b2 = groups["enophile"]
        for a,b,n in [(a1,a2,2),(a2,a3,1),(a3,h,1),(b1,b2,2),(a1,b1,0),(b2,h,0)]: bond(a,b,n)
        actions = [Shift(a=a2,b=a1,c=b1),Shift(a=b1,b=b2,c=h),Shift(a=h,b=a3,c=a2)]
    elif t in {"COPE_3_3","CLAISEN_3_3"}:
        c1,c2,c3,c4,c5,c6 = groups["ordered_chain"]
        for a,b,n in [(c1,c2,2),(c2,c3,1),(c3,c4,1),(c4,c5,1),(c5,c6,2),(c1,c6,0)]: bond(a,b,n)
        if t=="CLAISEN_3_3" and _in_cycle(source,c1,c2):
            raise Unsupported("Aromatic/cyclic-vinyl Claisen profile not implemented")
        if t=="COPE_3_3" and any(atoms[b.b if b.a==i else b.a].element=="O" for i in ids for b in source.bonds if i in {b.a,b.b}):
            raise Unsupported("Oxy-Cope requires its own profile")
        actions = [Shift(a=c2,b=c1,c=c6),Shift(a=c6,b=c5,c=c4),Shift(a=c4,b=c3,c=c2)]
    elif t.startswith("H_SHIFT"):
        chain, h = groups["ordered_chain"],groups["hydrogen"][0]
        bond(h,chain[0],1); bond(h,chain[-1],0)
        for i in range(len(chain)-1): bond(chain[i],chain[i+1],1 if i%2==0 else 2)
        actions = [Shift(a=chain[i+1],b=chain[i],c=chain[i-1]) for i in range(1,len(chain)-1,2)]
        actions.append(Shift(a=chain[0],b=h,c=chain[-1]))
    elif t.startswith("ELECTROCYCLIC"):
        chain = groups["ordered_chain"]
        for i in range(len(chain)-1): bond(chain[i],chain[i+1],2 if i%2==0 else 1)
        bond(chain[0],chain[-1],0)
        actions = [Shift(a=chain[1],b=chain[0],c=chain[-1])]
        actions.extend(Shift(a=chain[i],b=chain[i-1],c=chain[i-2]) for i in range(len(chain)-1,2,-2))
    elif t=="DIPOLAR_3_2":
        c1,n2,c3 = groups["dipole"]; c4,c5 = groups["dipolarophile"]
        for a,b,n in [(c1,n2,1),(n2,c3,2),(c4,c5,2),(c1,c5,0),(c3,c4,0)]: bond(a,b,n)
        if atoms[c1].formal_charge!=-1 or atoms[n2].formal_charge!=1 or atoms[c1].nonbonding_electrons!=2:
            raise Unsupported("DIPOLAR_3_2 implements the explicit azomethine ylide profile only")
        actions = [Pair(a=c1,b=c5),Shift(a=c5,b=c4,c=c3),Split(a=c3,b=n2,target=n2)]
    else:
        a1,a2 = groups["alkene_a"]; b1,b2 = groups["alkene_b"]
        for a,b,n in [(a1,a2,2),(b1,b2,2),(a1,b1,0),(a2,b2,0)]: bond(a,b,n)
        actions = [Shift(a=a2,b=a1,c=b1),Shift(a=b1,b=b2,c=a2)]
        _relative_alkene(source,a1,a2,report); _relative_alkene(source,b1,b2,report)
    for constraint in event.stereo_constraints:
        if constraint.get("kind") not in {"absolute","rotation","facial","relative"}:
            raise Unsupported("Unrecognized stereo constraint; cannot silently discard it")
        if constraint.get("kind")=="absolute":
            raise Unsupported("Absolute product stereochemistry is not assigned by v0.2 templates")
        if constraint.get("kind") in {"rotation","facial"}:
            expected = event.rotation_mode if constraint["kind"]=="rotation" else event.facial_mode
            if constraint.get("mode")!=expected:
                raise ConstraintError("Explicit stereo constraint contradicts template mode",report)
        if constraint.get("kind")=="relative":
            matches = [c for c in report["relative_constraints"] if c.get("kind")=="relative" and set(c["atoms"])==set(constraint.get("atoms",[]))]
            if not matches:
                raise Unsupported("This relative product assertion cannot be inferred by the implemented constraint solver")
            if any(c["relation"]!=constraint.get("relation") for c in matches):
                raise ConstraintError("Product relative configuration contradicts the supplied alkene stereo",report)
    event.stereo_constraints = event.stereo_constraints or list(report["relative_constraints"])
    return actions,report


def _relative_alkene(source,a,b,report):
    bond = next(bond for bond in source.bonds if {bond.a,bond.b}=={a,b})
    stereo = bond.stereo
    if stereo.get("status")=="specified" and stereo.get("rdkit_stereo") in {"STEREOE","STEREOZ","STEREOCIS","STEREOTRANS"}:
        report["relative_constraints"].append({"kind":"relative","atoms":stereo["reference_atoms"],
            "relation":"trans" if stereo["rdkit_stereo"] in {"STEREOE","STEREOTRANS"} else "cis",
            "basis":"specified_alkene_stereo_and_supra_template","assignment":"constraint_only"})
