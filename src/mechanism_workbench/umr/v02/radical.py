from ..chemistry import FormalError, Unsupported, atoms_by_id, order
from .electrons import atom_to_bond, bond_to_atom, bond_to_bond

RADICAL_REGISTRY = {
    "HOMOLYSIS": {"roles":["a","b"],"scope":"ordinary neutral single bond"},
    "RADICAL_ADD_PI": {"roles":["radical","attacked_atom","distal_atom"],"scope":"carbon alkene, explicitly chosen end"},
    "HAT": {"roles":["radical","hydrogen","donor"],"scope":"explicit H/D/T atom transfer"},
    "XAT": {"roles":["radical","halogen","donor"],"scope":"Cl/Br/I atom transfer"},
    "RADICAL_RECOMBINATION": {"roles":["a","b"],"scope":"two local radicals, spin checked separately"},
    "RADICAL_BETA_SCISSION": {"roles":["radical","adjacent","departing"],"scope":"C/O radical beta-scission"},
    "RADICAL_SUBSTITUTION_SH2": {"roles":["radical","sulfur","departing_sulfur"],"scope":"disulfide sulfur only"},
}
RADICAL_DEFERRED = ["DISPROPORTIONATION", "SET", "PCET", "DECARBOXYLATION", "ALKYNE_ADD", "GENERIC_CARBON_SH2"]


def compile_radical(source, event):
    if event.template_id not in RADICAL_REGISTRY:
        raise Unsupported(f"Unregistered radical template: {event.template_id}")
    if event.template_version != "0.2" or event.direction != "forward":
        raise Unsupported("Radical template version/direction not implemented")
    names = RADICAL_REGISTRY[event.template_id]["roles"]
    if set(event.bindings) != set(names) or any(type(event.bindings[n]) is not int for n in names):
        raise FormalError(f"Explicit scalar radical bindings required: {names}")
    ids = [event.bindings[n] for n in names]
    atoms = atoms_by_id(source)
    if len(set(ids)) != len(ids) or any(i not in atoms for i in ids):
        raise FormalError("Missing/repeated radical template atom bindings")
    if event.region_atom_ids and set(event.region_atom_ids) != set(ids):
        raise FormalError("Radical region differs from template bindings")
    if not event.region_atom_ids:
        event.region_atom_ids=list(ids)
    def bond(a,b,value=1):
        if order(source,a,b) != value:
            raise FormalError(f"Radical template requires bond {a}-{b} order {value}")
    def radical(i):
        if atoms[i].radical_electrons != 1:
            raise FormalError("Template requires the explicitly bound occupied radical centre")
    t = event.template_id
    if t == "HOMOLYSIS":
        a,b = ids
        if order(source,a,b)>1:
            raise Unsupported("Multiple-bond homolysis is not implemented")
        bond(a,b)
        if any(atoms[i].formal_charge or atoms[i].radical_electrons for i in ids):
            raise Unsupported("HOMOLYSIS profile requires neutral closed-shell endpoints")
        actions = [bond_to_atom(a,b,"sigma",0,a),bond_to_atom(a,b,"sigma",1,b)]
    elif t == "RADICAL_RECOMBINATION":
        a,b = ids
        radical(a); radical(b); bond(a,b,0)
        actions = [atom_to_bond(a,a,b),atom_to_bond(b,a,b)]
    elif t == "RADICAL_ADD_PI":
        r,a,b = ids
        if order(source,a,b)==3:
            raise Unsupported("Alkyne radical addition requires an independent profile")
        radical(r); bond(a,b,2); bond(r,a,0)
        if any(atoms[i].element != "C" or atoms[i].radical_electrons for i in [a,b]):
            raise Unsupported("RADICAL_ADD_PI accepts only a closed-shell carbon alkene")
        actions = [atom_to_bond(r,r,a),bond_to_bond(a,b,"pi1",0,r,a),bond_to_atom(a,b,"pi1",1,b)]
    elif t in {"HAT","XAT","RADICAL_SUBSTITUTION_SH2"}:
        r,x,d = ids
        radical(r); bond(d,x); bond(r,x,0)
        if t == "HAT" and atoms[x].element != "H":
            raise FormalError("HAT needs an explicitly mapped H/D/T atom")
        if t == "XAT" and atoms[x].element not in {"Cl","Br","I"}:
            raise Unsupported("XAT halogen profile is Cl/Br/I only")
        if t == "RADICAL_SUBSTITUTION_SH2" and any(atoms[i].element != "S" for i in (x,d)):
            raise Unsupported("SH2 profile is restricted to a disulfide S-S centre")
        actions = [atom_to_bond(r,r,x),bond_to_bond(d,x,"sigma",0,r,x),bond_to_atom(d,x,"sigma",1,d)]
    else:
        r,a,b = ids
        radical(r); bond(r,a); bond(a,b)
        if atoms[r].element not in {"C","O"} or atoms[a].element != "C":
            raise Unsupported("Beta scission supports C/O-centred radicals adjacent to carbon")
        actions = [atom_to_bond(r,r,a,"pi1"),bond_to_bond(a,b,"sigma",0,r,a,"pi1"),bond_to_atom(a,b,"sigma",1,b)]
    return actions, {"formal_execution":"passed","template_scope":"matched","orbital_symmetry":"not_applicable",
                     "geometry_status":"unverified","stereochemistry":"unknown",
                     "mechanism_evidence":event.mechanism_evidence,"kinetics_selectivity":"unknown",
                     "radical_polarity":event.radical_polarity,"radical_polarity_is_formal_charge":False}
