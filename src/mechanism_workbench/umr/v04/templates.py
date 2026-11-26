"""Registered finite reactive-center templates, compiled to existing PAIR/SHIFT/SPLIT/FISHHOOK."""
from ..ast import Pair, Shift, Split
from ..chemistry import FormalError, Unsupported, atoms_by_id, bond_valences, order
from ..v02.electrons import atom_to_bond, bond_to_atom, bond_to_bond
from .ledger import check_isolated_ring, occupation_profiles

RULE_VERSION="umr-reactive/0.4"
SOURCES={
    "hine1950":"https://doi.org/10.1021/ja01162a024",
    "methylene_ethene":"https://d-nb.info/1102811173/34",
    "methylene_insertion":"https://doi.org/10.1021/ja00331a007",
    "nitrene_ethene":"https://www.sciencedirect.com/science/article/abs/pii/S0166128007002473",
    "nitrene_small_molecules":"https://pmc.ncbi.nlm.nih.gov/articles/PMC12728564/",
    "benzyne1953":"https://doi.org/10.1021/ja01109a523",
    "CH_spin_caveat":"https://goldbook.iupac.org/terms/view/C00854",
}

TEMPLATES={
    "CARBENE_ALPHA_ELIMINATION":{"event_type":"POLAR","roles":["center","leaving_group"],"scope":"CCl3- -> paired CCl2 + Cl-; explicit trichloromethyl anion", "source":"hine1950"},
    "CARBENE_CYCLOPROPANATION":{"event_type":"POLAR","roles":["center","alkene_a","alkene_b"],"scope":"free paired carbene + isolated neutral C=C; one concerted [2+1] event", "source":"methylene_ethene"},
    "NITRENE_AZIRIDINATION":{"event_type":"POLAR","roles":["center","alkene_a","alkene_b"],"scope":"free paired nitrene + isolated neutral C=C; assumed concerted event", "source":"nitrene_ethene"},
    "CARBENE_CH_INSERTION":{"event_type":"POLAR","roles":["center","donor","hydrogen"],"scope":"paired free carbene into one explicit saturated neutral C-H bond", "source":"methylene_insertion"},
    "NITRENE_CH_INSERTION":{"event_type":"POLAR","roles":["center","donor","hydrogen"],"scope":"paired free nitrene into one explicit saturated neutral C-H bond; formal local model", "source":"nitrene_small_molecules"},
    "TRIPLET_CARBENE_FIRST_ADDITION":{"event_type":"RADICAL","roles":["center","alkene_a","alkene_b"],"scope":"two singly occupied carbene slots; only first C-C formation, retain 1,3-diradical", "source":"methylene_ethene"},
    "TRIPLET_NITRENE_FIRST_ADDITION":{"event_type":"RADICAL","roles":["center","alkene_a","alkene_b"],"scope":"two singly occupied nitrene slots; only first N-C formation, retain 1,3-diradical", "source":"nitrene_ethene"},
    "CARBYNE_DOUBLET_HAT":{"event_type":"RADICAL","roles":["center","donor","hydrogen"],"scope":"one-unpaired free carbyne + saturated C-H -> paired carbene + C radical; formal ledger example, not asserted dominant pathway", "source":"CH_spin_caveat"},
    "ARYNE_DEPROTONATION":{"event_type":"POLAR","roles":["ring","alpha_carbon","beta_carbon","hydrogen","leaving_group","base"],"scope":"six-carbon ortho-H chlorobenzene Kekule region + NH2-", "source":"benzyne1953"},
    "ARYNE_ELIMINATION":{"event_type":"POLAR","roles":["ring","alpha_carbon","beta_carbon","leaving_group"],"scope":"explicit ortho-chloro aryl carbanion -> formal-triple 1,2-aryne + Cl-", "source":"benzyne1953"},
    "ARYNE_NUCLEOPHILIC_ADDITION":{"event_type":"POLAR","roles":["aryne_a","aryne_b","nucleophile"],"scope":"declared formal-triple 1,2-aryne + NH2- -> mapped aryl carbanion", "source":"benzyne1953"},
    "ARYL_CARBA_PROTONATION":{"event_type":"POLAR","roles":["carbon","donor","hydrogen"],"scope":"explicit aryl carbanion + explicit NH3 H -> aryl C-H + NH2-", "source":"benzyne1953"},
}


def compile_template(view,event,*,_aryne_halides=frozenset({'Cl'}),_allow_spectator_metals=False):
    name=event.reactive_template
    if name not in TEMPLATES: raise Unsupported(f"Unregistered reactive-center template: {name}")
    spec=TEMPLATES[name]; atoms=atoms_by_id(view); val=bond_valences(view); profiles=occupation_profiles(view)
    spectator_only=_allow_spectator_metals and not view.metal_bindings and all(val[m.atom_id]==0 for m in view.metal_centers)
    if (view.metal_centers and not spectator_only) or event.metal_participation or event.metal_template:
        raise Unsupported("Free intermediate templates do not authorize metal carbenes/carbynes/nitrenes")
    if event.template_id:
        raise FormalError("Do not combine unrelated template_id and reactive_template")
    if event.event_type!=spec["event_type"]: raise FormalError("Event classification disagrees with registered reactive template")
    if not event.mechanism_assumption or not event.evidence:
        raise FormalError("Reactive template needs an explicit mechanism assumption and evidence provenance")
    if set(event.bindings)!=set(spec["roles"]): raise FormalError("Template requires exactly its registered atom roles")
    scalars=[]
    for role in spec["roles"]:
        value=event.bindings[role]
        if role=="ring":
            if not isinstance(value,list) or len(value)!=6 or len(set(value))!=6 or any(i not in atoms for i in value):
                raise FormalError("ring role needs six distinct mapped atoms")
        elif type(value) is not int or value not in atoms: raise FormalError("Missing/non-scalar role atom")
        else: scalars.append(value)
    if len(scalars)!=len(set(scalars)): raise FormalError("Template scalar roles must be distinct")
    region=set(scalars)|set(event.bindings.get("ring",[]))
    if event.region_atom_ids and set(event.region_atom_ids)!=region:
        raise FormalError("Reactive template region must match bound atom roles")
    event.region_atom_ids=sorted(region)
    b=event.bindings; created={}; aryne=None
    def require(condition,message):
        if not condition: raise FormalError(message)
    def neutral_carbon(i):
        return atoms[i].element=="C" and atoms[i].formal_charge==0 and atoms[i].nonbonding_electrons==0 and atoms[i].radical_electrons==0 and val[i]==4
    def amide(i):
        require(atoms[i].element=="N" and atoms[i].formal_charge==-1 and atoms[i].nonbonding_electrons==4 and atoms[i].radical_electrons==0 and val[i]==2,"Base/nucleophile must be NH2- in this finite template")
        require(atoms[i].hydrogens+sum(atoms[j].element=="H" for j in atoms if order(view,i,j))==2,"Only NH2- is registered")
    def carbon_h():
        x,h=b["donor"],b["hydrogen"]
        require(atoms[h].element=="H" and order(view,x,h)==1 and neutral_carbon(x),"Insertion/HAT needs designated explicit H on neutral saturated carbon")
        require(all(bond.order==1 for bond in view.bonds if x in {bond.a,bond.b}),"Donor must be saturated sp3 C-H")
        return x,h
    def center(profile):
        i=b["center"]
        require(profiles.get(i)==profile,"Input center occupation differs from required profile")
        if len(view.reactive_centers)!=1: raise Unsupported("Multiple special centers require separately registered coupling/template scope")
        if profile.endswith("paired") and view.electronic_state.character!="closed_shell":
            raise Unsupported("Paired-center concerted template requires declared closed-shell electronic state")
        if "two_singly" in profile and (view.electronic_state.character!="open_shell" or view.electronic_state.coupling!="high_spin"):
            raise Unsupported("Triplet template requires explicit high-spin correlation; no inference from slot count")
        return i
    if name=="CARBENE_ALPHA_ELIMINATION":
        c,x=b["center"],b["leaving_group"]
        ns=[i for i in atoms if order(view,c,i)]
        require(atoms[c].element=="C" and (val[c],atoms[c].nonbonding_electrons,atoms[c].radical_electrons,atoms[c].formal_charge)==(3,2,0,-1),"Alpha elimination requires explicit CCl3- precursor")
        require(len(ns)==3 and x in ns and all(atoms[i].element=="Cl" and order(view,c,i)==1 for i in ns),"Only trichloromethyl anion alpha elimination registered")
        actions=[Split(a=c,b=x,target=x)]; created[c]="carbene_paired"
    elif name in {"CARBENE_CYCLOPROPANATION","NITRENE_AZIRIDINATION","TRIPLET_CARBENE_FIRST_ADDITION","TRIPLET_NITRENE_FIRST_ADDITION"}:
        species="carbene" if "CARBENE" in name else "nitrene"; triplet=name.startswith("TRIPLET")
        c=center(species+("_two_singly" if triplet else "_paired")); a,z=b["alkene_a"],b["alkene_b"]
        require(order(view,a,z)==2 and neutral_carbon(a) and neutral_carbon(z),"Registered addition requires neutral tetravalent C=C")
        require(not order(view,c,a) and not order(view,c,z),"Center must not already bond to alkene")
        require(not any(x.order>1 and {x.a,x.b}!={a,z} and {x.a,x.b}&{a,z} for x in view.bonds),"Cumulated alkene scope unsupported")
        if triplet:
            actions=[atom_to_bond(c,c,a),bond_to_bond(a,z,"pi1",0,c,a),bond_to_atom(a,z,"pi1",1,z)]
        else: actions=[Pair(a=c,b=a),Shift(a=a,b=z,c=c)]
    elif name in {"CARBENE_CH_INSERTION","NITRENE_CH_INSERTION","CARBYNE_DOUBLET_HAT"}:
        x,h=carbon_h()
        if name=="CARBYNE_DOUBLET_HAT":
            c=center("carbyne_one_singly")
            if view.electronic_state.coupling!="one_unpaired": raise Unsupported("Carbyne HAT requires declared one-unpaired doublet model")
            actions=[atom_to_bond(c,c,h),bond_to_bond(x,h,"sigma",0,c,h),bond_to_atom(x,h,"sigma",1,x)]
            created[c]="carbene_paired"
        else:
            c=center("carbene_paired" if name.startswith("CARBENE") else "nitrene_paired")
            actions=[Pair(a=c,b=h),Shift(a=h,b=x,c=c)]
        require(not order(view,c,x) and not order(view,c,h),"Insertion/HAT sites must initially be disconnected")
    elif name in {"ARYNE_DEPROTONATION","ARYNE_ELIMINATION"}:
        a,z,x=b["alpha_carbon"],b["beta_carbon"],b["leaving_group"]; ring=b["ring"]
        check_isolated_ring(view,ring)
        require(a in ring and z in ring and all(atoms[i].element=="C" for i in ring),"Six-carbon ring roles required")
        require(order(view,a,z)==2 and atoms[x].element in _aryne_halides and order(view,a,x)==1,"Ortho elimination requires a registered aryl halide on the selected Kekule double-bond pair")
        orders=[order(view,i,ring[(j+1)%6]) for j,i in enumerate(ring)]
        require(orders in ([2,1,2,1,2,1],[1,2,1,2,1,2]),"Only alternating six-carbon aromatic ring registered")
        if name=="ARYNE_DEPROTONATION":
            h,q=b["hydrogen"],b["base"]; amide(q)
            require(atoms[h].element=="H" and order(view,z,h)==1 and neutral_carbon(z),"Ortho H must be explicit and attached to designated beta carbon")
            actions=[Pair(a=q,b=h),Split(a=h,b=z,target=z)]
        else:
            require((val[z],atoms[z].nonbonding_electrons,atoms[z].radical_electrons,atoms[z].formal_charge)==(3,2,0,-1),"Elimination requires a distinct explicit ortho aryl anion")
            actions=[Pair(a=z,b=a),Split(a=a,b=x,target=x)]
            aryne={"region_id":event.event_id+".aryne","ring_atom_ids":ring,"aryne_atom_ids":[a,z]}
    elif name=="ARYNE_NUCLEOPHILIC_ADDITION":
        a,z,q=b["aryne_a"],b["aryne_b"],b["nucleophile"]; amide(q)
        require(any(set(r.aryne_atom_ids)=={a,z} and r.representation=="ortho_formal_triple" for r in view.aryne_regions),"Addition requires declared and validated 1,2-aryne")
        require(not order(view,q,a),"Nucleophile already attached")
        actions=[Pair(a=q,b=a),Split(a=a,b=z,target=z)]
    else:
        c,n,h=b["carbon"],b["donor"],b["hydrogen"]
        require(atoms[c].element=="C" and (val[c],atoms[c].nonbonding_electrons,atoms[c].radical_electrons,atoms[c].formal_charge)==(3,2,0,-1),"Explicit carbon anion required")
        require(atoms[n].element=="N" and val[n]==3 and atoms[n].nonbonding_electrons==2 and not atoms[n].formal_charge and atoms[h].element=="H" and order(view,n,h)==1,"Proton donor must be explicit NH3 H")
        require(atoms[n].hydrogens+sum(atoms[j].element=="H" for j in atoms if order(view,n,j))==3,"Only NH3 proton donor registered")
        actions=[Pair(a=c,b=h),Split(a=h,b=n,target=n)]
    return actions,created,aryne,{"template_id":name,"rule_version":RULE_VERSION,"scope":spec["scope"],
        "source":SOURCES[spec["source"]],"evidence_status":"formal_example" if name=="CARBYNE_DOUBLET_HAT" else "template_assumption_with_literature_context",
        "geometry_status":"unknown","kinetics_selectivity":"unknown","chemical_feasibility":"unknown"}
