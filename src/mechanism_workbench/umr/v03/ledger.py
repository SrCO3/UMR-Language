"""Covalent ownership ledger and separate ionic-method coordination analysis.

Dative/eta pairs remain owned by the ligand in N; reserved pairs cannot be spent
again. They are counted ONCE globally, and independently in local metal VE analysis.
This is an explicit representation convention, not a claim about electron density.
"""
import json
from collections import Counter

from ..ast import Atom, Bond, LewisView as View01
from ..chemistry import (VALENCE, FormalError, Unsupported, atoms_by_id, bond_orders, bond_valences,
                         compare_views as compare01, electron_count, graph_signature as signature01)
from ..v02.ledger import validate_view as validate02

GROUP = {"Pd": 10, "Zn": 12}
# Other metals are preserved by import, but not silently assigned a supported ledger.
METALS = set("Li Na K Rb Cs Be Mg Ca Sr Ba Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg Al Ga In Tl Sn Pb La Ce Nd Sm U".split())


def charge(atom, valence):
    group = GROUP.get(atom.element, VALENCE.get(atom.element))
    if group is None or atom.nonbonding_electrons is None:
        raise Unsupported(f"Unsupported electron ownership at {atom.atom_id}:{atom.element}")
    return group - atom.nonbonding_electrons - valence


def metal_signature(view):
    return json.dumps({"centers": sorted((x.model_dump(mode="json") for x in view.metal_centers), key=lambda x:x["atom_id"]),
        "ligands": sorted((x.model_dump(mode="json") for x in view.ligands), key=lambda x:x["ligand_id"]),
        "bindings": sorted((x.model_dump(mode="json") for x in view.metal_bindings), key=lambda x:x["binding_id"]),
        "convention":view.electron_ownership_convention}, sort_keys=True)


def signature(view):
    return signature01(view), metal_signature(view)


def compare_views(actual, expected):
    compare01(actual, expected)
    if metal_signature(actual) != metal_signature(expected):
        raise FormalError("Metal coordination/oxidation state/ligand identity or geometry differs from asserted target")


def binding_checks(view):
    atoms, orders = atoms_by_id(view), bond_orders(view)
    centers = {m.atom_id:m for m in view.metal_centers}
    metal_ids = {i for i,a in atoms.items() if a.element in METALS-{"Li","Na","K"}}
    if set(centers) != metal_ids:
        raise FormalError("Every metal atom requires exactly one explicit MetalCenter")
    if any(atoms[i].element not in GROUP for i in metal_ids):
        raise Unsupported("Metal ledger supports Pd and Zn only; other metals retained but unsupported")
    ligands = {x.ligand_id:x for x in view.ligands}
    membership = set()
    for ligand in view.ligands:
        members, donors = set(ligand.atom_ids), set(ligand.donor_atom_ids)
        if len(members)!=len(ligand.atom_ids) or len(donors)!=len(ligand.donor_atom_ids):
            raise FormalError("Duplicate ligand atom/donor references")
        if not members<=atoms.keys() or members & metal_ids or not donors<=members:
            raise FormalError("Ligand membership/donor references are invalid")
        if membership & members:
            raise Unsupported("Overlapping/bridging ligand definitions are outside this adapter")
        membership |= members
        seen={next(iter(members))}
        while True:
            expanded=seen|{j for i in seen for j in members if tuple(sorted((i,j))) in orders}
            if expanded==seen: break
            seen=expanded
        if seen!=members:
            raise FormalError("A ligand must contain a connected complete organic fragment")
        if ligand.ionic_class!="unknown" and ligand.ionic_charge!={"L":0,"X":-1}[ligand.ionic_class]:
            raise FormalError("Ionic L/X class and declared ligand charge disagree")
    used_ligands, sites, covered, reserved = set(), set(), set(), Counter()
    for binding in view.metal_bindings:
        if binding.metal_atom_id not in centers or binding.ligand_id not in ligands:
            raise FormalError("Binding refers to an absent metal or ligand")
        ligand=ligands[binding.ligand_id]
        if binding.ligand_id in used_ligands:
            raise Unsupported("Chelation/bridging/multiple binding modes need a future adapter")
        used_ligands.add(binding.ligand_id)
        donors=set(binding.donor_atom_ids)
        if len(donors)!=len(binding.donor_atom_ids) or donors!=set(ligand.donor_atom_ids):
            raise FormalError("Binding donors must match the declared ligand donors")
        if binding.site_id:
            key=(binding.metal_atom_id,binding.site_id)
            if key in sites: raise FormalError("Two bindings occupy the same metal site")
            sites.add(key)
            if binding.site_id not in centers[binding.metal_atom_id].site_ids:
                raise FormalError("Binding site was not declared on the metal")
        ionic_fragment_charge=sum(atoms[i].formal_charge for i in ligand.atom_ids)-(1 if binding.mode=="sigma_covalent" else 0)
        if ligand.ionic_charge is not None and ligand.ionic_charge!=ionic_fragment_charge:
            raise Unsupported("Ligand ionic charge disagrees with this adapter's M-X ionic approximation; noninnocent/charged L ligands need another profile")
        if binding.mode=="sigma_covalent":
            if len(donors)!=1 or binding.hapticity!=1 or ligand.ionic_class!="X":
                raise Unsupported("Only monodentate X sigma-covalent ligands are executable")
            pair=tuple(sorted((binding.metal_atom_id,next(iter(donors)))))
            if orders.get(pair)!=1: raise FormalError("M-X binding must correspond to one ordinary single bond")
            covered.add(pair)
        else:
            if ligand.ionic_class!="L": raise Unsupported("Dative/eta adapter requires neutral L ligands")
            if any(tuple(sorted((binding.metal_atom_id,i))) in orders for i in ligand.atom_ids):
                raise FormalError("Dative/eta binding cannot also be counted as ordinary M-ligand bonds")
            if binding.mode=="dative":
                if len(donors)!=1 or binding.hapticity!=1:
                    raise Unsupported("Only one-donor 2e L dative bindings are supported")
                reserved[next(iter(donors))]+=2
            elif binding.hapticity!=2 or len(donors)!=2 or any(atoms[i].element!="C" for i in donors) or orders.get(tuple(sorted(donors)))!=2:
                raise Unsupported("Only one neutral eta2 alkene pi pair is supported")
        # A bound ligand is complete: ordinary organic neighbors cannot be hidden.
        for i in ligand.atom_ids:
            if any(j not in ligand.atom_ids and j!=binding.metal_atom_id for pair in orders if i in pair for j in pair if j!=i):
                raise FormalError("Bound ligand definition omits an attached organic atom")
    metal_bonds={p for p in orders if set(p)&metal_ids}
    if any(set(p)<=metal_ids for p in metal_bonds):
        raise Unsupported("Metal-metal bonding requires a separate multimetal adapter")
    if metal_bonds!=covered:
        raise FormalError("Every covalent metal bond requires a sigma_covalent ligand binding")
    binding_ids={b.binding_id for b in view.metal_bindings}
    for m in centers.values():
        if len(m.site_ids)!=len(set(m.site_ids)): raise FormalError("Duplicate metal sites")
        vacancies=m.vacant_site_ids
        if vacancies is not None and (len(vacancies)!=len(set(vacancies)) or not set(vacancies)<=set(m.site_ids) or any((m.atom_id,s) in sites for s in vacancies)):
            raise FormalError("Vacancy references an occupied, repeated or undeclared site")
        own={b.binding_id for b in view.metal_bindings if b.metal_atom_id==m.atom_id}
        for relation in m.relative_positions:
            if relation.first_binding_id==relation.second_binding_id or not {relation.first_binding_id,relation.second_binding_id}<=own:
                raise FormalError("Relative-position constraint references invalid bindings")
    for i,n in reserved.items():
        if atoms[i].nonbonding_electrons is None or atoms[i].nonbonding_electrons-atoms[i].radical_electrons<n:
            raise FormalError("Dative binding has no unspent donor-owned lone pair")
    return centers, ligands, reserved


def metal_counts(view):
    centers,ligands,_=binding_checks(view)
    atoms=atoms_by_id(view)
    reports=[]
    for i,m in centers.items():
        bindings=[b for b in view.metal_bindings if b.metal_atom_id==i]
        members={i}|{a for b in bindings for a in ligands[b.ligand_id].atom_ids}
        q=sum(atoms[a].formal_charge for a in members)
        known=all(ligands[b.ligand_id].ionic_charge is not None for b in bindings)
        os=q-sum(ligands[b.ligand_id].ionic_charge for b in bindings) if known else None
        d=GROUP[atoms[i].element]-os if os is not None else None
        ve=d+2*len(bindings) if d is not None else None
        for name,value in (("oxidation_state",os),("d_electrons",d),("valence_electrons",ve)):
            stored=getattr(m,name)
            if stored is not None and value is None: raise Unsupported(f"Cannot check declared {name} with unknown ligand assignment")
            if stored is not None and stored!=value: raise FormalError(f"Metal {i} {name} mismatch: declared {stored}, ionic calculation {value}")
        if d is not None and not 0<=d<=10: raise Unsupported("Formal d count outside supported d0..d10 metal profiles")
        reports.append({"metal_atom_id":i,"element":atoms[i].element,"complex_atom_ids":sorted(members),
            "complex_total_charge":q,"atomic_formal_charge":atoms[i].formal_charge,
            "oxidation_state":os,"d_electrons":d,"valence_electrons":ve,"counting_method":"ionic",
            "coordination_units":len(bindings),"eighteen_electron_rule":"advisory_only",
            "eighteen_electron_hint":"18e" if ve==18 else "non-18e is not a formal failure",
            "global_ledger_electrons_added":0,"geometry_status":"declared_only" if m.geometry!="unknown" else "unknown"})
    return reports


def validate_view(view, state=None):
    centers,_,_=binding_checks(view)
    if not centers:
        result=validate02(view,state)
        return {**result,"metal_counts":[],"rule_version":"umr-metal-ledger/0.3"}
    vals=bond_valences(view)
    for atom in view.atoms:
        if atom.formal_charge!=charge(atom,vals[atom.atom_id]):
            raise FormalError(f"Charge ledger mismatch at atom {atom.atom_id}")
        if atom.nonbonding_electrons<0: raise FormalError("Negative nonbonding electron count")
        if atom.atom_id in centers and (atom.radical_electrons or atom.nonbonding_electrons%2 or atom.hydrogens):
            raise Unsupported("Metal radicals, unspecified metal hydrogens and odd metal electron ownership are unsupported")
    # Validate organic valence with private bond-valence caps. No cap is saved in a
    # ChemicalState or counted globally; an M-X single bond supplies one valence.
    organic=[a.model_copy(deep=True) for a in view.atoms if a.atom_id not in centers]
    if organic:
        cap_atoms=[];cap_bonds=[];cap_id=max(a.atom_id for a in view.atoms)+1
        for a in organic:
            count=sum(b.order for b in view.bonds if a.atom_id in (b.a,b.b) and ({b.a,b.b}&centers.keys()))
            if a.element=="H" and count:
                # A numbered metal hydride is still H, not an H with an attached
                # unnumbered hydrogen. Use an explicit PRIVATE valence cap.
                if count!=1: raise FormalError("Metal-bound hydrogen must have one single bond")
                cap_atoms.append(Atom(atom_id=cap_id,element="H",formal_charge=0,nonbonding_electrons=0))
                cap_bonds.append(Bond(a=a.atom_id,b=cap_id,order=1));cap_id+=1
            else: a.hydrogens+=count
        shadow=View01(state_id=view.state_id,view_id=view.view_id,atoms=organic+cap_atoms,
                      bonds=[b for b in view.bonds if not {b.a,b.b}&centers.keys()]+cap_bonds)
        validate02(shadow)
    q=sum(a.formal_charge for a in view.atoms)
    radicals=sum(a.radical_electrons for a in view.atoms)
    if state:
        if state.state_id!=view.state_id or state.total_charge!=q: raise FormalError("State/view charge or identity mismatch")
        if state.electronic_mode!=("RADICAL" if radicals else "CLOSED_SHELL"): raise FormalError("State electronic_mode mismatch")
        if state.multiplicity is not None and state.multiplicity!= (2 if radicals==1 else 1):
            raise Unsupported("Only closed-shell metal plus singlet/one-organic-radical doublet spin profiles supported")
    return {"total_charge":q,"electrons":electron_count(view),"metal_counts":metal_counts(view),
            "chemical_feasibility":"unknown","rule_version":"umr-metal-ledger/0.3"}
