"""A measured organic endpoint with an atom-resolved, structurally unknown remainder."""
from collections import Counter
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from ..chemistry import FormalError, atoms_by_id, order
from ..v09.coordination import checks
from ..v09.ledger import validate_view
from ..v09.smiles import to_rdkit


def component(view, start):
    adjacency = {a.atom_id: set() for a in view.atoms}
    for bond in view.bonds:
        adjacency[bond.a].add(bond.b)
        adjacency[bond.b].add(bond.a)
    seen, pending = {start}, [start]
    while pending:
        new = adjacency[pending.pop()] - seen
        seen.update(new)
        pending.extend(new)
    return seen


def mol_subset(view, ids):
    """Construct an organic graph while preserving the atom-neighbour stereo convention."""
    mol = Chem.RWMol(to_rdkit(view))
    for atom in reversed(list(mol.GetAtoms())):
        if atom.GetAtomMapNum() not in ids:
            mol.RemoveAtom(atom.GetIdx())
    return mol


def element_inventory(atoms):
    counts = Counter()
    for atom in atoms:
        counts[(atom.element, atom.isotope or 0)] += 1
        if atom.hydrogens:
            counts[("H", 0)] += atom.hydrogens
    return {((str(mass) if mass else "") + el): n for (el, mass), n in sorted(counts.items())}


def evaluate(view, observation):
    validate_view(view)
    a = atoms_by_id(view)
    c, o, ch = observation.carbonyl, observation.carbonyl_oxygen, observation.methylene
    zn = list(observation.zinc_ends)
    ti = observation.titanium
    selected = [c, o, ch, *zn, ti]
    if len(set(selected)) != 6 or any(i not in a for i in selected):
        raise FormalError("Net observation needs six distinct atoms in the initial inventory")
    if [a[i].element for i in selected] != ["C", "O", "C", "Zn", "Zn", "Ti"]:
        raise FormalError("Expected carbonyl C/O, bridging methylene, two Zn and a Ti activator")
    counts = {m["metal_atom_id"]: m for m in checks(view)[1]}
    if not {ti,*zn} <= counts.keys():
        raise FormalError("Declare both Zn(II) donors and the Ti(IV) activator as metal centers")
    neighbours = lambda i: {b.b if b.a == i else b.a for b in view.bonds if i in (b.a, b.b)}
    if (order(view, c, o) != 2 or neighbours(o) != {c} or len(neighbours(c)) != 3
            or any(a[i].element != "C" for i in neighbours(c) - {o})):
        raise FormalError("A nonprotonated ketone carbonyl with two carbon substituents is required")
    if (a[ch].hydrogens != 2 or neighbours(ch) != set(zn)
            or any(order(view, ch, i) != 1 or counts[i]["oxidation_state"] != 2 for i in zn)):
        raise FormalError("The methylene must bridge the two declared Zn(II) centers")
    cluster = component(view, ch)
    zinc = {i for i in cluster if a[i].element == "Zn"}
    bridges = {i for i in cluster if a[i].element == "C" and a[i].hydrogens == 2
               and len(neighbours(i)) == 2 and neighbours(i) <= zinc}
    if (len(zinc) != 3 or len(bridges) != 2 or sum(a[i].element == "Br" for i in cluster) != 2
            or not any(a[i].element == "O" and len({b.metal_atom_id for b in view.metal_bindings
                       if b.mode == "dative" and i in b.donor_atom_ids} & zinc) == 2 for i in a)):
        raise FormalError("Expected the THF-bridged trizinc Nysted reagent")
    if (counts[ti]["oxidation_state"] != 4 or len(neighbours(ti)) != 4
            or any(a[i].element != "Cl" or order(view, ti, i) != 1 for i in neighbours(ti))):
        raise FormalError("The recorded activator is titanium(IV) chloride")
    substrate = component(view, c)
    if substrate & cluster:
        raise FormalError("The recorded ketone and methylene donor must be distinct components")
    product_ids = (substrate - {o}) | {ch}
    remainder = list(observation.unassigned_atom_ids)
    if len(set(remainder)) != len(remainder) or set(remainder) != set(a) - product_ids:
        raise FormalError("Incomplete or duplicated atom inventory across product and unresolved remainder")
    mol = mol_subset(view, product_ids)
    indices = {atom.GetAtomMapNum(): atom.GetIdx() for atom in mol.GetAtoms()}
    mol.AddBond(indices[c], indices[ch], Chem.BondType.DOUBLE)
    Chem.SanitizeMol(mol)
    mapped = Chem.MolToSmiles(mol, isomericSmiles=True)
    product_charge = Chem.GetFormalCharge(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return dict(id=observation.observation_id, source_view=view.view_id, profile=observation.profile,
                status="organic_net_matched_inventory", microscopic_steps=0,
                mapped_product=mapped, product=Chem.MolToSmiles(mol, isomericSmiles=True),
                formula=rdMolDescriptors.CalcMolFormula(mol), product_atom_ids=sorted(product_ids),
                carbonyl_carbon=c, transferred_carbon=ch, carbonyl_oxygen=o, donor_zinc=zn, activator=ti,
                remainder=dict(atom_ids=sorted(remainder), elements=element_inventory(a[i] for i in remainder),
                               charge=sum(atom.formal_charge for atom in a.values()) - product_charge,
                               chemical_species=None, oxygen_atom_id=o),
                evidence=[e.model_dump(mode="json") for e in observation.evidence])
