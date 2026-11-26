"""Mapped SMILES is derived; the AST retains information SMILES cannot carry."""

from rdkit import Chem, rdBase

from .ast import Atom, Bond, ChemicalState, LewisView
from .chemistry import FormalError, Unsupported, VALENCE

RDKIT_VERSION = rdBase.rdkitVersion


def import_mapped_smiles(smiles: str, state_id: str = "s0", view_id: str = "v0",
                         labels: dict[int, str] | None = None, *, _preserve_lewis=False) -> tuple[ChemicalState, LewisView]:
    params = Chem.SmilesParserParams()
    params.removeHs = False
    params.sanitize = not _preserve_lewis
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None:
        raise FormalError("Invalid SMILES")
    if _preserve_lewis:
        # Do not re-aromatize an explicitly supplied Kekule contributor before
        # reading it: that can swap all alternating bonds while preserving maps.
        Chem.SanitizeMol(mol,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    maps = [a.GetAtomMapNum() for a in mol.GetAtoms()]
    if any(x <= 0 for x in maps) or len(maps) != len(set(maps)):
        raise FormalError("Every explicit atom, including transferred H, must have a unique positive map/atom_id")
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as exc:
        raise Unsupported("Aromatic graph could not be represented as an integer Lewis view") from exc
    atoms = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        hydrogens = atom.GetTotalNumHs(includeNeighbors=False)
        bond_sum = sum(b.GetBondTypeAsDouble() for b in atom.GetBonds())
        if any(b.GetBondType() not in {Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE}
               for b in atom.GetBonds()):
            raise Unsupported("Coordination/dative/fractional bonds are not implemented")
        nonbonding = VALENCE[symbol] - atom.GetFormalCharge() - int(bond_sum) - hydrogens if symbol in VALENCE else None
        if nonbonding is not None and nonbonding < 0:
            raise FormalError("Negative nonbonding electron count in imported structure")
        tag = str(atom.GetChiralTag())
        stereo = ({"status": "specified", "rdkit_tag": tag,
                   "neighbor_order": [n.GetAtomMapNum() for n in atom.GetNeighbors()]}
                  if tag != "CHI_UNSPECIFIED" else {"status": "unknown"})
        atoms.append(Atom(atom_id=atom.GetAtomMapNum(), element=symbol, isotope=atom.GetIsotope(),
                          formal_charge=atom.GetFormalCharge() if symbol != "*" else None,
                          nonbonding_electrons=nonbonding, hydrogens=hydrogens,
                          radical_electrons=atom.GetNumRadicalElectrons(), stereo=stereo,
                          label=(labels or {}).get(atom.GetAtomMapNum(), "R" if symbol == "*" else None)))
    bonds = []
    for bond in mol.GetBonds():
        stereo = {"status": "unknown"}
        if bond.GetStereo() != Chem.BondStereo.STEREONONE:
            stereo = {"status": "specified", "rdkit_stereo": str(bond.GetStereo()),
                      "reference_atoms": [maps[i] for i in bond.GetStereoAtoms()]}
        bonds.append(Bond(a=maps[bond.GetBeginAtomIdx()], b=maps[bond.GetEndAtomIdx()],
                          order=int(bond.GetBondTypeAsDouble()), stereo=stereo))
    view = LewisView(view_id=view_id, state_id=state_id, atoms=atoms, bonds=bonds,
                     annotations={"imported_from_mapped_smiles": smiles, "rdkit_version": RDKIT_VERSION})
    state = ChemicalState(state_id=state_id, primary_view_id=view_id,
                          total_charge=sum(a.GetFormalCharge() for a in mol.GetAtoms()))
    return state, view


def to_rdkit(view: LewisView, *, _sanitize_ops=None):
    rw = Chem.RWMol()
    indices = {}
    for atom in view.atoms:
        rd_atom = Chem.Atom(atom.element)
        rd_atom.SetAtomMapNum(atom.atom_id)
        rd_atom.SetIsotope(atom.isotope)
        rd_atom.SetFormalCharge(atom.formal_charge or 0)
        rd_atom.SetNumExplicitHs(atom.hydrogens)
        rd_atom.SetNoImplicit(True)
        if atom.radical_electrons is None:
            raise Unsupported("Unknown occupation requires the JSON sidecar; deterministic spin depiction unavailable")
        rd_atom.SetNumRadicalElectrons(atom.radical_electrons)
        indices[atom.atom_id] = rw.AddAtom(rd_atom)
    for bond in view.bonds:
        rw.AddBond(indices[bond.a], indices[bond.b],
                   {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[bond.order])
    mol = rw.GetMol()
    if _sanitize_ops is not None:
        # Version adapters may preserve an explicitly registered Lewis graph.
        # They must validate it themselves; this never changes legacy defaults.
        Chem.SanitizeMol(mol,sanitizeOps=_sanitize_ops)
    elif hasattr(view,"reactive_centers"):
        # RDKit FINDRADICALS otherwise changes a declared paired carbene into a
        # two-radical encoding. Keep authority's explicit occupation and H counts.
        Chem.SanitizeMol(mol,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_FINDRADICALS ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        if any(a.GetNumRadicalElectrons()!=original.radical_electrons or a.GetNumImplicitHs()!=0
               for a,original in zip(mol.GetAtoms(),view.atoms)):
            raise Unsupported("RDKit changed authoritative radical/H occupation")
    else:
        Chem.SanitizeMol(mol)
    for atom in view.atoms:
        if atom.stereo.get("status") != "specified":
            continue
        tag = atom.stereo.get("rdkit_tag")
        if tag not in {"CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW"}:
            raise Unsupported("Only tetrahedral atom stereo is implemented by the SMILES adapter")
        ref = atom.stereo["neighbor_order"]
        current = [a.GetAtomMapNum() for a in mol.GetAtomWithIdx(indices[atom.atom_id]).GetNeighbors()]
        if set(ref) != set(current):
            raise Unsupported("Stereo neighbor references changed")
        permutation = [ref.index(x) for x in current]
        parity = sum(permutation[i] > permutation[j] for i in range(len(permutation)) for j in range(i+1, len(permutation))) % 2
        if parity:
            tag = "CHI_TETRAHEDRAL_CCW" if tag == "CHI_TETRAHEDRAL_CW" else "CHI_TETRAHEDRAL_CW"
        mol.GetAtomWithIdx(indices[atom.atom_id]).SetChiralTag(getattr(Chem.ChiralType, tag))
    for bond in view.bonds:
        if bond.stereo.get("status") == "specified":
            rd_bond = mol.GetBondBetweenAtoms(indices[bond.a], indices[bond.b])
            refs = bond.stereo.get("reference_atoms", [])
            if len(refs) == 2:
                rd_bond.SetStereoAtoms(indices[refs[0]], indices[refs[1]])
                rd_bond.SetStereo(getattr(Chem.BondStereo, bond.stereo["rdkit_stereo"]))
    Chem.AssignStereochemistry(mol, cleanIt=False, force=True)
    # AST double-bond stereo alone is insufficient for an isomeric SMILES writer:
    # neighboring single-bond directions must be regenerated (audit sample 053).
    # This repairs a derived serialization loss; no UMR AST semantics are changed.
    Chem.SetDoubleBondNeighborDirections(mol)
    return mol


def export_mapped_smiles(view: LewisView) -> str:
    return Chem.MolToSmiles(to_rdkit(view), canonical=False, isomericSmiles=True, kekuleSmiles=True)
