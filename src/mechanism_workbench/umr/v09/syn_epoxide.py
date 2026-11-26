"""Relative facial constraint for the concerted peracid oxygen transfer only."""
from rdkit import Chem
from .smiles import to_rdkit
from ..chemistry import FormalError

KEY='peracid_syn_centers_v1'

def parity(atom, ordering):
    tag=atom.GetChiralTag()
    if tag not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW,Chem.ChiralType.CHI_TETRAHEDRAL_CCW):return None
    local=[a.GetIdx() for a in atom.GetNeighbors()]
    if len(local)==3:local.append(-1)
    rest=set(local)-set(ordering)
    if len(local)!=4 or len(rest)!=1:return None
    perm=[local.index(i) for i in [*ordering,*rest]]
    inv=sum(a>b for j,a in enumerate(perm) for b in perm[j+1:])
    return (1 if tag==Chem.ChiralType.CHI_TETRAHEDRAL_CW else -1)*(-1)**inv

def check(source,target,roles):
    if not isinstance(roles,list) or len(roles)!=3 or len(set(roles))!=3:raise FormalError('Three distinct epoxide roles required')
    before,after=to_rdkit(source),to_rdkit(target)
    for m in (before,after):Chem.SanitizeMol(m)
    ib={a.GetAtomMapNum():a.GetIdx() for a in before.GetAtoms()}
    ia={a.GetAtomMapNum():a.GetIdx() for a in after.GetAtoms()}
    if not all(i in ib and i in ia for i in roles):raise FormalError('Unknown syn epoxide role')
    a,b,o=roles
    bond=before.GetBondBetweenAtoms(ib[a],ib[b])
    def order(i,j):
        z=after.GetBondBetweenAtoms(ia[i],ia[j]);return z.GetBondTypeAsDouble() if z else 0
    if not bond or bond.GetBondTypeAsDouble()!=2 or any(order(i,j)!=1 for i,j in [(a,b),(b,o),(a,o)]):
        raise FormalError('Syn epoxide constraint requires a C=C to epoxide transaction')
    refs=list(bond.GetStereoAtoms())
    st=bond.GetStereo()
    if refs and st in (Chem.BondStereo.STEREOE,Chem.BondStereo.STEREOZ,Chem.BondStereo.STEREOCIS,Chem.BondStereo.STEREOTRANS):
        cis=st in (Chem.BondStereo.STEREOZ,Chem.BondStereo.STEREOCIS)
        if bond.GetBeginAtomIdx()!=ib[a]:refs.reverse()
    else:
        rings=[r for r in before.GetRingInfo().AtomRings() if 3<=len(r)<=7 and ib[a] in r and ib[b] in r
               and abs(r.index(ib[a])-r.index(ib[b])) in (1,len(r)-1)]
        if not rings:return 'input_geometry_unknown'
        r=min(rings,key=len)
        refs=[next(j for j in (r[r.index(i)-1],r[(r.index(i)+1)%len(r)]) if j!=other)
              for i,other in ((ib[a],ib[b]),(ib[b],ib[a]))];cis=True
    idx=[ia[before.GetAtomWithIdx(i).GetAtomMapNum()] for i in refs]
    signs=[parity(after.GetAtomWithIdx(ia[p]),[ia[q],ia[o],ref]) for p,q,ref in ((a,b,idx[0]),(b,a,idx[1]))]
    if None in signs:return 'product_relative_geometry_unspecified'
    if cis!=(signs[0]!=signs[1]):raise FormalError('Product epoxide contradicts syn addition to the input alkene')
    return 'relative_geometry_consistent'

