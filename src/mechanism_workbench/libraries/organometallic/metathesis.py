"""Source-pinned catalyst identities; structure, initiation and scope are separate.

The graphs are explicit finite molecular models. A family name or a metal
symbol is never sufficient to assign these ligands or stereoselectivity.
"""
from copy import deepcopy
from rdkit import Chem,rdBase

VERSION='umr-metathesis-catalysts/1.0'
PCY='P(C1CCCCC1)(C1CCCCC1)C1CCCCC1'
NHC='[C]1N(c2c(C)cc(C)cc2C)CCN1c1c(C)cc(C)cc1C'
PY='n1cc(Br)ccc1'
SCHROCK='CC(C)(c1ccccc1)C=[M](=Nc1c(C(C)C)cccc1C(C)C)(OC(C)(C(F)(F)F)C(F)(F)F)OC(C)(C(F)(F)F)C(F)(F)F'


def ru(p=False,nhc=True,chelate=False,nitro=False,pyridines=False):
    benzyl='=Cc1cc([N+](=O)[O-])ccc1OC(C)C' if nitro else '=Cc1ccccc1OC(C)C' if chelate else '=Cc1ccccc1'
    raw='Cl[Ru](Cl)('+benzyl+')'
    ligands=([PCY] if p else [])+[NHC if nhc else PCY]+([PY,PY] if pyridines else [])
    raw+=''.join('(<-'+ligand+')' for ligand in ligands)
    with rdBase.BlockLogs():m=Chem.MolFromSmiles(raw)
    if m is None:raise ValueError('Invalid finite Ru catalyst graph')
    if chelate:
        w=Chem.RWMol(m);metal=next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=='Ru')
        donors=[a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()=='O' and len(a.GetNeighbors())==2 and
                any(n.GetIsAromatic() for n in a.GetNeighbors())]
        if len(donors)!=1:raise ValueError('Chelated styrenyl ether must have one actual O donor')
        w.AddBond(donors[0],metal,Chem.BondType.DATIVE);m=w.GetMol();m.UpdatePropertyCache(strict=False)
    return Chem.MolToSmiles(m,isomericSmiles=True)


def records():
    entries=[
        ('Grubbs I',['Grubbs-I','GrubbsI','G-I',"Grubbs' I"],ru(p=True,nhc=False),'phosphine_dissociation','10.1021/ja952676d'),
        ('Grubbs II',['Grubbs-II','GrubbsII','G-II',"Grubbs' II","Grubbs'II",'Grubbs 2nd generation'],ru(p=True),'phosphine_dissociation','10.1021/ol990909q'),
        ('Hoveyda-Grubbs I',['Hoveyda-Grubbs-I','HG-I','HG I','Hoveyda-Grubbs(I)'],ru(nhc=False,chelate=True),'ether_olefin_exchange','10.1021/ja983222u'),
        ('Hoveyda-Grubbs II',['Hoveyda-Grubbs-II','HG-II','HG II','Hoveyda-Grubbs(II)'],ru(chelate=True),'ether_olefin_exchange','10.1021/ja001179g'),
        ('Grela catalyst',['Grela','nitro-Grela'],ru(chelate=True,nitro=True),'ether_olefin_exchange','10.1021/ja048794v'),
        ('Grubbs III',['Grubbs-III','G-III','Grubbs third generation'],ru(pyridines=True),'pyridine_dissociation','10.1039/C9FD00043G'),
        ('Schrock Mo imido',['Schrock Mo','Mo imido alkylidene'],SCHROCK.replace('[M]','[Mo]'),'active_alkylidene','10.1021/ja00166a023'),
        ('Schrock W imido',['W imido alkylidene'],SCHROCK.replace('[M]','[W]'),'active_alkylidene','10.1021/om00158a025'),
    ]
    out=[]
    for ident,aliases,graph,init,doi in entries:
        m=Chem.MolFromSmiles(graph,sanitize=False)
        element=next(a.GetSymbol() for a in m.GetAtoms() if a.GetSymbol() in {'Ru','Mo','W'})
        out.append({'catalyst_id':ident,'aliases':aliases,'structure_smiles':graph,'metal':element,
            'initiation_profile':init,'carbene_convention':'neutral_L' if element=='Ru' else 'X2',
            'precatalyst_oxidation_state':2 if element=='Ru' else 6,'total_spin':'unknown',
            'supported_events':['olefin_association','metallacyclobutane_formation','metallacyclobutane_cleavage','product_release'],
            'executability':'finite_olefin_model','selectivity':'unknown unless separately evidenced',
            'sources':[{'doi':doi,'scope':'catalyst identity or family; not a universal initiation rate law'},
                       {'doi':'10.1021/ja010624k' if element=='Ru' else '10.1016/S0040-4020(99)00304-X','scope':'mechanistic family'}],
            'initiation_uncertainty':'Hoveyda exchange is one evidence-compatible grouped coordination channel; associative/interchange/dissociative weights are not predicted.' if init=='ether_olefin_exchange' else 'kinetics and detailed solvent coordination unknown',
            'structure_evidence':'source-family finite molecular model; exact ligand/rotamer state is an explicit modelling assumption' if element=='W' else 'named finite molecular graph; coordination/spin is a declared Lewis model',
            'human_approved':False})
    for ident,scope,doi in [
        ('Z-stereoretentive Ru dithiolate','Requires exact dithiolate/NHC structure and a registered stereo-specific profile','10.1038/s41467-024-52876-4'),
        ('Ru enyne metathesis','Dedicated alkene/alkyne electron topology required; ordinary olefin rule does not cover it','10.1021/ja050174m'),
        ('Mo/W alkylidyne metathesis','Alkyne X3 rules remain distinct from olefin X2/L rules','10.1021/ja00166a023')]:
        out.append({'catalyst_id':ident,'aliases':[],'structure_smiles':None,'executability':'requires_specific_structure_or_profile',
                    'scope':scope,'sources':[{'doi':doi,'scope':'family/context only'}],'human_approved':False})
    return out


def get(identity):
    for rec in records():
        if identity==rec['catalyst_id'] or identity in rec['aliases']:return deepcopy(rec)
    raise ValueError('Unregistered or ambiguous metathesis catalyst identity: '+str(identity))


def active_methylidene_graph(identity):
    """A separately declared finite active-species MODEL, not a reaction that
    deletes precursor atoms. Real initiation is always executed elsewhere.
    """
    rec=get(identity);m=Chem.MolFromSmiles(rec['structure_smiles']);w=Chem.RWMol(m)
    metal=next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol()==rec['metal'])
    alk=next(b.GetOtherAtomIdx(metal) for b in m.GetAtomWithIdx(metal).GetBonds()
             if b.GetBondType()==Chem.BondType.DOUBLE and b.GetOtherAtom(m.GetAtomWithIdx(metal)).GetSymbol()=='C')
    def organic(root):
        ids={root};stack=[root]
        while stack:
            i=stack.pop()
            for n in m.GetAtomWithIdx(i).GetNeighbors():
                j=n.GetIdx()
                if j!=metal and j not in ids:ids.add(j);stack.append(j)
        return ids
    removed=organic(alk)
    labile=[b.GetBeginAtomIdx() for b in m.GetBonds() if b.GetBondType()==Chem.BondType.DATIVE and b.GetEndAtomIdx()==metal
            and b.GetBeginAtom().GetSymbol() in {'P','N'}]
    if rec['initiation_profile']=='phosphine_dissociation':
        for d in labile[:1]:removed|=organic(d)
    elif rec['initiation_profile']=='pyridine_dissociation':
        for d in labile:removed|=organic(d)
    new_metal=metal-sum(i<metal for i in removed)
    for i in sorted(removed,reverse=True):w.RemoveAtom(i)
    a=Chem.Atom('C');a.SetNumExplicitHs(2);a.SetNoImplicit(True)
    c=w.AddAtom(a);w.AddBond(new_metal,c,Chem.BondType.DOUBLE)
    out=w.GetMol();out.UpdatePropertyCache(strict=False)
    return Chem.MolToSmiles(out,isomericSmiles=True)


def register(entries):
    for rec in records():
        if not rec.get('structure_smiles'):continue
        ident=rec['catalyst_id'];old=entries.get(ident,{})
        entries[ident]={**old,'entity_id':ident,'aliases':list(dict.fromkeys([ident,*old.get('aliases',[]),*rec['aliases']])),
            'structure_smiles':rec['structure_smiles'],'roles':['metathesis_catalyst'],'kind':'concrete_reagent',
            'identity_status':'registered_identity','structure_status':'source_pinned_finite_metathesis_graph',
            'sources':rec['sources'],'model_version':VERSION}


def decorate(view,metal,*,convention=None):
    """Classify an independently present complete graph, preserving its atoms.

    Actual electron ownership is untouched, except that the explicitly named
    bound NHC's paired occupation is recorded with evidence. Total spin stays
    unknown. No free carbene template permissions are granted.
    """
    from ...umr.v09.ast import Ligand,MetalBinding
    from ...umr.v04.ast import ReactiveCenter
    from ...umr.chemistry import atoms_by_id,bond_orders
    ats=atoms_by_id(view);edges=bond_orders(view);mids={x.atom_id for x in view.metal_centers}
    bound_lids={b.ligand_id for b in view.metal_bindings if b.metal_atom_id==metal}
    if any(b.metal_atom_id!=metal and b.ligand_id in bound_lids for b in view.metal_bindings):raise ValueError('Metathesis catalyst cannot silently replace a bridging ligand')
    old_datives=[i for b in view.metal_bindings if b.metal_atom_id==metal and b.mode=='dative' for i in b.donor_atom_ids]
    view.ligands=[l for l in view.ligands if l.ligand_id not in bound_lids]
    view.metal_bindings=[b for b in view.metal_bindings if b.metal_atom_id!=metal]
    contacts=[(b if a==metal else a,o) for (a,b),o in edges.items() if metal in (a,b)]
    roots={i for i,_ in contacts}|set(old_datives);seen=set()
    for root in sorted(roots):
        if root in seen:continue
        ids={root};stack=[root]
        while stack:
            i=stack.pop()
            for a,b in edges:
                if i not in (a,b):continue
                j=b if i==a else a
                if j not in mids|ids:ids.add(j);stack.append(j)
        seen|=ids;cov=[(i,n) for i,n in contacts if i in ids];dat=[i for i in old_datives if i in ids]
        alk=[i for i,n in cov if n==2 and ats[i].element=='C']
        neutral=bool(alk and ats[metal].element=='Ru' and convention!='X2')
        q=sum(ats[i].formal_charge for i in ids)-sum(n for i,n in cov)+(2 if neutral else 0)
        klass={(0,0):'L',(0,1):'L',(0,2):'L2',(-1,0):'X',(-2,0):'X2'}.get((q,len(dat)+(1 if neutral else 0)))
        if klass is None:raise ValueError('Metathesis ligand partition outside finite model')
        typ='Ru_alkylidene' if neutral else 'Schrock_alkylidene' if alk else 'generic'
        for i in dat:
            ns=[ats[j].element for edge in edges if i in edge for j in edge if j!=i]
            if ats[i].element=='C' and sorted(ns)==['N','N']:
                if ats[i].nonbonding_electrons!=2 or ats[i].formal_charge:raise ValueError('Named NHC pair does not match graph')
                ats[i].radical_electrons=0;typ='NHC'
                if not any(c.atom_id==i for c in view.reactive_centers):
                    view.reactive_centers.append(ReactiveCenter(center_id='metathesis.nhc.'+str(i),atom_id=i,
                        species='carbene',profile='carbene_paired',environment='metal_bound',existence_basis='formal_model',
                        reaction_context='Named coordinated NHC ligand, not an independently generated free carbene',
                        evidence=[{'level':'proposed_completion','source':'doi:10.1021/ol990909q',
                            'note':'Bound neutral NHC with two N neighbours and one donor-owned pair; formal finite ligand assignment, total spin unknown.'}]))
        lid='om.ligand.'+str(min(ids));donors=sorted({i for i,_ in cov}|set(dat))
        view.ligands.append(Ligand(ligand_id=lid,atom_ids=sorted(ids),donor_atom_ids=donors,
            ionic_class=klass,ionic_charge=q,ligand_type=typ,assignment_basis=VERSION+'; explicit carbene convention, total spin unknown'))
        for i,n in cov:
            view.metal_bindings.append(MetalBinding(binding_id=lid+'.sigma.'+str(i),ligand_id=lid,metal_atom_id=metal,
                donor_atom_ids=[i],mode='sigma_covalent',donated_electrons=2 if neutral else 2*n,ionic_allocation=q))
        for i in dat:
            view.metal_bindings.append(MetalBinding(binding_id=lid+'.dative.'+str(i),ligand_id=lid,metal_atom_id=metal,
                donor_atom_ids=[i],mode='dative',donated_electrons=2,ionic_allocation=0))
    for c in view.metal_centers:
        if c.atom_id==metal:c.oxidation_state=None;c.d_electrons=None;c.valence_electrons=None;c.relative_positions=[];c.vacant_site_ids=None
    view.annotations['metathesis_representation']={'version':VERSION,'metal_atom_id':metal,
        'carbene_convention':convention or ('neutral_L' if ats[metal].element=='Ru' else 'X2'),'spin':'unknown'}
