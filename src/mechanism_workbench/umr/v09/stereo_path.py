"""Event-local geometry and explicit, outcome-isolated path constraints.

No R/S label is toggled by reaction name. Parity is relative to stable mapped
neighbours; the source and executed states remain the authority.
"""
from copy import deepcopy
from ..chemistry import FormalError,atoms_by_id,bond_orders

VERSION='umr-stereo-path/1.0'
TAGS={'CHI_TETRAHEDRAL_CW','CHI_TETRAHEDRAL_CCW'}


def potential_centers(view):
    from rdkit import Chem,rdBase
    from .smiles import to_rdkit
    try:
        with rdBase.BlockLogs():mol=to_rdkit(view)
        ids=[a.GetAtomMapNum() for a in mol.GetAtoms()]
        for a in mol.GetAtoms():a.SetAtomMapNum(0)
        return {ids[s.centeredOn] for s in Chem.FindPotentialStereo(mol,cleanIt=False,flagPossible=True)
                if s.type==Chem.StereoType.Atom_Tetrahedral}
    except (ValueError,RuntimeError):return {a.atom_id for a in view.atoms if a.stereo.get('status')=='specified'}


def potential_double_bonds(view):
    from rdkit import Chem,rdBase
    from .smiles import to_rdkit
    try:
        with rdBase.BlockLogs():mol=to_rdkit(view)
        ids=[a.GetAtomMapNum() for a in mol.GetAtoms()]
        for a in mol.GetAtoms():a.SetAtomMapNum(0)
        out=set()
        for s in Chem.FindPotentialStereo(mol,cleanIt=False,flagPossible=True):
            if s.type==Chem.StereoType.Bond_Double:
                b=mol.GetBondWithIdx(s.centeredOn)
                out.add(tuple(sorted((ids[b.GetBeginAtomIdx()],ids[b.GetEndAtomIdx()]))))
        return out
    except (ValueError,RuntimeError):
        return {tuple(sorted((b.a,b.b))) for b in view.bonds if b.order==2 and b.stereo.get('status')=='specified'}


def adjacent(view,i):
    return {b.b if b.a==i else b.a:b.order for b in view.bonds if i in (b.a,b.b)}


def parity(stereo,view,i):
    if stereo.get('status')!='specified' or stereo.get('rdkit_tag') not in TAGS:return None
    refs=list(stereo.get('neighbor_order',[]));ns=adjacent(view,i)
    if not refs or len(set(refs))!=len(refs) or not set(refs)<=ns.keys():return None
    if len(refs)==3:
        extra=set(ns)-set(refs)
        if len(extra)==1:refs.append(next(iter(extra)))
        elif not extra and atoms_by_id(view)[i].hydrogens==1:refs.append(-i)
    if len(refs)!=4:return None
    return sum(a>b for j,a in enumerate(refs) for b in refs[j+1:])%2 ^ (stereo['rdkit_tag']=='CHI_TETRAHEDRAL_CW')


def stable(source,target,i):
    a,b=atoms_by_id(source)[i],atoms_by_id(target)[i]
    return (a.element=='C' and b.element=='C' and adjacent(source,i)==adjacent(target,i) and
        (a.hydrogens,a.formal_charge,a.nonbonding_electrons,a.radical_electrons)==
        (b.hydrogens,b.formal_charge,b.nonbonding_electrons,b.radical_electrons) and
        a.formal_charge==0 and a.radical_electrons==0 and
        all(n==1 for n in adjacent(source,i).values()))


def reconcile(source,target,event):
    """Restore proven unchanged local geometry; reject unproved remote flips.

    This does not fill newly created, planarised, radical, charged or otherwise
    changed centres. A named SN2 keeps the existing local inversion rule.
    """
    old,new=atoms_by_id(source),atoms_by_id(target);notes=[]
    for i,a in old.items():
        if i not in new or not stable(source,target,i):continue
        p=parity(a.stereo,source,i)
        if p is None:continue
        q=parity(new[i].stereo,target,i)
        if q is None:
            previous=deepcopy(new[i].stereo);new[i].stereo=deepcopy(a.stereo)
            notes.append({'atom_id':i,'effect':'retention','basis':'unchanged mapped sigma neighbourhood and local occupation',
                'replaced_uncertain_annotation':previous,'rule_version':VERSION})
        elif p!=q:raise FormalError('Unproved stereochemical inversion at unchanged atom '+str(i))
    # Untouched E/Z constraints are geometrical relationships, not frozen CIP
    # labels. A changed double bond or its direct substituent set is reassessed.
    for bond in source.bonds:
        if bond.order!=2 or bond.stereo.get('status')!='specified':continue
        a,b=bond.a,bond.b
        after=next((x for x in target.bonds if {x.a,x.b}=={a,b}),None)
        if after is None or after.order!=2 or any(adjacent(source,i)!=adjacent(target,i) for i in (a,b)):continue
        if after.stereo.get('status')!='specified':after.stereo=deepcopy(bond.stereo)
        elif after.stereo!=bond.stereo:
            # Endpoint order alone may reverse the stored references.
            left=after.stereo.get('reference_atoms',[])
            expected=bond.stereo.get('reference_atoms',[])
            if after.a!=bond.a:left=list(reversed(left))
            if left!=expected or after.stereo.get('rdkit_stereo')!=bond.stereo.get('rdkit_stereo'):
                raise FormalError('Unproved E/Z change at unchanged bond '+str((a,b)))
    return notes


def relation(source,target,event,i):
    a,b=atoms_by_id(source)[i],atoms_by_id(target)[i]
    ns,nt=adjacent(source,i),adjacent(target,i)
    if event.event_type=='SN2' and event.bindings.get('electrophile')==i:
        return {'effect':'inversion','replacement':{str(event.bindings['leaving_group']):event.bindings['nucleophile']},
                'basis':'declared SN2 backside substitution; relative parity, not CIP inversion'}
    if stable(source,target,i):return {'effect':'retention','replacement':{},'basis':'unchanged local sigma geometry'}
    if a.element!='C':return {'effect':'unknown','basis':'non-carbon stereodynamics not generalised'}
    before=sum(ns.values())+a.hydrogens;after=sum(nt.values())+b.hydrogens
    if any(v>1 for v in nt.values()) or after<4 or b.radical_electrons:
        return {'effect':'destroyed_or_labile','basis':'unsaturation, low coordination or open-shell centre; no inversion implied'}
    if b.formal_charge:return {'effect':'unknown','basis':'charged centre requires a configurational-stability model'}
    if ns!=nt or before!=after:return {'effect':'created_or_changed','basis':'new local neighbourhood; facial selectivity is separate'}
    return {'effect':'unknown','basis':'occupation or geometry model incomplete'}


def transported(stereo,rel,reverse=False):
    value=deepcopy(stereo);replace={int(a):b for a,b in rel.get('replacement',{}).items()}
    if reverse:replace={b:a for a,b in replace.items()}
    value['neighbor_order']=[replace.get(i,i) for i in value.get('neighbor_order',[])]
    if rel['effect']=='inversion' and value.get('rdkit_tag') in TAGS:
        value['rdkit_tag']='CHI_TETRAHEDRAL_CCW' if value['rdkit_tag']=='CHI_TETRAHEDRAL_CW' else 'CHI_TETRAHEDRAL_CW'
    return value


def analyze_path(document,*,reference_view=None,known_product_diagnostic=False):
    """Report original uncertainty and backwards constraints without mutating.

    A mapped reference is accepted only by the explicit known-product mode.
    Inferred labels are patch proposals, never input to forward generation.
    """
    if reference_view is not None and not known_product_diagnostic:
        raise ValueError('A product reference is restricted to known-product diagnosis')
    views={v.view_id:v for v in document.views};rows=[];uncertain=[];potential={k:potential_centers(v) for k,v in views.items()}
    pi={k:potential_double_bonds(v) for k,v in views.items()}
    for e in document.events:
        if e.reactant_view_id not in views or e.product_view_id not in views:continue
        s,t=views[e.reactant_view_id],views[e.product_view_id];effects=[]
        for i,a in atoms_by_id(s).items():
            if i not in atoms_by_id(t):continue
            b=atoms_by_id(t)[i]
            if a.stereo.get('status')=='specified' or b.stereo.get('status')=='specified' or i in potential[s.view_id]|potential[t.view_id]:
                rel=relation(s,t,e,i)
                if a.element=='C':
                    effects.append({'atom_id':i,**rel,'before':a.stereo,'after':b.stereo})
                    if rel['effect'] not in {'retention','inversion'} or b.stereo.get('status')!='specified':
                        uncertain.append({'event_id':e.event_id,'atom_id':i,'effect':rel['effect'],'reason':rel['basis'],
                            'uncertainty_kind':'geometry_unspecified' if b.stereo.get('status')!='specified' else 'selectivity_not_independently_proven'})
        bonds=[]
        sb={tuple(sorted((b.a,b.b))):b for b in s.bonds};tb={tuple(sorted((b.a,b.b))):b for b in t.bonds}
        for pair in sorted(pi[s.view_id]|pi[t.view_id]):
            a,b=sb.get(pair),tb.get(pair)
            same=bool(a and b and a.order==b.order==2 and all(adjacent(s,i)==adjacent(t,i) for i in pair))
            effect='retention' if same else 'created_or_changed' if b and b.order==2 else 'destroyed_or_labile'
            item={'atom_ids':list(pair),'effect':effect,'before':a.stereo if a else {},'after':b.stereo if b else {},
                'basis':'unchanged mapped pi neighbourhood' if same else 'pi geometry lost or generated; E/Z preference needs evidence'}
            bonds.append(item)
            if effect!='retention' or not b or b.stereo.get('status')!='specified':
                uncertain.append({'event_id':e.event_id,'atom_ids':list(pair),'effect':effect,
                    'reason':item['basis'],'uncertainty_kind':'double_bond_geometry_or_selectivity'})
        rows.append({'event_id':e.event_id,'before':e.reactant_view_id,'after':e.product_view_id,'effects':effects,'bond_effects':bonds})
    proposals=[];conflicts=[]
    if reference_view is not None and document.events:
        end=views[document.events[-1].product_view_id];known=atoms_by_id(end)
        for i,a in atoms_by_id(reference_view).items():
            if i not in known:raise ValueError('Reference atom is absent from executed endpoint')
            if (a.element,a.isotope)!=(known[i].element,known[i].isotope):raise ValueError('Reference identity differs')
            if a.stereo.get('status')!='specified':continue
            if adjacent(reference_view,i)!=adjacent(end,i) or a.hydrogens!=known[i].hydrogens:
                raise ValueError('Reference stereocentre neighbourhood differs at atom '+str(i))
            value=deepcopy(a.stereo)
            for e in reversed(document.events):
                if e.reactant_view_id not in views or e.product_view_id not in views:break
                s,t=views[e.reactant_view_id],views[e.product_view_id]
                if i not in atoms_by_id(t):break
                current=atoms_by_id(t)[i];p=parity(value,t,i);q=parity(current.stereo,t,i)
                if p is not None and q is not None and p!=q:
                    conflicts.append({'event_id':e.event_id,'view_id':t.view_id,'atom_id':i,
                        'reason':'known product conflicts with retained/executed geometry','executed':current.stereo,'reference_constraint':value})
                elif p is not None and q is None:
                    proposals.append({'event_id':e.event_id,'view_id':t.view_id,'atom_id':i,'stereo':deepcopy(value),
                        'basis':'known_product_reverse_constraint','level':'inferred','human_approved':False})
                rel=relation(s,t,e,i)
                if rel['effect'] not in {'retention','inversion'}:
                    proposals.append({'event_id':e.event_id,'atom_id':i,'inference_boundary':rel['effect'],
                        'reason':'cannot propagate geometry across destruction/creation/unknown selectivity','human_approved':False});break
                value=transported(value,rel,reverse=True)
            else:
                initial=views[document.events[0].reactant_view_id];current=atoms_by_id(initial)[i]
                p=parity(value,initial,i);q=parity(current.stereo,initial,i)
                if p is not None and q is not None and p!=q:
                    conflicts.append({'view_id':initial.view_id,'atom_id':i,
                        'reason':'known product conflicts with the initial geometry across only retention/inversion steps',
                        'executed':current.stereo,'reference_constraint':value,
                        'allowed_stereo_change_events':[],'requires_source_review':True})
                elif p is not None and q is None:
                    proposals.append({'view_id':initial.view_id,'atom_id':i,'stereo':deepcopy(value),
                        'basis':'known_product_reverse_constraint','level':'inferred','human_approved':False})
    return {'version':VERSION,'operation':'STEREOCHEMISTRY','status':'completed_with_uncertainty' if uncertain else 'complete',
        'formal_execution_blocked_by_unknown_stereo':False,'events':rows,'uncertain_events':uncertain,
        'reference_mode':'known_product_diagnostic' if reference_view is not None else 'input_and_executed_states_only',
        'reference_used_for_generation':False,'inferred_patch_proposals':proposals,'conflicts':conflicts,
        'authority_mutated':False,'human_approved':False}
