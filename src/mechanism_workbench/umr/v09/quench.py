"""Explicit proposed-completion branches; never rewrite the source trajectory.

Water is a declared material admission, not an implicit source of untracked H.
The finite protonation rules predict candidates, not isolated-product identity.
"""
import hashlib,json
from .ast import MechanismDocument,ReactionEvent,ReactionPhase,ChemicalState
from .smiles import import_mapped_smiles,export_mapped_smiles
from ..ast import Pair,Split
from ..chemistry import atoms_by_id,bond_orders,bond_valences,Unsupported,FormalError

LABEL='推定淬灭'

def propose(result,max_candidates=4):
    base={'label':LABEL,'origin':'proposed_completion','human_approved':False,
          'counts_as_source_mechanism':False,'chemical_feasibility':'not_reviewed',
          'assumption':'原文未给出充分后处理细节；条件化假设最后加入水。',
          'scope':'有限游离负离子质子化或单核 Mg/Zn–C/O 水解；不预测分离收率、pH或全部水解反应。'}
    if not result.execution_complete:return {**base,'status':'not_proposed_incomplete_reaction','candidates':[]}
    phases=result.document.phases
    if any(p.kind in {'workup','quench'} and p.event_ids and p.status=='specified' for p in phases):
        return {**base,'status':'explicit_workup_present','candidates':[]}
    doc=result.document
    viewid=doc.events[-1].product_view_id if doc.events else next(s.primary_view_id for s in doc.states if s.state_id==doc.initial_state_id)
    v=next(v for v in doc.views if v.view_id==viewid);ats=atoms_by_id(v);edges=bond_orders(v);vals=bond_valences(v)
    protected={i for b in v.metal_bindings for i in b.donor_atom_ids}
    special={c.atom_id for c in v.reactive_centers}
    sites=[]
    for i,a in sorted(ats.items()):
        if i in protected|special or a.radical_electrons!=0 or a.formal_charge!=-1:continue
        neighbors=[j for edge in edges if i in edge for j in edge if j!=i]
        if a.element=='C' and vals[i]==3 and a.nonbonding_electrons==2:sites.append((i,None))
        if a.element=='O' and len(neighbors)==1 and ats[neighbors[0]].element=='C' and vals[i]==1:
            # Carboxylate and hetero-oxo salts are not default water-protonation targets.
            carbon=neighbors[0]
            if not any(carbon in edge and o>=2 and ats[next(j for j in edge if j!=carbon)].element in {'O','S','N'} for edge,o in edges.items()):sites.append((i,None))
    for b in v.metal_bindings:
        if b.mode!='sigma_covalent' or len(b.donor_atom_ids)!=1:continue
        i=b.donor_atom_ids[0];m=b.metal_atom_id
        if ats[m].element in {'Mg','Zn'} and ats[i].element in {'C','O'} and edges.get(tuple(sorted((i,m))))==1:sites.append((i,m))
    # Different sites are explicit alternatives, never collapsed into one answer.
    sites=sorted(set(sites),key=lambda p:(p[0],p[1] or 0));rows=[]
    for k,(site,metal) in enumerate(sites[:max_candidates]):
        try:
            sv=v.model_copy(deep=True);sv.state_id='quench.s0';sv.view_id='quench.v0'
            first=max(ats)+1;o,h,h2=first,first+1,first+2
            _,water=import_mapped_smiles(f'[O:{o}]([H:{h}])[H:{h2}]')
            sv.atoms.extend(water.atoms);sv.bonds.extend(water.bonds)
            source_state=next(s for s in doc.states if s.state_id==v.state_id).model_copy(deep=True)
            source_state.state_id=sv.state_id;source_state.primary_view_id=sv.view_id
            evidence=[{'level':'proposed_completion','source':'UMR v0.9 default aqueous-workup policy (user-requested)',
                       'note':LABEL+'；水为新增物料，非原文已给出的试剂。水解方向与分离产物仍需化学审核。'}]
            common=dict(event_id='quench.e1',reactant_view_id=sv.view_id,product_state_id='quench.s1',product_view_id='quench.v1',
                phase_id='quench',evidence=evidence,explanation=LABEL+'：显式水提供 H，保留羟基/金属羟基副产物。',conditions={'quench_medium':'water','conditions_origin':'proposed_completion'})
            if metal is None:
                event=ReactionEvent(**common,event_type='PROTON_TRANSFER',actions=[Pair(a=site,b=h),Split(a=h,b=o,target=o)],bindings={'base':site,'hydrogen':h,'acid':o})
            else:
                event=ReactionEvent(**common,event_type='POLAR',general_event='AQUEOUS_METAL_QUENCH',general_roles={'metal':metal,'organic':site,'water_o':o,'hydrogen':h},mechanism_assumption='Finite aqueous hydrolysis model; inferred workup, not an asserted elementary transition state')
            digest=hashlib.sha256(json.dumps(v.model_dump(mode='json'),sort_keys=True).encode()).hexdigest()
            candidate=MechanismDocument(initial_state_id=sv.state_id,states=[source_state],views=[sv],events=[event],
                phases=[ReactionPhase(phase_id='quench',kind='quench',event_ids=[event.event_id],activate_atom_ids=[o,h,h2],conditions=event.conditions,evidence=evidence)],
                metadata={'human_approved':False,'example_scope':'inferred_quench_branch','display_label':LABEL,
                  'parent_view_id':v.view_id,'parent_view_sha256':digest,'material_admissions':[{'species':'water','atom_ids':[o,h,h2],'origin':'proposed_completion'}],
                  'source_mechanism_verified':False,'alternative_sites':len(sites),'complete_isolated_product_claim':False})
            from .document import run_document
            qr=run_document(candidate,max_views=4,max_expansions=16,propose_quench=False)
            row={'site':site,'metal':metal,'status':qr.status,'errors':qr.errors,'input_document':candidate.model_dump(mode='json'),
                 'result':qr.to_dict(),'material_admission':candidate.metadata['material_admissions'][0]}
            if qr.execution_complete:
                product=next(w for w in qr.document.views if w.view_id=='quench.v1')
                try:row['product_mapped_smiles']=export_mapped_smiles(product)
                except (Unsupported,ValueError):row['product_mapped_smiles']=None
                row['product_view']=product.model_dump(mode='json')
            rows.append(row)
        except (FormalError,Unsupported,ValueError) as exc:rows.append({'site':site,'metal':metal,'status':'unsupported','reason':str(exc)})
    return {**base,'status':'candidates_generated' if rows else 'no_supported_workup_change',
       'candidates':rows,'alternative_sites':len(sites),'truncated':len(sites)>max_candidates,
       'no_change_note':None if rows else '保留末态作为推定淬灭后的候选结构；当前规则没有可执行的水淬灭变化，不宣称真实后处理不会改变结构。',
       'unchanged_product_view':v.model_dump(mode='json') if not rows else None}
