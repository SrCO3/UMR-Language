"""Evidence-labelled two-event candidates; never silently rewrite a source event."""
from ..ast import Pair,Split,Evidence
from ..chemistry import FormalError,Unsupported,atoms_by_id,order,bond_valences
from .ast import ReactionEvent,MechanismDocument
from .smiles import import_mapped_smiles

FAMILIES={
 'ACYL':{'pi_element':'O','order':2,'source':'https://openstax.org/books/organic-chemistry/pages/21-2-nucleophilic-acyl-substitution-reactions'},
 'IMIDOYL':{'pi_element':'N','order':2,'source':'UMR authored imidoyl addition/elimination candidate; original-case and primary mechanistic evidence pending'},
 'CYANOGEN_HALIDE':{'pi_element':'N','order':3,'source':'UMR authored halogen-C#N addition/elimination candidate; original-case and primary mechanistic evidence pending'},
 'SNAR':{'pi_element':None,'order':2,'source':'https://doi.org/10.1038/s41557-018-0079-7'},
 'VINYL_ACTIVATED':{'pi_element':'C','order':2,'source':'UMR authored activated-vinyl addition/elimination candidate; substrate-specific primary evidence pending'},
}


def aromatic_ring(view,c,partner=None):
    atoms=atoms_by_id(view);found=[]
    def visit(path):
        if len(path)==6:
            if order(view,path[-1],path[0])==1:found.append(path)
            return
        for b in view.bonds:
            if path[-1] not in {b.a,b.b} or b.order!=(2 if len(path)%2 else 1):continue
            j=b.b if b.a==path[-1] else b.a
            if j not in path and atoms[j].element in {'C','N'} and atoms[j].radical_electrons==0:
                visit(path+[j])
    visit([c])
    return next((p for p in found if partner is None or p[1]==partner),None)


def require_roles(view,roles):
    if set(roles)!={'nucleophile','carbon','pi_partner','leaving'} or len(set(roles.values()))!=4:
        raise FormalError('Substitution requires four distinct explicit roles')
    if not set(roles.values())<=atoms_by_id(view).keys():raise FormalError('Missing substitution role atom')


def check_input(view,family,roles,conditions):
    require_roles(view,roles)
    if family not in FAMILIES:raise Unsupported('Unregistered unsaturated substitution family')
    nu,c,x,lg=[roles[k] for k in ('nucleophile','carbon','pi_partner','leaving')];a=atoms_by_id(view);spec=FAMILIES[family]
    if a[c].element!='C' or order(view,c,x)!=spec['order'] or order(view,c,lg)!=1 or order(view,nu,c):
        raise FormalError('Substitution center / pi partner / leaving bond mismatch')
    if spec['pi_element'] and a[x].element!=spec['pi_element']:raise FormalError('Wrong pi-partner element')
    if a[nu].radical_electrons!=0 or a[nu].nonbonding_electrons is None or a[nu].nonbonding_electrons<2:
        raise FormalError('Paired nucleophile required')
    if a[lg].element not in ({'Cl','Br','I'} if family=='CYANOGEN_HALIDE' else {'F','Cl','Br','I','O','S'}):
        raise Unsupported('Leaving group outside registered scope')
    if family=='CYANOGEN_HALIDE' and (bond_valences(view)[x]!=3 or a[x].formal_charge!=0):raise Unsupported('Only halogen-C#N cyanation scope')
    if family=='SNAR' and not aromatic_ring(view,c,x):raise Unsupported('SNAr scope is alternating six-membered C/pyridine-N ring')
    if family=='VINYL_ACTIVATED' and not conditions.get('activation_evidence'):
        raise Unsupported('Vinylic substitution needs explicit substrate/activation evidence')


def compile_stage(view,event):
    family,stage=event.rule_template.split('.',1);roles=event.template_roles
    require_roles(view,roles);nu,c,x,lg=[roles[k] for k in ('nucleophile','carbon','pi_partner','leaving')]
    if stage=='addition':
        check_input(view,family,roles,event.conditions)
        actions=[Pair(a=nu,b=c),Split(a=c,b=x,target=x)]
    elif stage=='elimination':
        if family not in FAMILIES:raise Unsupported('Unregistered substitution family')
        if order(view,c,x)!=FAMILIES[family]['order']-1 or order(view,c,nu)!=1 or order(view,c,lg)!=1:
            raise FormalError('Elimination does not match the declared addition-state bonds')
        actions=[Pair(a=x,b=c),Split(a=c,b=lg,target=lg)]
    elif stage=='concerted' and family=='SNAR':
        check_input(view,family,roles,event.conditions)
        if event.mechanism_assumption!='concerted' or not event.evidence:raise Unsupported('Concerted SNAr requires explicit evidence and assumption')
        actions=[Pair(a=nu,b=c),Split(a=c,b=lg,target=lg)]
    else:raise Unsupported('Unregistered substitution stage')
    return actions,{'scope':family,'candidate_expansion':event.expansion,'mechanism_confirmed':False}


def expand(document,event_id,*,family,roles,explanation,source):
    """Returns a NEW candidate branch. No original event/approval is overwritten."""
    data=document.model_dump(mode='json');index=next(i for i,e in enumerate(document.events) if e.event_id==event_id)
    original=document.events[index];view=next(v for v in document.views if v.view_id==original.reactant_view_id)
    if original.mechanism_assumption=='concerted':raise ValueError('Explicit concerted source cannot be automatically split')
    if original.exception_ids:raise ValueError('Remove/review existing event exceptions in a separate revision before expansion')
    check_input(view,family,roles,original.conditions)
    mid_s=f'{original.product_state_id}.addition';mid_v=f'{original.product_view_id}.addition'
    if any(s.state_id==mid_s for s in document.states) or any(v.view_id==mid_v for v in document.views):raise ValueError('Candidate IDs already exist')
    evidence=Evidence(level='proposed_completion',source=source,note=explanation)
    common={'event_type':'POLAR','template_roles':roles,'conditions':original.conditions,
        'mechanism_assumption':'stepwise_candidate','mechanism_evidence':'assumed_by_template',
        'explanation':explanation,'evidence':[evidence],
        'expansion':{'original_event_id':event_id,'source_event':original.model_dump(mode='json'),
            'status':'proposed_completion','intermediate_existence':'unverified','source_steps':1,'candidate_steps':2}}
    first=ReactionEvent(event_id=event_id+'.addition',reactant_view_id=original.reactant_view_id,
        product_state_id=mid_s,product_view_id=mid_v,rule_template=family+'.addition',depends_on_event=original.depends_on_event,**common)
    second=ReactionEvent(event_id=event_id,reactant_view_id=mid_v,product_state_id=original.product_state_id,
        product_view_id=original.product_view_id,rule_template=family+'.elimination',depends_on_event=first.event_id,**common)
    data['events'][index:index+1]=[first.model_dump(mode='json'),second.model_dump(mode='json')]
    data['record_id']=document.record_id+'.expanded';data['revision']=1;data['analyses']=[]
    data['metadata']={**data['metadata'],'parent_record_id':document.record_id,'human_approved':False,
        'branch_kind':'correction','source_event_preserved':original.model_dump(mode='json')}
    return MechanismDocument.model_validate(data)


def example(family,smiles,roles,*,conditions=None):
    s,v=import_mapped_smiles(smiles)
    e=ReactionEvent(event_id='e1',event_type='POLAR',reactant_view_id='v0',product_state_id='s2',product_view_id='v2',
        conditions=conditions or {})
    doc=MechanismDocument(record_id=family.lower(),initial_state_id='s0',states=[s],views=[v],events=[e])
    return expand(doc,'e1',family=family,roles=roles,explanation='形式两步候选；来源与证据状态见source，不证明加成态存在、寿命或优势机理。',source=FAMILIES[family]['source'])
