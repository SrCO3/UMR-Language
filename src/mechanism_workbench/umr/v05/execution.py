"""Element adapters around the existing atomic electron transaction engine."""
from collections import Counter
from ..ast import Split,Shift
from ..chemistry import FormalError,Unsupported,atoms_by_id,order,electron_count
from ..engine import apply_actions,check_event_rule,action_key
from ..v04.execution import execute_event as execute04,product_spin as spin04,make_state as state04
from .ast import LewisAdduct
from .elements import charge,METAL_GROUP,MAIN
from .coordination import checks
from .ledger import validate_view

TEMPLATES={
 'BORON_L_ASSOCIATION':{'roles':['donor','boron'],'scope':'one N/O/P/S lone pair to neutral trivalent B; donor-owned convention'},
 'BORON_L_DISSOCIATION':{'roles':['donor','boron'],'scope':'inverse of the explicitly stored B donor reservation'},
 'METAL_X_HETEROLYSIS':{'roles':['metal','leaving_group'],'scope':'single monodentate sigma X bond; electrons move to ligand'},
 'GRIGNARD_CARBONYL_TRANSFER':{'roles':['metal','nucleophile','carbonyl_carbon','oxygen'],'scope':'monomeric RMgX formal concerted carbon transfer to C=O; aggregation/solvation/kinetics unknown'},
}

def execute_event(source,event,*,_validator=validate_view):
    validate_view=_validator
    validate_view(source)
    if event.reactant_view_id!=source.view_id or event.product_state_id==source.state_id:
        raise FormalError('Event requires input view and a new chemical state')
    if event.event_type=='ISC':raise Unsupported('requires_spin_transition: spin conversion not implemented')
    atoms=atoms_by_id(source);reserved,_=checks(source)
    extended=source.metal_centers or source.lewis_adducts or any(a.element in MAIN for a in source.atoms)
    if not event.element_template and (event.reactive_template or event.template_id or event.event_type in {'RADICAL','PERICYCLIC','SET','COORDINATION'}):
        if extended:raise Unsupported('This legacy reactive/radical/pericyclic/metal template has no v0.5 element adapter; retain its original schema version')
        return execute04(source,event)
    remove_bindings=set();adduct=None
    if event.element_template:
        if event.element_template not in TEMPLATES:raise Unsupported('Unregistered element template')
        roles=TEMPLATES[event.element_template]['roles']
        try:ids=[event.bindings[r] for r in roles]
        except KeyError as exc:raise FormalError('Missing element-template role') from exc
        if any(not isinstance(i,int) or i not in atoms for i in ids) or len(set(ids))!=len(ids):raise FormalError('Invalid/distinct role binding required')
        if not event.evidence:raise FormalError('Element template requires source/model-assumption evidence')
        if event.element_template.startswith('BORON_L_'):
            if event.event_type!='COORDINATION':raise FormalError('B association/dissociation uses COORDINATION event classification')
            d,b=ids;expected=[]
            if event.element_template=='BORON_L_ASSOCIATION':
                if any(x.acceptor_atom_id==b for x in source.lewis_adducts):raise FormalError('B acceptor already occupied')
                adduct=LewisAdduct(binding_id='adduct.'+event.event_id,donor_atom_id=d,acceptor_atom_id=b,evidence=event.evidence)
            elif not any(x.donor_atom_id==d and x.acceptor_atom_id==b for x in source.lewis_adducts):raise FormalError('Declared B adduct is absent')
        else:
            m,x=ids[:2]
            if atoms[m].element not in METAL_GROUP or order(source,m,x)!=1:raise FormalError('Expected one supported metal-X sigma bond')
            binding=next((b for b in source.metal_bindings if b.metal_atom_id==m and b.donor_atom_ids==[x] and b.mode=='sigma_covalent'),None)
            if binding is None:raise FormalError('Missing sigma binding')
            ligand=next(l for l in source.ligands if l.ligand_id==binding.ligand_id)
            if ligand.ionic_class!='X':raise Unsupported('Only explicit monodentate X ligand partition registered')
            remove_bindings.add(binding.binding_id)
            if event.element_template=='METAL_X_HETEROLYSIS':
                if event.event_type!='HETEROLYSIS':raise FormalError('Metal X cleavage must be HETEROLYSIS')
                expected=[Split(a=m,b=x,target=x)]
            else:
                c,o=ids[2:]
                if event.event_type!='POLAR' or atoms[m].element!='Mg' or atoms[x].element!='C' or atoms[c].element!='C' or atoms[o].element!='O' or order(source,c,o)!=2:
                    raise FormalError('Grignard carbonyl-transfer input roles/profile mismatch')
                others=[b for b in source.metal_bindings if b.metal_atom_id==m and b.binding_id!=binding.binding_id]
                if len(others)!=1 or len(others[0].donor_atom_ids)!=1 or atoms[others[0].donor_atom_ids[0]].element not in {'Cl','Br','I'}:
                    raise Unsupported('Grignard template limited to monomeric RMgX without aggregation/extra donors')
                expected=[Shift(a=m,b=x,c=c),Split(a=c,b=o,target=o)]
        if event.actions is not None and Counter(action_key(a) for a in event.actions)!=Counter(action_key(a) for a in expected):raise FormalError('Actions disagree with registered element template')
        event.actions=expected
        report={'template_scope':TEMPLATES[event.element_template]['scope'],'mechanism_assumption':'formal localized model; not proof of actual pathway','conditions':event.conditions}
    else:
        if not event.actions:raise FormalError('Polar event requires an action group')
        touched={getattr(a,k) for a in event.actions for k in ('a','b','c') if hasattr(a,k)}
        if any(atoms[i].element in METAL_GROUP for i in touched):raise Unsupported('Metal reaction requires an explicit registered element_template')
        report=check_event_rule(source,event)
    spent=Counter(a.a for a in event.actions if a.op=='PAIR')
    for i,count in spent.items():
        if 2*count+reserved[i]>atoms[i].nonbonding_electrons-(atoms[i].radical_electrons or 0):raise FormalError('Repeated consumption of a coordination-reserved pair')
    def commit_validator(view):
        if view.state_id==event.product_state_id:
            view.metal_bindings=[b for b in view.metal_bindings if b.binding_id not in remove_bindings]
            # Occupied-site geometry becomes unknown after a binding change.
            for m in view.metal_centers:
                if remove_bindings:
                    m.relative_positions=[r for r in m.relative_positions if r.first_binding_id not in remove_bindings and r.second_binding_id not in remove_bindings]
                    m.vacant_site_ids=None
            if event.element_template=='BORON_L_ASSOCIATION':view.lewis_adducts.append(adduct)
            if event.element_template=='BORON_L_DISSOCIATION':view.lewis_adducts=[x for x in view.lewis_adducts if (x.donor_atom_id,x.acceptor_atom_id)!=tuple(ids)]
        return validate_view(view)
    target=apply_actions(source,event.actions,event.product_state_id,event.product_view_id,_validator=commit_validator,_charge=charge)
    report.update(event_id=event.event_id,event_type=event.event_type,element_template=event.element_template,
        rule_version='umr-elements/0.5',atomic_commit=True,formal_status='passed',ledger=validate_view(target),
        compiled_actions=[a.model_dump(mode='json') for a in event.actions],human_approved=False,chemical_feasibility='unknown')
    return target,report

def product_spin(origin,event,target):
    if any(a.radical_electrons is None for a in target.atoms):
        if origin.multiplicity is not None or event.product_multiplicity is not None:raise Unsupported('Unknown metal occupation cannot establish multiplicity')
        return None,{'status':'unknown','spin_transition_executed':False}
    return spin04(origin,event,target)

def make_state(origin,event,target,multiplicity,report):
    state=state04(origin,event,target,multiplicity,report)
    if any(a.radical_electrons is None for a in target.atoms):state.electronic_mode='UNKNOWN'
    return state
