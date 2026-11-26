"""Donor-owned B coordination: no covalent electron mutation or unknown filling."""
from ..chemistry import FormalError, Unsupported, atoms_by_id, identity
from ..v05.ast import LewisAdduct
from ..v07.ledger import symbols
from .ledger import validate_view
from .ast import RULE_VERSION

def execute(source,event):
    validate_view(source)
    if event.event_type!='COORDINATION' or set(event.bindings)!={'donor','boron'}:
        raise FormalError('B coordination requires exact donor/boron roles')
    d,b=[event.bindings[k] for k in ('donor','boron')];ats=atoms_by_id(source)
    if d==b or not {d,b}<=ats.keys():raise FormalError('Distinct existing B coordination endpoints required')
    if {d,b}&symbols(source):raise Unsupported('Opaque R/n cannot be an active B coordination endpoint')
    if ats[b].element!='B':raise FormalError('B acceptor role must be boron')
    if event.actions or not event.evidence:raise FormalError('Donor-owned coordination requires evidence and no covalent arrows')
    target=source.model_copy(deep=True)
    matches=[a for a in target.lewis_adducts if (a.donor_atom_id,a.acceptor_atom_id)==(d,b)]
    if event.element_template=='BORON_L_DISSOCIATION':
        if len(matches)!=1:raise FormalError('Exactly one existing B adduct required')
        target.lewis_adducts.remove(matches[0])
    else:
        if any(a.acceptor_atom_id==b for a in target.lewis_adducts):raise FormalError('B acceptor already occupied')
        target.lewis_adducts.append(LewisAdduct(binding_id='adduct.'+event.event_id,donor_atom_id=d,acceptor_atom_id=b,evidence=event.evidence))
    target.state_id=event.product_state_id;target.view_id=event.product_view_id
    report=validate_view(target)
    if identity(source)!=identity(target) or source.atoms!=target.atoms or source.bonds!=target.bonds:
        raise FormalError('Donor-owned association changed covalent electron inventory')
    return target,{'event_id':event.event_id,'event_type':event.event_type,'element_template':event.element_template,
        'rule_version':RULE_VERSION,'atomic_commit':True,'formal_status':'passed','execution_mode':'electron_actions',
        'compiled_actions':[],'coordination_change':{'donor':d,'boron':b},'ledger':report,
        'electron_ownership_changed':False,'chemical_feasibility':'unknown','human_approved':False}
