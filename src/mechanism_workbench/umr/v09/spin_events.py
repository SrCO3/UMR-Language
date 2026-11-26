"""Explicit finite electronic-state models; never an ISC kinetics calculation.

These are real state transitions, not resonance transforms. The carbonyl model
uses a C/O diradical Lewis occupation; it is not a wavefunction of an n-pi* state.
"""
from collections import Counter
import json
from ..chemistry import FormalError,Unsupported,atoms_by_id,order,identity,electron_count
from ..v02.electrons import bond_to_atom
from ..v04.ast import ElectronicState
from ..v07.ledger import project_symbols
from .ledger import validate_view
from .ast import RULE_VERSION

SPECS={
 'PHOTO_CARBONYL_DIRADICAL':{'roles':['carbon','oxygen'],'source':'source-declared triplet C/O diradical Lewis model; not n-pi* orbital simulation'},
 'DIRADICAL_ISC_TO_SINGLET':{'roles':['a','b'],'source':'evidenced two-center triplet-to-open-shell-singlet declaration; rates unknown'},
}

def execute(source,event):
    from .execution import transact
    validate_view(source)
    name=event.general_event;r=event.general_roles;ats=atoms_by_id(source)
    if event.event_type!='ISC' or event.electron_transfers:raise FormalError('Explicit ISC event required')
    if set(r)!=set(SPECS[name]['roles']) or len(set(r.values()))!=len(r) or not set(r.values())<=ats.keys():
        raise FormalError('Exact distinct electronic transition roles required')
    if not event.evidence or not event.mechanism_assumption:raise FormalError('Electronic-state model requires evidence and an explicit assumption')
    # Opaque spectator R groups remain unknown; no statement about their spin.
    v=project_symbols(source);known=atoms_by_id(v)
    if any(a.radical_electrons is None for a in v.atoms):raise Unsupported('Unresolved active electronic occupation')
    if any(a.element not in {'H','C','N','O','F','S','Cl','Br','I'} for a in v.atoms):
        raise Unsupported('Finite organic diradical transition scope')
    if event.electronic_outcome is not None:raise FormalError('Outcome is derived from the registered transition, not overwritten')
    if name=='PHOTO_CARBONYL_DIRADICAL':
        c,o=r['carbon'],r['oxygen']
        if [ats[i].element for i in [c,o]]!=['C','O'] or order(source,c,o)!=2:raise FormalError('Explicit carbonyl required')
        if any(a.radical_electrons for a in v.atoms) or any(ats[i].formal_charge for i in [c,o]):raise Unsupported('Neutral closed-shell carbonyl occupation required')
        if event.activation!='photochemical' and event.conditions.get('activation')!='photochemical':raise FormalError('Photochemical condition required')
        if event.product_multiplicity!=3:raise FormalError('Registered carbonyl model requires a declared triplet')
        actions=[bond_to_atom(c,o,'pi1',0,c),bond_to_atom(c,o,'pi1',1,o)]
        outcome=ElectronicState(character='open_shell',coupling='high_spin',label='declared triplet carbonyl diradical Lewis model',evidence=event.evidence)
        def post(t):t.electronic_state=outcome.model_copy(deep=True)
        target=transact(source,actions,event.product_state_id,event.product_view_id,after=post)
    else:
        a,b=r['a'],r['b']
        radicals={i for i,x in known.items() if x.radical_electrons}
        if radicals!={a,b} or any(ats[i].radical_electrons!=1 for i in [a,b]):raise Unsupported('Exactly two resolved one-electron centers required')
        if source.electronic_state.coupling!='high_spin' or source.electronic_state.character!='open_shell':raise Unsupported('Declared triplet precursor required')
        if event.product_multiplicity!=1:raise FormalError('Explicit open-shell singlet outcome required')
        # Limited to one connected organic diradical, not arbitrary spin pairs.
        seen={a}
        while True:
            nxt=seen|{i for bond in source.bonds if {bond.a,bond.b}&seen for i in (bond.a,bond.b)}
            if nxt==seen:break
            seen=nxt
        if b not in seen:raise Unsupported('Intermolecular spin coupling needs its own model')
        actions=[];target=source.model_copy(deep=True)
        target.state_id=event.product_state_id;target.view_id=event.product_view_id
        target.electronic_state=ElectronicState(character='open_shell_singlet',coupling='unknown',label='declared coupled singlet diradical; no wavefunction computed',evidence=event.evidence)
        validate_view(target)
    key=lambda x:json.dumps(x.model_dump(mode='json'),sort_keys=True)
    if event.actions is not None and Counter(map(key,event.actions))!=Counter(map(key,actions)):raise FormalError('Actions disagree with electronic-state model')
    before,after=project_symbols(source),project_symbols(target)
    if identity(source)!=identity(target) or electron_count(before)!=electron_count(after) or sum(x.formal_charge for x in before.atoms)!=sum(x.formal_charge for x in after.atoms):raise FormalError('Electronic transition ledger changed')
    event.actions=actions
    return target,{'event_id':event.event_id,'event_type':'ISC','general_event':name,'rule_version':RULE_VERSION,
       'execution_mode':'declared_electronic_transition','formal_status':'passed','atomic_commit':True,
       'compiled_actions':[x.model_dump(mode='json') for x in actions],
       'spin_transition_executed':True,'spin_dynamics_calculated':False,'chemical_feasibility':'unknown',
       'guarantee':'local occupation, references and conservation only; electronic state is source/model declared',
       'ledger':validate_view(target)}

def product_spin(origin,event,target):
    if event.general_event=='PHOTO_CARBONYL_DIRADICAL' and origin.multiplicity not in {None,1}:
        raise FormalError('Carbonyl precursor multiplicity must be singlet or unknown')
    if event.general_event=='DIRADICAL_ISC_TO_SINGLET' and origin.multiplicity!=3:
        raise FormalError('Triplet-to-singlet transition needs a triplet source state')
    return event.product_multiplicity,{'status':'declared_electronic_transition','spin_transition_executed':True,
        'spin_dynamics_calculated':False,'global_spin_verified':False}
