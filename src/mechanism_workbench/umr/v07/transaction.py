"""Version adapter to the shared pair/fishhook engines; no serial arrow commits."""
from collections import Counter
from ..engine import apply_actions
from ..v02.electrons import apply_fishhooks
from ..chemistry import FormalError,Unsupported,atoms_by_id
from .ledger import validate_view,project_symbols,restore_symbols,symbols
from .elements import charge,coordination


def touched(actions):
    ids=set()
    for a in actions:
        if a.op=='FISHHOOK':
            for e in (a.source,a.destination):ids.update([e.atom_id] if hasattr(e,'atom_id') else [e.a,e.b])
        else:ids.update(getattr(a,k) for k in ('a','b','c') if hasattr(a,k))
    return ids


def transact(source,actions,state_id,view_id,*,validator=validate_view,preserve_stereo=False,after=None,
             charge_fn=charge,coordination_fn=coordination,_unknown_pair_receivers=frozenset(),
             _declared_pair_donors=frozenset()):
    validator(source);active=touched(actions)
    if active&symbols(source):raise Unsupported('Electron actions cannot touch an opaque R/n fragment')
    reserved,_=coordination_fn(project_symbols(source));atoms=atoms_by_id(source)
    debits=Counter(a.a for a in actions if a.op=='PAIR')
    lp_slots=Counter()
    for a in actions:
        if a.op=='FISHHOOK' and a.source.kind=='lone_pair':lp_slots[a.source.atom_id]+=1
    for i in set(debits)|set(lp_slots):
        r=atoms[i].radical_electrons
        if r is None:
            if i not in _declared_pair_donors or lp_slots[i] or debits[i]!=1:
                raise Unsupported('Active donor needs resolved occupation')
            # The caller must validate a finite, explicitly declared paired
            # channel. This does not assert the remaining occupation or spin.
            r=0
        if 2*debits[i]+lp_slots[i]+reserved[i]>atoms[i].nonbonding_electrons-r:
            raise FormalError('Repeated consumption of coordination-reserved electrons')
    shadow=project_symbols(source);unknown={a.atom_id for a in source.atoms if a.radical_electrons is None and a.atom_id not in symbols(source)}
    if active&unknown and any(a.op=='FISHHOOK' for a in actions):raise Unsupported('Single-electron participant has unresolved spin occupation')
    for a in shadow.atoms:
        if a.atom_id in unknown:a.radical_electrons=0
    def restore(v):
        # Product binding reconstruction works on the private fixed-boundary
        # projection. Unknown R/n charges are restored before validation/export.
        if after and symbols(source) and v.state_id==state_id and v.view_id==view_id:
            after(v)
        v=restore_symbols(source,v)
        for a in v.atoms:
            if a.atom_id in unknown:
                change=a.nonbonding_electrons-atoms[a.atom_id].nonbonding_electrons
                paired_only=all(x.op!='FISHHOOK' for x in actions)
                permitted=(a.atom_id in _unknown_pair_receivers and change>0 and change%2==0 or
                           a.atom_id in _declared_pair_donors and change==-2)
                if change and not(permitted and paired_only):
                    raise Unsupported('Unresolved-spin center cannot change nonbonding occupation')
                a.radical_electrons=None
        if after and not symbols(source) and v.state_id==state_id and v.view_id==view_id:after(v)
        return v
    def check(v):return validator(restore(v.model_copy(deep=True)))
    if all(a.op=='FISHHOOK' for a in actions):
        result=apply_fishhooks(shadow,actions,state_id,view_id,preserve_stereo=preserve_stereo,
            _validator=check,_charge=charge_fn,_allow_lone_pairs=True)
    elif any(a.op=='FISHHOOK' for a in actions):raise Unsupported('Mixed pair/fishhook transaction unsupported')
    else:result=apply_actions(shadow,actions,state_id,view_id,preserve_stereo=preserve_stereo,_validator=check,_charge=charge_fn)
    return restore(result)
