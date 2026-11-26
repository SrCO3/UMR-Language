"""Finite bromide opening of a positively charged three-membered C–Br–C ring."""
from ..chemistry import FormalError,atoms_by_id,bond_orders
from ..v06.stereo import apply_event
from .execution import transact
from .ledger import validate_view

def opening(source,event,waivers):
    if waivers or event.event_type!='POLAR' or not event.actions:
        raise FormalError('Bromonium opening requires a polar action group')
    roles=event.bindings
    if set(roles)!={'nucleophile','electrophile','leaving_group'}:raise FormalError('Bromonium roles incomplete')
    n,c,l=(roles[k] for k in ('nucleophile','electrophile','leaving_group'))
    a=atoms_by_id(source);edges=bond_orders(source)
    neighbors=lambda i:{v if u==i else u:o for (u,v),o in edges.items() if i in (u,v)}
    br=neighbors(l)
    other=set(br)-{c}
    if (len({n,c,l})!=3 or a[n].element!='Br' or a[n].formal_charge!=-1 or neighbors(n)
        or a[l].element!='Br' or a[l].formal_charge!=1 or br.get(c)!=1 or len(other)!=1):
        raise FormalError('Expected free bromide and a bridged bromonium ion')
    second=next(iter(other))
    if (any((a[i].element,a[i].formal_charge)!=('C',0) for i in (c,second)) or br[second]!=1 or neighbors(c).get(second) not in (1,2)
        or any(a[i].radical_electrons!=0 for i in (n,c,l,second))):
        raise FormalError('Expected closed-shell three-membered C-Br-C topology')
    actual=[x.model_dump(mode='json') for x in event.actions]
    expected=[{'op':'PAIR','a':n,'b':c},{'op':'SPLIT','a':c,'b':l,'target':l}]
    if len(actual)!=2 or any(x not in actual for x in expected):raise FormalError('Bromonium opening arrows disagree with bound roles')
    target=transact(source,event.actions,event.product_state_id,event.product_view_id)
    return target,{'event_id':event.event_id,'execution_mode':'electron_actions','formal_status':'passed','atomic_commit':True,
       'compiled_actions':actual,'ledger':validate_view(target),'stereochemistry':apply_event(source,target,event),
       'chemical_feasibility':'unknown','source':'https://art.torvergata.it/handle/2108/11752','human_approved':False}

