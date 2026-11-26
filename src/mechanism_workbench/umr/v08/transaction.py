from ..v07.transaction import transact as shared
from .elements import charge,coordination
from .ledger import validate_view


def transact(source,actions,state_id,view_id,**kwargs):
    return shared(source,actions,state_id,view_id,validator=validate_view,charge_fn=charge,
                  coordination_fn=coordination,**kwargs)
