import sys
from ..v02.document import run_document as shared_run
from .execution import execute_event,make_state,product_spin
from .ledger import compare_views,signature,validate_view
from .resonance import resonance_scan,validate_resonance

def run_document(document,*,max_views=16,max_expansions=64):
    return shared_run(document,max_views=max_views,max_expansions=max_expansions,_backend=sys.modules[__name__])
