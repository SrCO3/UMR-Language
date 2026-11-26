"""Version adapter to the shared v0.2 document orchestrator; no parallel executor."""
import sys

from ..v02.document import run_document as shared_run
from .execution import execute_event, make_state, product_spin
from .ledger import compare_views, signature, validate_view
from .resonance import resonance_scan, validate_resonance


def run_document(document,*,max_views=32,max_expansions=128):
    result=shared_run(document,max_views=max_views,max_expansions=max_expansions,_backend=sys.modules[__name__])
    # No ISC transition has been executed. Requests remain visible in event/error records.
    result.document.metadata["execution_counters"]={"real_reaction_steps":len(result.events),
        "resonance_derivations":len(result.document.resonance_transforms),"spin_transitions":0,
        "requested_spin_transitions":sum(e.event_type=="ISC" for e in document.events)}
    return result
