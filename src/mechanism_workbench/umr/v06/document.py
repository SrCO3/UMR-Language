import sys
from ..v02.document import run_document as shared_run
from .ledger import compare_views,signature,validate_view
from .execution import execute_event,make_state,product_spin,prepare_source
from .resonance import resonance_scan,validate_resonance
from .ast import AnalysisRecord
from .stereo import analyze


def run_document(document,*,max_views=16,max_expansions=64):
    result=shared_run(document,max_views=max_views,max_expansions=max_expansions,_backend=sys.modules[__name__])
    existing={a.analysis_id for a in result.document.analyses}
    for state in result.document.states:
        view=next(v for v in result.document.views if v.view_id==state.primary_view_id)
        key='stereo.'+view.view_id
        if key not in existing:
            report=analyze(view)
            result.document.analyses.append(AnalysisRecord(analysis_id=key,operation='STEREOCHEMISTRY',
                state_id=state.state_id,view_id=view.view_id,status='complete' if report['status']=='complete' else 'unsupported',output=report))
    return result
