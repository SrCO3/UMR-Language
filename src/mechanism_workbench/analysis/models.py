"""Typed analysis record accepted by the published UMR document schema."""
from ..umr.v06.ast import AnalysisRecord as BaseAnalysis


class AnalysisRecord(BaseAnalysis):
    # The envelope is compatible with ANALYZE; operations are namespaced.
    operation: str
    rule_version: str = 'umr-analysis/0.1'
