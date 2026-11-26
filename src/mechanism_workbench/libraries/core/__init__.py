"""Shared infrastructure. Keep chemistry optional until a native service is used."""

__all__ = ['ChemicalState', 'LewisView', 'ReactionEvent', 'MechanismDocument',
           'run_document', 'execute_event', 'resonance_scan']


def __getattr__(name):
    from importlib import import_module
    modules = {'ChemicalState': 'ast', 'LewisView': 'ast', 'ReactionEvent': 'ast',
        'MechanismDocument': 'ast', 'run_document': 'document',
        'execute_event': 'execution', 'resonance_scan': 'resonance'}
    if name not in modules:
        raise AttributeError(name)
    return getattr(import_module('mechanism_workbench.umr.v09.' + modules[name]), name)
