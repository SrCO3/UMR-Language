"""Availability/condition chronology. Stored future inventory is not yet mixed."""
from ..chemistry import FormalError,Unsupported
from ..v07.transaction import touched

class PhasePlan:
    def __init__(self,doc):
        self.phases=doc.phases;self.by_event={};self.ids={a.atom_id for v in doc.views for a in v.atoms};self.activation={}
        if len({p.phase_id for p in self.phases})!=len(self.phases):raise FormalError('Duplicate phase ID')
        for k,p in enumerate(self.phases):
            for i in p.activate_atom_ids:
                if i in self.activation or i not in self.ids:raise FormalError('Duplicate/absent staged inventory atom')
                self.activation[i]=k
            for e in p.event_ids:
                if e in self.by_event:raise FormalError('Event belongs to more than one phase')
                self.by_event[e]=(k,p)
        if self.phases:
            if set(self.by_event)!={e.event_id for e in doc.events}:raise FormalError('Phases must assign every event exactly once')
            order=[self.by_event[e.event_id][0] for e in doc.events]
            if order!=sorted(order):raise FormalError('Event sequence violates condition/stage order')
            for e in doc.events:
                if e.phase_id!=self.by_event[e.event_id][1].phase_id:raise FormalError('Event/phase references disagree')

    def before(self,source,event):
        if not self.phases:return event
        k,p=self.by_event[event.event_id]
        if p.status=='missing_details':raise Unsupported('Phase conditions/material incomplete; no automatic water or reagent insertion')
        active=touched(event.actions or [])|set(event.general_roles.values())
        for value in (event.bindings or {}).values():
            if type(value)==int:active.add(value)
            elif isinstance(value,list):active.update(x for x in value if type(x)==int)
        for tx in event.electron_transfers:active.update([tx.donor_atom_id,tx.acceptor_atom_id])
        future={i for i,n in self.activation.items() if n>k}
        if active&future:raise FormalError('Event consumes reagent before its declared introduction phase')
        for b in source.bonds:
            if bool(b.a in future)!=bool(b.b in future):raise FormalError('Future inventory already bonded to current reaction material')
        conflicts={key for key in p.conditions.keys()&event.conditions.keys() if p.conditions[key]!=event.conditions[key]}
        if conflicts:raise FormalError('Phase/event condition conflict: '+','.join(sorted(conflicts)))
        trial=event.model_copy(deep=True);trial.conditions={**p.conditions,**event.conditions}
        return trial

    def report(self,reports):
        executed={r['event_id'] for r in reports};rows=[]
        for p in self.phases:
            complete=bool(p.event_ids) and set(p.event_ids)<=executed and p.status=='specified'
            rows.append({'phase_id':p.phase_id,'kind':p.kind,'status':'completed' if complete else 'pending_workup' if p.kind in {'workup','quench'} else 'incomplete_phase',
             'remaining_events':[i for i in p.event_ids if i not in executed],'automatic_water_added':False})
        return rows
