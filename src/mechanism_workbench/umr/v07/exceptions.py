"""Append-only human decisions, separate from model-supplied mechanism data."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4
from .ast import RULE_VERSION


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,ensure_ascii=False,allow_nan=False,separators=(',',':')).encode()).hexdigest()


def document_digest(document):
    return digest(document.model_dump(mode='json'))


class ExceptionStore:
    """Local single-writer trust boundary; approval is only invoked by human UI/CLI.

    Not a cryptographic multiuser authentication system. Payload flags are never
    consulted. Filesystem access to the review store is a trusted local capability.
    """
    def __init__(self,root,*,rule_version=RULE_VERSION,runner=None):
        self.root=Path(root); self.rule_version=rule_version; self.runner=runner

    def _rows(self):
        return [json.loads(p.read_text(encoding='utf-8')) for p in sorted(self.root.glob('*.json'))] if self.root.exists() else []

    def decisions(self,document):
        key=document_digest(document); latest={}
        for row in self._rows():
            if row['record_id']==document.record_id:latest[row['exception_id']]=row
        result={}
        for spec in document.evidence_exceptions:
            row=latest.get(spec.exception_id)
            state='not_approved' if row is None else row['decision']
            if row and (row['document_digest']!=key or row['revision']!=document.revision or
                        row['exception_digest']!=digest(spec.model_dump(mode='json')) or
                        row['rule_version']!=self.rule_version or spec.rule_version!=self.rule_version):state='invalidated'
            result[spec.exception_id]={'status':state,'decision':row}
        return result

    def approve(self,document,exception_id,*,reviewer,reason,human_confirm=False,record_binding=None):
        if human_confirm is not True or not reviewer.strip() or not reason.strip():
            raise ValueError('Explicit human confirmation, reviewer and reason are required')
        spec=next((x for x in document.evidence_exceptions if x.exception_id==exception_id),None)
        if spec is None or spec.rule_version!=self.rule_version:raise ValueError('Missing exception or outdated rule version')
        from .document import run_document
        strict=(self.runner or run_document)(document)
        event_report=next((x for x in strict.events if x['event_id']==spec.event_id),None)
        row={'schema_version':'umr-exception-decision/0.1','record_id':document.record_id,'revision':document.revision,
             'document_digest':document_digest(document),'exception_id':exception_id,
             'exception_digest':digest(spec.model_dump(mode='json')),'rule_version':self.rule_version,
             'target_digest':digest({'state':spec.target_state.model_dump(mode='json') if spec.target_state else None,
                                     'view':spec.target_view.model_dump(mode='json') if spec.target_view else None}),
             'decision':'approved','actor_kind':'human','reviewer':reviewer,'reason':reason,
             'record_binding':record_binding,'strict_before_decision':{'status':strict.status,'errors':strict.errors,
                 'event_report':event_report,'failed_event_report':strict.failed_event_report}}
        return self._append(row)

    def revoke(self,document,exception_id,*,reviewer,reason,human_confirm=False):
        if not human_confirm or not reviewer.strip() or not reason.strip():raise ValueError('Human confirmation required')
        row=self.decisions(document).get(exception_id,{}).get('decision')
        if not row:raise ValueError('No decision to revoke')
        return self._append({**row,'decision':'revoked','reviewer':reviewer,'reason':reason})

    def _append(self,row):
        row={**row,'at':datetime.now(timezone.utc).isoformat(),'decision_id':str(uuid4())}
        self.root.mkdir(parents=True,exist_ok=True)
        filename=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')+'-'+row['decision_id']+'.json'
        with (self.root/filename).open('x',encoding='utf-8') as f:json.dump(row,f,ensure_ascii=False,indent=2)
        return row
