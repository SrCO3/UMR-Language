"""Small recursive reader: JSON blocks + typed arrow grammar, never eval."""

import json
import re

from .ast import (AnalysisRecord, ChemicalState, LewisView, MechanismDocument,
                  Pair, ReactionEvent, ResonanceTransform, Shift, Split)

IDENT = r"[A-Za-z_][A-Za-z0-9_.:~\-]*"


def _no_constant(value):
    raise ValueError(f"Non-finite JSON constant is forbidden: {value}")


class Parser:
    def __init__(self, text: str):
        self.text, self.pos = text, 0
        self.decoder = json.JSONDecoder(parse_constant=_no_constant)

    def skip(self):
        while True:
            match = re.match(r"(?:\s+|#[^\n]*(?:\n|$))", self.text[self.pos:])
            if not match:
                return
            self.pos += match.end()

    def expect(self, pattern):
        self.skip()
        match = re.match(pattern, self.text[self.pos:])
        if not match:
            line = self.text[:self.pos].count("\n") + 1
            raise ValueError(f"DSL syntax error at line {line}, offset {self.pos}: expected {pattern}")
        self.pos += match.end()
        return match

    def json_object(self):
        self.skip()
        value, end = self.decoder.raw_decode(self.text, self.pos)
        if not isinstance(value, dict):
            raise ValueError("A declaration body must be a JSON object")
        self.pos = end
        return value

    def fishhook_endpoint(self, source):
        self.skip()
        if self.version != "1.0" or self.text[self.pos:self.pos+1] == "{":
            return self.json_object()
        kind = self.expect(r"(RAD|BOND)\b" if source else r"(NEW|ATOM)\b").group()
        self.expect(r"\(")
        first = int(self.expect(r"[1-9][0-9]*").group())
        if kind in {"RAD", "ATOM"}:
            endpoint = {"kind": "radical", "atom_id": first, "slot": 0} if source else {"kind": "atom", "atom_id": first}
        else:
            self.expect(",")
            second = int(self.expect(r"[1-9][0-9]*").group())
            endpoint = {"kind": "bond", "a": first, "b": second, "pair": "sigma"}
            if source:
                self.expect(",")
                endpoint["slot"] = int(self.expect(r"[01]").group())
        self.expect(r"\)")
        return endpoint

    def actions(self):
        self.expect(r"DO\b")
        self.expect(r"\{")
        values = []
        while True:
            self.skip()
            if self.text[self.pos:self.pos+1] == "}":
                self.pos += 1
                return values
            op = self.expect(r"(PAIR|SPLIT|SHIFT|FISHHOOK)\b").group(1)
            self.expect(r"\(")
            if op == "FISHHOOK":
                from .v02.ast import Fishhook
                if self.version in {'0.7','0.8','0.9','1.0'}:
                    from .v07.ast import Fishhook
                source = self.fishhook_endpoint(source=True)
                self.expect(",")
                destination = self.fishhook_endpoint(source=False)
                self.expect(r"\)"); self.expect(";")
                values.append(Fishhook(source=source,destination=destination))
                continue
            a = int(self.expect(r"[1-9][0-9]*").group())
            self.expect(",")
            b = int(self.expect(r"[1-9][0-9]*").group())
            if op == "PAIR":
                value = Pair(a=a, b=b)
            elif op == "SPLIT":
                self.expect("->")
                target = int(self.expect(r"[1-9][0-9]*").group())
                value = Split(a=a, b=b, target=target)
            else:
                self.expect(",")
                c = int(self.expect(r"[1-9][0-9]*").group())
                value = Shift(a=a, b=b, c=c)
            self.expect(r"\)")
            self.expect(";")
            values.append(value)

    def parse(self) -> MechanismDocument:
        self.expect(r"UMR\b")
        version = self.expect(r"(?:0\.[123456789]|1\.0)\b").group()
        self.version=version
        wire_version=version
        if version=="1.0":version="0.9"
        self.expect(r"INITIAL\b")
        initial = self.expect(IDENT).group()
        self.expect(";")
        sections = {"states": [], "views": [], "events": [], "resonance_transforms": [], "analyses": []}
        classes = {"STATE": ("states", "state_id", ChemicalState), "VIEW": ("views", "view_id", LewisView),
                   "REACT": ("events", "event_id", ReactionEvent),
                   "RESONANCE": ("resonance_transforms", "transform_id", ResonanceTransform),
                   "ANALYZE": ("analyses", "analysis_id", AnalysisRecord)}
        document_class=MechanismDocument
        if version=="0.2":
            from .v02 import ast as ast02
            classes.update({"STATE":("states","state_id",ast02.ChemicalState),
                            "REACT":("events","event_id",ast02.ReactionEvent),
                            "RESONANCE":("resonance_transforms","transform_id",ast02.ResonanceTransform),
                            "ANALYZE":("analyses","analysis_id",ast02.AnalysisRecord)})
            document_class=ast02.MechanismDocument
        if version in {"0.3","0.4","0.5","0.6","0.7","0.8","0.9"}:
            if version=='0.9':
                from .v09 import ast as ast03
            elif version=='0.8':
                from .v08 import ast as ast03
            elif version in {'0.7','0.8','0.9'}:
                from .v07 import ast as ast03
            elif version=="0.6":
                from .v06 import ast as ast03
            elif version=="0.5":
                from .v05 import ast as ast03
            elif version=="0.4":
                from .v04 import ast as ast03
            else:
                from .v03 import ast as ast03
            classes.update({"STATE":("states","state_id",ast03.ChemicalState),
                            "VIEW":("views","view_id",ast03.LewisView),
                            "REACT":("events","event_id",ast03.ReactionEvent),
                            "RESONANCE":("resonance_transforms","transform_id",ast03.ResonanceTransform),
                            "ANALYZE":("analyses","analysis_id",ast03.AnalysisRecord)})
            document_class=ast03.MechanismDocument
        metadata = {};options={}
        if version in {'0.7','0.8','0.9'}:
            classes['EXCEPTION']=('evidence_exceptions','exception_id',ast03.EvidenceException)
            sections['evidence_exceptions']=[]
        if version=='0.9':
            classes['CYCLE']=('cycles','cycle_id',ast03.CycleRecord)
            sections['cycles']=[]
            classes['PHASE']=('phases','phase_id',ast03.ReactionPhase)
            sections['phases']=[]
            classes['REAGENT']=('reagents','reagent_id',ast03.ReagentDeclaration)
            sections['reagents']=[]
            for kind, section, cls in [('VALENCE','valence_rules',ast03.ValenceRule),
                                       ('INTERMEDIATE','intermediate_rules',ast03.IntermediateRule)]:
                classes[kind]=(section,'rule_id',cls)
                sections[section]=[]
        if wire_version=="1.0":
            from .v10.ast import LewisView as View10, MechanismDocument as Document10, MethylenationObservation
            document_class=Document10
            classes["VIEW"]=("views","view_id",View10)
            classes["OBSERVE"]=("observations","observation_id",MethylenationObservation)
            sections["observations"]=[]
        while True:
            self.skip()
            if self.pos == len(self.text):
                break
            kind = self.expect(r"(STATE|VIEW|REACT|RESONANCE|ANALYZE|META|EXCEPTION|RECORD|CYCLE|PHASE|VALENCE|INTERMEDIATE|REAGENT|OBSERVE)\b").group(1)
            if kind == "META":
                metadata.update(self.json_object())
            elif kind=='RECORD':
                if version not in {'0.7','0.8','0.9'} or options:raise ValueError('One RECORD declaration is supported in UMR 0.7 only')
                options=self.json_object()
                if set(options)-{'record_id','revision','symbolic_mode'}:raise ValueError('Unknown RECORD fields')
            else:
                if kind not in classes:raise ValueError(f'{kind} is unavailable in UMR {version}')
                if kind=='EXCEPTION' and version not in {'0.7','0.8','0.9'}:raise ValueError('EXCEPTION requires UMR 0.7')
                name = self.expect(IDENT).group()
                data = self.json_object()
                section, key, cls = classes[kind]
                if key in data or "actions" in data:
                    raise ValueError("IDs belong in declaration headers; arrows belong in DO blocks")
                data[key] = name
                if kind in {"REACT", "RESONANCE"}:
                    self.skip()
                    if version in {"0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"} and kind=="REACT" and self.text[self.pos:self.pos+1]==";":
                        data["actions"]=None
                    else:
                        data["actions"] = self.actions()
                sections[section].append(cls.model_validate(data))
            self.expect(";")
        return document_class(initial_state_id=initial, metadata=metadata, **sections,**options)


def parse_dsl(text: str) -> MechanismDocument:
    return Parser(text).parse()


def action_text(action) -> str:
    if action.op=="FISHHOOK":
        return "FISHHOOK("+json.dumps(action.source.model_dump(),ensure_ascii=False)+","+json.dumps(action.destination.model_dump(),ensure_ascii=False)+")"
    if isinstance(action, Pair):
        return f"PAIR({action.a},{action.b})"
    if isinstance(action, Split):
        return f"SPLIT({action.a},{action.b}->{action.target})"
    return f"SHIFT({action.a},{action.b},{action.c})"


def dump_dsl(document: MechanismDocument) -> str:
    version=document.schema_version.split("/")[1]
    lines = [f"UMR {version} INITIAL {document.initial_state_id};"]
    if version=="1.0":
        version="0.9"
        for obs in document.observations:
            body=obs.model_dump(mode="json");name=body.pop("observation_id")
            lines.append("OBSERVE "+name+" "+json.dumps(body,ensure_ascii=False)+";")
    if version in {'0.7','0.8','0.9'}:
        lines.append('RECORD '+json.dumps({k:getattr(document,k) for k in ('record_id','revision','symbolic_mode')},ensure_ascii=False)+';')
        for item in document.evidence_exceptions:
            body=item.model_dump(mode='json');name=body.pop('exception_id')
            lines.append('EXCEPTION '+name+' '+json.dumps(body,ensure_ascii=False,allow_nan=False)+';')
    for kind, key, sequence in (("STATE", "state_id", document.states), ("VIEW", "view_id", document.views),
                                ("RESONANCE", "transform_id", document.resonance_transforms),
                                ("REACT", "event_id", document.events), ("ANALYZE", "analysis_id", document.analyses)):
        for item in sequence:
            body = item.model_dump(mode="json")
            name = body.pop(key)
            actions = body.pop("actions", None)
            declaration = f"{kind} {name} " + json.dumps(body, ensure_ascii=False, allow_nan=False)
            if actions is not None:
                declaration += " DO { " + " ".join(action_text(x) + ";" for x in item.actions) + " }"
            lines.append(declaration + ";")
    if version=='0.9':
        for kind, key_field, rules in [('VALENCE','rule_id',document.valence_rules),
                                       ('INTERMEDIATE','rule_id',document.intermediate_rules),
                                       ('REAGENT','reagent_id',document.reagents)]:
            for rule in rules:
                data=rule.model_dump(mode='json');key=data.pop(key_field)
                lines.append(f'{kind} {key} '+json.dumps(data,ensure_ascii=False,allow_nan=False)+';')
        for phase in document.phases:
            data=phase.model_dump(mode='json');key=data.pop('phase_id')
            lines.append('PHASE '+key+' '+json.dumps(data,ensure_ascii=False)+';')
        for cycle in document.cycles:
            data=cycle.model_dump(mode='json');key=data.pop('cycle_id')
            lines.append('CYCLE '+key+' '+json.dumps(data,ensure_ascii=False)+';')
    lines.append("META " + json.dumps(document.metadata, ensure_ascii=False, allow_nan=False) + ";")
    return "\n".join(lines) + "\n"
