"""Typeset declarations attached to the views actually represented in a figure."""
import re
from .svg import element, label, text_box
from .palette import INK, PURPLE, TEAL_TEXT
from ..umr.v09.local_rules import ValenceRule, active_rules, active_reagents
from ..umr.dsl import action_text


_NAMES = {'element':'element', 'bond_order_sum':'ΣBO', 'formal_charge':'q',
          'nonbonding_electrons':'N', 'radical_electrons':'r',
          'oxidation_state':'OS', 'occupied_sites':'occupied sites'}


def atom_line(atom_id, fields):
    return f'#{atom_id}   ' + '   '.join(f'{_NAMES[k]}={v if v is not None else "?"}' for k,v in fields.items())


def rule_lines(rule):
    if isinstance(rule, ValenceRule):
        values = {key:getattr(rule,key) for key in _NAMES
                  if key not in {'oxidation_state','occupied_sites'} or getattr(rule,key) is not None}
        yield atom_line(rule.atom_id, values)
    else:
        for test in rule.atoms:
            yield atom_line(test.atom_id,test.model_dump(exclude={'atom_id'}, exclude_none=True))
        if rule.bonds:
            yield '   '.join(f'BO(#{b.a}, #{b.b})={b.order}' for b in rule.bonds)
        if rule.total_charge is not None:yield f'Q = {rule.total_charge}'
    if rule.source:yield rule.source


def wrap(text, width, size):
    current = ''
    for token in re.findall(r'[A-Za-z0-9_./:#=+\-]+|\s+|.', text):
        if current and text_box(current + token,0,0,size)[2] > width:
            yield current.rstrip()
            current = ''
        for char in token:
            if current and text_box(current+char,0,0,size)[2] > width:
                yield current.rstrip()
                current = ''
            if current or not char.isspace():current += char
    if current:yield current.rstrip()


def append_rules(root, view_ids, events=()):
    views = list(dict.fromkeys(view_ids))
    rules = [r for vid in views for r in active_rules(vid)]
    blocks=[]
    for vid in views:
        for reagent in active_reagents(vid):
            field = 'Atoms ' + ', '.join(f'#{i}' for i in reagent.atom_ids) + '   ·   ' + reagent.role
            if reagent.phase_id: field += f'   ·   PHASE {reagent.phase_id}'
            fields = [field] + ([reagent.source] if reagent.source else [])
            blocks.append(('data_reagent_id',reagent.reagent_id,vid,f'REAGENT {reagent.reagent_id}  ·  {reagent.description}',fields))
    for event in events:
        if getattr(event,'general_event',None)!='AUTHORED_ELECTRON_ACTIONS':continue
        scope=f'{event.reactant_view_id} → {event.product_view_id}'
        fields=[v for v in (event.mechanism_assumption,event.explanation) if v]
        fields.extend(action_text(a) for a in event.actions or [])
        fields.extend(f'1e: #{t.donor_atom_id} → #{t.acceptor_atom_id}' for t in event.electron_transfers)
        blocks.append(('data_rule_id',event.event_id,scope,f'REACT {event.event_id} · {event.event_type}',fields))
    for rule in rules:
        kind='VALENCE' if isinstance(rule,ValenceRule) else 'INTERMEDIATE'
        blocks.append(('data_rule_id',rule.rule_id,rule.view_id,f'{kind} {rule.rule_id}  ·  {rule.description}',rule_lines(rule)))
    if not blocks:return []
    _,_,width,old_height = map(float,root.get('viewBox').split())
    scale = float(root.get('height'))/old_height
    section=element('g',id='local-rule-legend')
    root.append(section)
    y=old_height+26
    label(section,'自定义声明 / Author declarations',22,y,size=16,fill=PURPLE,font_weight='bold')
    for field,key,scope,title,fields in blocks:
        group=element('g',data_view_id=scope,**{field:key})
        section.append(group)
        for row in wrap(f'{scope}  ·  {title}',width-58,14):
            y+=21
            label(group,row,23,y,size=14,fill=TEAL_TEXT,font_weight='bold')
        for field in fields:
            for row in wrap(field,width-68,13):
                y+=18
                label(group,row,35,y,size=13,fill=INK)
        y+=6
    height=y+12
    root.set('viewBox',f'0 0 {width:g} {height:g}')
    root.set('height',str(round(height*scale)))
    return [r.rule_id for r in rules]
