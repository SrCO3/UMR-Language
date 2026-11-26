"""Export an already-executed result; never re-run it without its approval context."""
import json
import hashlib,re
from pathlib import Path
from html import escape
from xml.etree import ElementTree as ET
from ..exporting import write_json
from ..dsl import dump_dsl
from .smiles import export_mapped_smiles,sidecar
from .exceptions import document_digest


def render_result(result,directory,config=None):
    from ..v09.local_rules import use_rules
    with use_rules(result.document):
        return _render_result(result,directory,config)

def _render_result(result,directory,config):
    from ...drawing.models import RenderConfig
    from ...drawing.render import draw_panel
    from ...drawing.svg import svg_to_png
    from ...drawing.materials import phase_inventory
    from ...drawing.hydrogens import document_hydrogens
    config=config or RenderConfig(max_relayouts=0,routing_seconds=0.4,max_resonance_panels=2)
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    doc=result.document;views={v.view_id:v for v in doc.views};panels=[]
    reports={r['event_id']:r for r in result.events}
    events={e.event_id:e for e in doc.events}
    hydrogens=document_hydrogens(doc,result.events)
    jobs=[]
    if not doc.events:
        v=views[next(s.primary_view_id for s in doc.states if s.state_id==doc.initial_state_id)]
        jobs.append((v,None,[],v.view_id,'初始结构 · 形式检查',{},False,True,'structure'))
    for e in doc.events:
        if e.reactant_view_id not in views:continue
        report=reports.get(e.event_id);mode=report.get('execution_mode','electron_actions') if report else 'not_executed'
        target=views.get(e.product_view_id) if report else None
        label={'manual_bridge':'人工状态衔接 · 电子动作未验证','external_event':'外部领域事件 · 非化学执行',
            'not_executed':'未执行／未通过','electron_actions':'真实事件 · 账本通过',
            'declared_electronic_transition':'证据声明的电子态转换 · 非自旋动力学计算',
            'net_transformation':'受限净转化 · 物料闭合，微观电子路径未展开'}[mode]
        if report and report.get('waived_checks'):label='人工规则例外 · 账本通过'
        if result.status=='passed_symbolic' and mode=='electron_actions':label='真实事件 · 局部账本通过，R/n 总账未知'
        if report and report.get('depends_on_exceptions') and not report.get('waived_checks') and mode=='electron_actions':
            label+=' · 依赖前序例外'
        if doc.metadata.get('example_scope')=='inferred_quench_branch':label='推定淬灭 · 候选产物，非原文转写'
        phase=next((p for p in getattr(doc,'phases',[]) if p.phase_id==getattr(e,'phase_id',None)),None)
        if phase and phase.kind=='quench' and phase.conditions.get('label_zh')=='推定淬灭':
            label='推定淬灭 · 候选产物，非原文转写'
        actions=(report.get('compiled_actions',e.actions) or []) if mode=='electron_actions' and report else []
        conditions={**(report.get('effective_conditions',{}) if report else {}),**e.conditions,'说明':e.explanation or '条件或依据未声明'}
        if report and report.get('occupation_assumption') and getattr(e,'general_event',None)=='MN_HYPOCHLORITE_OXO':
            conditions['供对声明']='显式二电子通道；轨道占据未知，不绘制虚构孤对电子点'
            conditions['_declared_pair_donors']=[e.general_roles['m']]
        if report and report.get('declared_pair_donor_ids'):
            conditions['_declared_pair_donors']=report['declared_pair_donor_ids']
            conditions['供对声明']='有来源约束的二电子分子通道；金属未知局域占据不补画孤对电子点'
        transfers=getattr(e,'electron_transfers',[])
        if transfers and report:
            conditions['电子归属转移']='; '.join(f'{t.donor_atom_id} → {t.acceptor_atom_id}: {t.electron_count}e' for t in transfers)+'（局域轨道未知，未虚构轨道弯箭头）'
        jobs.append((views[e.reactant_view_id],target,actions,e.event_id,label+' · '+(getattr(e,'general_event',None) or e.rule_template or e.event_type),
            conditions,False,False,mode))
    for t in doc.resonance_transforms[:config.max_resonance_panels]:
        if t.source_view_id in views and t.target_view_id in views:
            jobs.append((views[t.source_view_id],views[t.target_view_id],t.actions,t.transform_id,'共振表示 · 不增加真实反应步',{},True,False,'resonance'))
    for source,target,actions,key,title,conditions,resonance,state_only,mode in jobs:
        try:
            declared_pair_donors=conditions.pop('_declared_pair_donors',[])
            filename=re.sub(r'[^A-Za-z0-9_-]','_',key)+'_'+hashlib.sha256(key.encode()).hexdigest()[:8]
            panel=draw_panel(source,target,actions,key,title,conditions,directory/filename,config,resonance=resonance,state_only=state_only,
                spin_transition_request=mode=='declared_electronic_transition',
                spin_transition_executed=mode=='declared_electronic_transition',
                declared_pair_donors=declared_pair_donors,event=events.get(key),
                required_hydrogens=hydrogens,inactive=phase_inventory(doc,key),
                released_fragments=reports.get(key,{}).get('released_fragments',[]),
                definitions=doc.metadata.get('drawing_abbreviations',[]))
            # Clearly distinguish non-electron transitions without inventing curves.
            if mode in {'manual_bridge','external_event','net_transformation','declared_electronic_transition'}:
                svg=directory/panel['svg']
                if svg and svg.exists():
                    tree=ET.parse(svg)
                    for node in tree.iter():
                        if node.attrib.get('id')=='reaction-flow':
                            node.set('data-transition-kind',mode)
                            for line in node:
                                if line.tag.endswith('line'):line.set('stroke-dasharray','8 5');line.set('stroke','#a75b17')
                    tree.write(svg,encoding='unicode')
                    svg_to_png(svg,svg.with_suffix('.png'))
            panels.append({**panel,'event_id':key,'execution_mode':mode,'executed':key in reports})
        except Exception as exc:
            panels.append({'event_id':key,'status':'render_failed','reason':str(exc),'execution_mode':mode})
    status='render_failed' if any(p['status']=='render_failed' for p in panels) else 'layout_warning' if any(p['status']!='rendered' for p in panels) else 'rendered'
    from ...drawing.chain import draw_chain
    chain=draw_chain(directory,panels,doc)
    report={'status':status,'panels':panels,'chain':chain,'original_execution_status':result.status,'authority_not_modified':True}
    write_json(directory/'rendering.json',report)
    return report


def export_bundle(result,directory,*,draw=True,_smiles=None,_companion=None):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=False)
    doc=result.document;errors=[];derivatives=[]
    export_mapped_smiles=_smiles or globals()['export_mapped_smiles']
    companion=_companion or sidecar
    if doc.schema_version=='umr/0.8':
        from ..v08.smiles import sidecar as companion
    write_json(directory/'input.json',result.input_document.model_dump(mode='json'))
    write_json(directory/'mechanism.json',doc.model_dump(mode='json'))
    write_json(directory/'result.json',result.to_dict())
    (directory/'mechanism.jsonl').write_text(json.dumps(doc.model_dump(mode='json'),ensure_ascii=False)+'\n',encoding='utf-8')
    (directory/'mechanism.umr').write_text(dump_dsl(doc),encoding='utf-8')
    states={s.state_id:s for s in doc.states}
    for v in doc.views:
        try:
            smi=export_mapped_smiles(v);folder=directory/'smiles';folder.mkdir(exist_ok=True)
            (folder/(v.view_id+'.smi')).write_text(smi+'\n',encoding='utf-8')
            write_json(folder/(v.view_id+'.metadata.json'),companion(states[v.state_id],v,smi));derivatives.append(v.view_id)
        except Exception as exc:errors.append({'view_id':v.view_id,'reason':str(exc)})
    for e in doc.events:
        write_json(directory/'mechsmiles'/(e.event_id+'.metadata.json'),{'event':e.model_dump(mode='json'),
            'schema_version':doc.schema_version,'authority':'../mechanism.json','lossless_without_sidecar':False})
        report=next((r for r in result.events if r['event_id']==e.event_id),{})
        unsupported=report.get('execution_mode') in {'manual_bridge','external_event','net_transformation','declared_electronic_transition'} or any(a.op=='FISHHOOK' for a in e.actions or [])
        if unsupported or e.actions is None:
            write_json(directory/'mechsmiles'/(e.event_id+'.unsupported.json'),{'status':'unsupported','reason':'Requires UMR typed arrows/coordination/exception metadata'})
        else:
            arrows=[]
            for a in e.actions:
                arrows.append(f'({a.a},{a.b})' if a.op=='PAIR' else f'(({a.a},{a.b}),{a.target})' if a.op=='SPLIT' else f'(({a.a},{a.b}),{a.c})')
            (directory/'mechsmiles'/(e.event_id+'.arrows.txt')).write_text(';'.join(arrows)+'\n',encoding='utf-8')
    rendering=render_result(result,directory/'drawing') if draw else None
    manifest={'schema_version':'umr-export/0.7','status':result.status,'input_digest':document_digest(result.input_document),
        'structures_exported':len(derivatives),'derivative_errors':errors,'rendering':rendering,
        'strict_formal_pass':result.status=='passed','human_approved':False,
        'full_electron_path_verified':result.to_dict().get('full_electron_path_verified'),
        'dataset_eligibility':'derivative bundle only; human dataset export requires live approvals and review',
        'reaction_steps':len(result.events),'representation_steps':len(doc.resonance_transforms)}
    write_json(directory/'manifest.json',manifest)
    images=[]
    if rendering:
        from ...drawing.materials import pool_details
        if rendering.get('chain'):
            images.append('<figure><img src="drawing/mechanism_chain.svg"><figcaption>连续机理</figcaption></figure>')
        # Job order is chemical/path order. Lexicographic filenames put om10
        # before om2 in long catalyst paths and must not determine chronology.
        for panel in rendering['panels']:
            if not panel.get('svg'):continue
            p=directory/'drawing'/panel['svg']
            if p.is_file():
                images.append('<figure><img src="'+escape(p.relative_to(directory).as_posix())+'"><figcaption>'+escape(panel['event_id'])+'</figcaption>'+
                    pool_details(p.parent,panel,'drawing/')+'</figure>')
    (directory/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>UMR 0.7</title><style>body{font:16px system-ui;margin:2em;background:#f4f7fa;color:#173148}img{width:100%;height:auto}figure{background:white;padding:1em;border-radius:12px}pre{white-space:pre-wrap}</style><h1>'+escape(doc.record_id)+'</h1><p>'+escape(result.status)+' · 最小形式样例／人工例外草稿，非原案例复测、非已批准机理。</p>'+''.join(images)+'<p><a href="input.json">输入</a> · <a href="result.json">完整执行及严格结果</a> · <a href="mechanism.umr">DSL</a></p>',encoding='utf-8')
    return manifest
