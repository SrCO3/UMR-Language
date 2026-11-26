from pathlib import Path
from html import escape
from ..v07.exporting import export_bundle as shared
from ..exporting import write_json
from .smiles import export_mapped_smiles,sidecar

def cycle_svg(document,report,cycle,path):
    from ...drawing.cycle import draw_cycle
    from ...drawing.models import RenderConfig
    return draw_cycle(document,report,cycle,Path(path),RenderConfig())

def export_bundle(result,directory,*,draw=True):
    manifest=shared(result,directory,draw=draw,_smiles=export_mapped_smiles,_companion=sidecar)
    directory=Path(directory)
    manifest.update(schema_version='umr-export/0.9',cycles=result.cycle_reports,formal_path_status=result.status,
       full_named_mechanism_claim=False,source_scope=result.input_document.metadata.get('example_scope','source_scope_requires_review'))
    stereo=getattr(result,'stereo_path',{})
    manifest['stereochemistry']=stereo
    write_json(directory/'stereo-path.json',stereo)
    if draw:
        from ...drawing.cycle import draw_cycles
        from ...drawing.models import RenderConfig
        manifest['cycle_diagrams']=draw_cycles(result,directory/'drawing',RenderConfig())
    write_json(directory/'manifest.json',manifest)
    index=directory/'index.html'
    html=index.read_text(encoding='utf-8').replace('UMR 0.7','UMR 0.9').replace('最小形式样例／人工例外草稿，非原案例复测、非已批准机理。','执行输出；来源范围见 manifest.json。形式通过不代表来源已审核或完整人名机理成立。')
    if stereo:
        warning=stereo.get('uncertain_events',[])
        message=('形式执行通过；下列中间态的构型或选择性仍不明确，可以保留继续执行。' if result.execution_complete else '形式执行未完成，原因见事件报告；下列立体不确定性另行保留。') if warning else '已记录各事件的局部立体关系；CIP名称变化不自动等于翻转。'
        rows=''.join('<tr><td>'+escape(str(x['event_id']))+'</td><td>'+str(x.get('atom_id',x.get('atom_ids')))+'</td><td>'+escape(x['effect'])+'</td><td>'+escape(x['reason'])+'</td></tr>' for x in warning)
        block='<section style="padding:16px;background:#fff5d9"><h2>立体化学审查</h2><p>'+message+'</p><table><tr><th>事件</th><th>原子</th><th>关系</th><th>待定原因</th></tr>'+rows+'</table><p>未知构型不代替化学证据；产物反推仅保存为独立推断记录。</p></section>'
        html=html.replace('<h1>',block+'<h1>',1)
    for p in (directory/'drawing').glob('cycle-*.svg'):html+='<figure><img src="drawing/'+p.name+'" alt="循环事件图"></figure>'
    proposal=getattr(result,'inferred_quench',None)
    if proposal and proposal['status'] not in {'explicit_workup_present','not_proposed_incomplete_reaction'}:
        write_json(directory/'inferred-quench.json',proposal)
        html+='<section><h2>推定淬灭</h2><p>最后加入水的条件化候选；不是原文已给出的步骤，未经人工审核。原始路径与推定结果分开保存。</p>'
        from .document import run_document
        from .ast import MechanismDocument
        for k,row in enumerate(proposal['candidates']):
            if not row.get('input_document'):continue
            qr=run_document(MechanismDocument.model_validate(row['input_document']),max_views=4,max_expansions=16,propose_quench=False)
            sub=directory/f'inferred-quench-{k+1}';export_bundle(qr,sub,draw=draw)
            html+=f'<p><a href="inferred-quench-{k+1}/index.html">推定淬灭候选 {k+1}（位点 {row["site"]}）及产物</a></p>'
            if draw:
                for p in (sub/'drawing').glob('*.svg'):html+='<figure><img src="'+p.relative_to(directory).as_posix()+'" alt="推定淬灭及产物"></figure>'
        if not proposal['candidates']:html+='<p>'+escape(proposal['no_change_note'])+'</p>'
        html+='</section>'
    index.write_text(html,encoding='utf-8');return manifest
