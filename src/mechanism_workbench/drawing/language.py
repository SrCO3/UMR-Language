"""Chinese reviewer-facing labels generated from actual actions and provenance."""
EVENT_NAMES={'SN2':'亲核取代','POLAR':'双电子极性反应','PROTON_TRANSFER':'质子转移','CARBONYL_ADDITION':'羰基加成',
 'HETEROLYSIS':'化学键异裂','E2':'协同消除','E1_IONIZATION':'离去基团解离','E1_ELIMINATION':'碳正离子消除',
 'E1CB_DEPROTONATION':'脱质子生成碳负离子','E1CB_ELIMINATION':'碳负离子消除','RADICAL':'自由基反应','PERICYCLIC':'周环反应',
 'COORDINATION':'配体配位','OXIDATIVE_ADDITION':'氧化加成','REDUCTIVE_ELIMINATION':'还原消除','TRANSMETALLATION':'转金属化'}
TEMPLATES={'HAT':'自由基夺氢','HOMOLYSIS':'化学键均裂','RADICAL_ADD_PI':'自由基加成','DA_4_2':'狄尔斯–阿尔德环加成',
 'CYCLOADDITION_2_2':'光照［2+2］环加成','L_COORDINATION':'配体配位','L_DISSOCIATION':'配体解离'}
REACTIVE_TEMPLATES={'CARBENE_ALPHA_ELIMINATION':'α-消除生成卡宾 · Hine 机理类比',
 'CARBENE_CYCLOPROPANATION':'单重态卡宾 · 协同环丙烷化','NITRENE_AZIRIDINATION':'单重态氮宾 · 协同氮杂环丙烷化',
 'CARBENE_CH_INSERTION':'单重态卡宾 · 指定 C—H 插入','NITRENE_CH_INSERTION':'单重态氮宾 · 指定 C—H 插入',
 'TRIPLET_CARBENE_FIRST_ADDITION':'三重态卡宾 · 首次加成，保留双自由基',
 'TRIPLET_NITRENE_FIRST_ADDITION':'三重态氮宾 · 首次加成，保留双自由基',
 'CARBYNE_DOUBLET_HAT':'双重态卡拜夺氢 · 仅形式账本模型',
 'ARYNE_DEPROTONATION':'苯炔路径 · 指定邻位氢脱质子', 'ARYNE_ELIMINATION':'消除生成苯炔 · Roberts 路径类比',
 'ARYNE_NUCLEOPHILIC_ADDITION':'苯炔 · 氨基负离子亲核加成','ARYL_CARBA_PROTONATION':'芳基碳负离子 · 显式氢质子转移'}
SOURCES={'rmechdb_data':'RMechDB 自由基反应机理数据集','flowER':'FlowER 电子流机理数据集',
 'ReactionMechanismReasoning-main':'反应机理推理数据集','pmechrp_datasets':'pMechRP 机理预测数据集',
 'oMeBench-main':'oMeBench 有机机理评测集','ChRIMP-main':'ChRIMP 反应机理数据'}
ELEMENTS={'C':'碳','O':'氧','N':'氮','H':'氢','Br':'溴','Cl':'氯','F':'氟','S':'硫','P':'磷','Pd':'钯'}


METATHESIS_NAMES={
    "OM_DONOR_DISSOCIATION":"给体配体解离","OM_OLEFIN_ASSOCIATION":"烯烃配位",
    "OM_METALLACYCLE_FORMATION":"金属杂环丁烷形成","OM_METALLACYCLE_CLEAVAGE":"金属杂环丁烷裂解",
    "OM_PRODUCT_RELEASE":"烯烃产物解离",
}


def event_name(event):
    if hasattr(event,'model_dump'):event=event.model_dump(mode='json')
    return event.get('conditions',{}).get('step_title') or METATHESIS_NAMES.get(event.get('general_event')) or REACTIVE_TEMPLATES.get(event.get('reactive_template')) or TEMPLATES.get(event.get('template_id')) or EVENT_NAMES.get(event.get('event_type'),'电子转移')


def source_name(entry):
    source=entry.get('source',{})
    if source.get('fixture') or entry.get('conversion_nature','').startswith('hand-authored'):
        return '本项目教学示例 · 人工编写与形式回放'
    return SOURCES.get(source.get('dataset_id'),'来源信息待补充')


def conditions_text(conditions):
    if not conditions:return '条件尚未注明'
    names={'activation':'激活方式','solvent':'溶剂','temperature_kelvin':'温度','catalyst':'催化剂',
        'photochemical':'光照','thermal':'热反应','unknown':'未知','step_title':'步骤'}
    def present(k,v):
        if k=='step_title' or isinstance(v,(dict,list)):return False
        if ('spin' in k.lower() or '自旋' in k or k=='multiplicity') and str(v).lower() in {'unknown','unspecified','none','未知','待定'}:return False
        return True
    return ' · '.join(f"{names.get(k,k)}：{names.get(str(v),str(v))}" for k,v in conditions.items() if present(k,v)) or '条件尚未注明'


def describe_actions(view,actions):
    atoms={a.atom_id:a for a in view.atoms} if hasattr(view,'atoms') else {a['atom_id']:a for a in view.get('atoms',[])}
    def name(i):
        a=atoms.get(i);element=a.element if hasattr(a,'element') else a.get('element','?') if a else '?'
        return ELEMENTS.get(element,element)+'原子'
    result=[]
    for a in actions:
        a=a.model_dump(mode='json') if hasattr(a,'model_dump') else a;op=a.get('op')
        if op=='PAIR':text=f"{name(a['a'])}提供一对孤对电子，向{name(a['b'])}成键。"
        elif op=='SPLIT':text=f"{name(a['a'])}与{name(a['b'])}之间的键断开一对电子，电子归于{name(a['target'])}。"
        elif op=='SHIFT':text=f"原有键上的一对电子经{name(a['b'])}转向新的成键位置，旧键键级降低、新键键级升高。"
        elif op=='FISHHOOK':
            s,d=a['source'],a['destination'];origin=f"{name(s['atom_id'])}的未配对电子" if s['kind']=='radical' else '原有键中的一个电子'
            target=f"{name(d['atom_id'])}" if d['kind']=='atom' else '新的成键位置'
            text=f'{origin}转移至{target}，使用单电子鱼钩箭头。'
        else:text='此电子动作尚未支持绘图，原始记录已保留。'
        result.append(text)
    return result
