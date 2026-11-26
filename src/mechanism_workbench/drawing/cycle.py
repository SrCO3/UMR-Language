"""Draw one declared turnover around an orbit, with its actual mapped intermediates."""
import hashlib
import math
import re

from .electron_layout import place_electrons
from .language import event_name,METATHESIS_NAMES
from .materials import components, display_plan, phase_inventory
from .structure import structure,place_layout
from .svg import element,label,text_box,write_svg,svg_to_png
from ..umr.exporting import write_json
from ..umr.v09.cycles import mapped_fragment
from .figure_keys import append_keys


def cycle_materials(cycle,source,target,plan,touched,executed):
    exchange={"input":[],"output":[]}
    if not source or not target or not executed:
        return exchange
    for side,view,other,opposite in (("input",source,target,"output"),("output",target,source,"input")):
        core=set(cycle.persistent_core_ids)
        metals={m.atom_id for m in view.metal_centers}
        carrier=core & metals or core
        if cycle.kind=="radical_propagation" and not carrier:
            radicals={a.atom_id for a in view.atoms if a.radical_electrons}
            carrier=radicals or set(cycle.carrier_mapping if side=="input" else cycle.carrier_mapping.values())
        other_groups={tuple(c["atom_ids"]) for c in plan.get(opposite,[])}
        for c in plan.get(side,[]):
            ids=set(c["atom_ids"])
            if not ids & touched or ids & carrier or c["reason"]=="not_yet_added":
                continue
            if tuple(c["atom_ids"]) in other_groups and mapped_fragment(view,ids)==mapped_fragment(other,ids):
                def metal_states(v):
                    return {m.atom_id:(m.oxidation_state,m.d_electrons,m.spin_state) for m in v.metal_centers if m.atom_id in ids}
                if metal_states(view)==metal_states(other):
                    continue
            exchange[side].append(c)
    return exchange


def flow_arrows(parent,event_id,exchange,angle,point,step):
    """Attach each material stream to its turnover event, with an explicit tip."""
    radial=(math.cos(angle),math.sin(angle))
    tangent=(-radial[1],radial[0])
    px,py=point(angle)
    both=all(exchange.values())
    rows=[]
    for side,materials in exchange.items():
        if not materials:
            continue
        incoming=side=="input"
        offset=(-15 if incoming else 15) if both else 0
        def position(r,t):
            return (px+radial[0]*r+tangent[0]*t,py+radial[1]*r+tangent[1]*t)
        inner=position(5,offset)
        outer=position(82,offset*2)
        start,end=(outer,inner) if incoming else (inner,outer)
        vx,vy=end[0]-start[0],end[1]-start[1]
        length=math.hypot(vx,vy)
        ux,uy=vx/length,vy/length
        color="#006b65" if incoming else "#660874"
        all_ids=[i for c in materials for i in c["atom_ids"]]
        g=element("g",data_event_id=event_id,data_flow=side,data_atom_ids=" ".join(map(str,all_ids)))
        parent.append(g)
        g.append(element("line",x1=start[0],y1=start[1],x2=end[0],y2=end[1],stroke=color,stroke_width=2.3))
        wing1=(end[0]-10*ux-4.3*uy,end[1]-10*uy+4.3*ux)
        wing2=(end[0]-10*ux+4.3*uy,end[1]-10*uy-4.3*ux)
        g.append(element("polygon",points=" ".join(f"{x:.2f},{y:.2f}" for x,y in (end,wing1,wing2)),fill=color))
        words=[]
        mark=[]
        for c in materials:
            q=c["charge"]
            charge="" if q==0 else "?" if q is None else (str(abs(q)) if abs(q)>1 else "")+("+" if q>0 else "−")
            if mark:
                words.append((" + ",15,0))
            mark.append(c["formula"]+charge)
            for token in re.findall(r"\d+|[^\d]+",c["formula"]):
                words.append((token,11,4) if token.isdigit() else (token,15,0))
            if charge:
                words.append((charge,11,-5))
        x,y=position(96,offset*2)
        caption=" + ".join(mark)
        # Explicit text positions keep formula digits readable in Qt's
        # offscreen renderer, including systems without Unicode subscripts.
        widths=[text_box(word,0,0,size)[2]-2 for word,size,_ in words]
        if radial[0]<-.45:x-=sum(widths)
        elif radial[0]<=.45:x-=sum(widths)/2
        boxes=[]
        for (word,size,dy),w in zip(words,widths):
            label(g,word,x,y+4+dy,size=size,fill=color)
            boxes.append(text_box(word,x,y+4+dy,size))
            x+=w
        text_bounds=[min(b[0] for b in boxes),min(b[1] for b in boxes),max(b[2] for b in boxes),max(b[3] for b in boxes)]
        rows.append(dict(event_id=event_id,step=step,side=side,atom_groups=[c["atom_ids"] for c in materials],
                         start=start,end=end,orbit_point=[px,py],label=caption,label_box=text_bounds,
                         arrow_box=[min(start[0],end[0])-5,min(start[1],end[1])-5,max(start[0],end[0])+5,max(start[1],end[1])+5]))
    return rows


EXCHANGE_EVENTS={
    "OM_PRODUCT_RELEASE","OM_DONOR_DISSOCIATION","OM_OLEFIN_ASSOCIATION",
    "COORDINATION","DISSOCIATION","LIGAND_ASSOCIATION","LIGAND_DISSOCIATION",
    "LIGAND_COORDINATION","L_COORDINATION","L_DISSOCIATION",
    "LA_L_ASSOCIATION","LA_L_DISSOCIATION","BORON_L_ASSOCIATION","BORON_L_DISSOCIATION",
    "AL_CYANIDE_ASSOCIATION","ANION_METAL_ASSOCIATION",
}


def turnover_name(event,materials):
    code=event.general_event or event.metal_template or event.element_template
    if (code or event.template_id or event.event_type) in EXCHANGE_EVENTS and any(materials.values()):
        return ""
    if not code:
        return event_name(event)
    return (event.conditions.get("step_title") or METATHESIS_NAMES.get(code) or
            {"OXIDATIVE_ADDITION_2E":"氧化加成","REDUCTIVE_ELIMINATION_2E":"还原消除",
             "TRANSMETALATION_PD_ZN":"转金属化","LIGAND_SUBSTITUTION":"配体取代",
             "CC_BASE_LIGAND_EXCHANGE":"碱辅助的配体交换","CC_ANION_LIGAND_EXCHANGE":"阴离子配体交换"}.get(code) or
            (event_name(event) if code==event.event_type else code.replace("_"," ")))


def orbit_label(parent,step,angle,point,center,radius,occupied,bounds):
    """Reserve an unoccupied horizontal text box near the actual arc arrowhead."""
    title=step["label"]
    if not title:
        return True
    cx,cy=center
    px,py=point(angle)
    ux,uy=math.cos(angle),math.sin(angle)
    for size in (15,13,11):
        lines=[""]
        for ch in title:
            if lines[-1] and text_box(lines[-1]+ch,0,0,size)[2]>160:
                lines.append("")
            lines[-1]+=ch
        font_boxes=[text_box(line,0,0,size) for line in lines]
        w=max(b[2] for b in font_boxes)+4
        line_height=size+4
        h=line_height*len(lines)+4
        half=abs(ux)*w/2+abs(uy)*h/2
        options=[]
        for clearance in (14,25,42,64,90,120,150):
            for sign in (-1,1):
                r=sign*(half+clearance)
                for along in (0,-14,14,-32,32,-58,58,-85,85):
                    x,y=px+r*ux-along*uy,py+r*uy+along*ux
                    box=(x-w/2,y-h/2,x+w/2,y+h/2)
                    a,b,c,d=box
                    if a<12 or b<90 or c>bounds[0]-12 or d>bounds[1]-12:
                        continue
                    if any(a<z+6 and c>u-6 and b<t+6 and d>v-6 for u,v,z,t in occupied):
                        continue
                    near=math.hypot(max(a-cx,0,cx-c),max(b-cy,0,cy-d))
                    far=max(math.hypot(i-cx,j-cy) for i in (a,c) for j in (b,d))
                    if far>radius-10 and near<radius+10:
                        continue
                    options.append((math.hypot(r,along)+(2 if sign>0 else 0),x,y,box))
        if options:
            _,x,y,box=min(options)
            group=element("g",data_event_label=step["event_id"])
            parent.append(group)
            for k,line in enumerate(lines):
                label(group,line,x,y-(len(lines)-1)*line_height/2+k*line_height+size*.32,
                      size=size,text_anchor="middle",fill="#006b65")
            occupied.append(box)
            step["label_box"]=box
            step["arrow_tip"]=[px,py]
            return True
    return False


def molecule_in_box(parent,view,ids,touched,config,definitions,box,caption):
    x,y,w,h=box
    # Interrupt the turnover curve behind the structure and its caption.
    parent.append(element("rect",x=x,y=y,width=w,height=h,fill="white"))
    label(parent,caption,x+w/2,y+21,size=14,text_anchor="middle",fill="#660874")
    if not view or not ids:
        label(parent,"结构待核对",x+w/2,y+h/2,size=14,text_anchor="middle")
        return []
    drawing=structure(view,config,visible_ids=ids,touched=touched,definitions=definitions)
    body=element("g",data_view_id=view.view_id)
    numbers=element("g",style="" if config.show_atom_ids else "display:none")
    coords,obs,_=place_layout(drawing,body,(0,0),view,"cycle",caption,numbers,config.show_atom_ids)
    body.append(numbers)
    extra,_=place_electrons(drawing["display_view"],[],coords,body,obs,config,side="cycle")
    boxes=[o["box"] for o in obs+extra]+[[px-8,py-8,px+8,py+8] for px,py in coords.values()]
    if not boxes:
        return []
    xmin=min(b[0] for b in boxes);ymin=min(b[1] for b in boxes)
    width=max(b[2] for b in boxes)-xmin;height=max(b[3] for b in boxes)-ymin
    scale=min((w-28)/max(width,1),(h-43)/max(height,1),1.15)
    tx=x+(w-width*scale)/2;ty=y+30+(h-36-height*scale)/2
    body.set("transform",f"translate({tx} {ty}) scale({scale}) translate({-xmin} {-ymin})")
    parent.append(body)
    return drawing["labels"]


def draw_cycle(document,report,cycle,path,config,executed=None,*,spacing=1.0):
    from ..umr.v09.local_rules import use_rules
    with use_rules(document):
        return _draw_cycle(document,report,cycle,path,config,executed,spacing=spacing)

def _draw_cycle(document,report,cycle,path,config,executed,*,spacing):
    views={v.view_id:v for v in document.views}
    events={e.event_id:e for e in document.events}
    closed=report.get("status")=="closed_one_turnover" and bool(report.get("compressible"))
    stages=[(eid,events[eid].reactant_view_id if eid in events else None) for eid in cycle.event_ids]
    if not closed:
        stages.append((None,cycle.final_view_id))
    n=len(stages)
    if not n:
        raise ValueError("Cycle contains no events or states to depict")
    tw,th=266,168
    slots=max(3,n+(0 if closed else 1))
    radius=max(310,365/(2*math.sin(math.pi/slots)))*spacing
    cx=radius+tw/2+140;cy=radius+th/2+170
    width=math.ceil(2*cx);bottom=math.ceil(cy+radius+th/2+92)
    height=bottom+(35 if cycle.termination_event_ids else 0)
    # Retain vector detail for long cycles while keeping the PNG pixel budget.
    screen=min(1,math.sqrt(12_000_000/(width*height)))
    root=element("svg",width=round(width*screen),height=round(height*screen),viewBox=f"0 0 {width} {height}")
    root.append(element("rect",width="100%",height="100%",fill="white"))
    orbit=element("g",id="cycle-orbit",data_closure=report["status"])
    root.append(orbit)
    point=lambda theta:(cx+radius*math.cos(theta),cy+radius*math.sin(theta))
    angles=[-math.pi/2+2*math.pi*k/slots for k in range(n)]
    if closed:
        orbit.append(element("circle",cx=cx,cy=cy,r=radius,fill="none",stroke="#018b83",stroke_width=3))
    else:
        # A missing turnover has a visible gap, even if the same metal recurs.
        ox,oy=point(angles[0]);ex,ey=point(angles[-1])
        orbit.append(element("path",d=f"M {ox} {oy} A {radius} {radius} 0 {int(angles[-1]-angles[0]>math.pi)} 1 {ex} {ey}",
                             fill="none",stroke="#a28352",stroke_width=2,stroke_dasharray="9 6"))
    nodes=element("g",id="cycle-intermediates");root.append(nodes)
    streams=element("g",id="cycle-material-flows");root.append(streams)
    notes=element("g",id="cycle-labels");root.append(notes)
    title="自由基传播循环" if cycle.kind=="radical_propagation" else "催化循环"
    label(notes,title,24,32,size=25,fill="#660874")
    label(notes,"一次周转已核验" if closed else "循环尚未核验闭合 · "+str(report["status"]),24,57,size=16)
    if cycle.activation_event_ids:
        label(notes,"环外活化／引发  "+" → ".join(cycle.activation_event_ids)+" → "+cycle.initial_view_id,
              24,84,size=14,fill="#536761")
    defs=document.metadata.get("drawing_abbreviations",[])
    records=[];steps=[];flows=[];all_folds=[]
    occupied=[(cx-130,cy-86,cx+130,cy+135)]
    for k,(eid,vid) in enumerate(stages):
        view=views.get(vid)
        e=events.get(eid)
        target=views.get(e.product_view_id) if e else None
        touched=set()
        plan={}
        if view and e:
            plan,touched,_=display_plan(view,target,e.actions or [],e,[],phase_inventory(document,eid),config,False)
        ids=set()
        if view:
            cores=set(cycle.persistent_core_ids)
            metals={m.atom_id for m in getattr(view,"metal_centers",[])} & cores
            carrier=metals or cores or (set(cycle.carrier_mapping.values()) if not eid else set(cycle.carrier_mapping))
            groups=components(view)
            ids=set().union(*(g for g in groups if g & carrier))
            if not ids:
                ids=set().union(*(g for g in groups if g & touched))
        x,y=point(angles[k])
        occupied.append((x-tw/2,y-th/2,x+tw/2,y+th/2))
        card=element("g",data_view_id=vid or "",data_event_id=eid or "")
        nodes.append(card)
        known=eid and (executed is None or eid in executed)
        folds=molecule_in_box(card,view,ids,touched,config,defs,(x-tw/2,y-th/2,tw,th),
                              ("中间体 " if known else "声明视图 ")+(vid or "未知"))
        records.append(dict(view_id=vid,event_id=eid,atom_ids=sorted(ids),abbreviations=folds))
        all_folds.extend(folds)
        if not eid:
            continue
        mid=angles[k]+math.pi/slots;px,py=point(mid)
        ux,uy=-math.sin(mid),math.cos(mid)
        pending=executed is not None and eid not in executed
        color="#a48f72" if pending or not closed else "#018b83"
        orbit.append(element("path",d=f"M {px-ux*9-uy*5} {py-uy*9+ux*5} L {px} {py} L {px-ux*9+uy*5} {py-uy*9-ux*5}",
                             fill="none",stroke=color,stroke_width=3,data_event_id=eid))
        verified=eid in executed if executed is not None else closed
        exchange=cycle_materials(cycle,view,target,plan,touched,verified)
        flows.extend(flow_arrows(streams,eid,exchange,mid,point,k+1))
        steps.append(dict(event_id=eid,input=vid,output=e.product_view_id if e else None,materials=exchange,
                          label=turnover_name(e,exchange) if e else "缺少事件",
                          executed=None if executed is None else not pending))
    final=views.get(cycle.final_view_id)
    products=set(i for group in cycle.product_atom_groups for i in group)
    if final and products and products<={a.atom_id for a in final.atoms}:
        all_folds.extend(molecule_in_box(notes,final,products,set(),config,defs,(cx-125,cy-81,250,153),
                        "周转产物" if closed else "声明的产物"))
    exchange=[f"{i}→{j}" for i,j in cycle.carrier_mapping.items() if i!=j]
    if closed:
        label(notes,f"{cycle.final_view_id} → {cycle.initial_view_id}  载体再生",cx,cy+101,size=15,text_anchor="middle")
        if exchange:label(notes,"可交换原子  "+", ".join(exchange),cx,cy+125,size=13,text_anchor="middle")
    elif report.get("reason"):
        label(notes,"缺口见循环记录",cx,cy+108,size=14,text_anchor="middle")
    if cycle.termination_event_ids:
        label(notes,"环外终止／失活  "+" → ".join(cycle.termination_event_ids),24,bottom+8,size=14)
    for f in flows:
        occupied.extend((f["arrow_box"],f["label_box"]))
    for k,s in enumerate(steps):
        if not orbit_label(notes,s,angles[k]+math.pi/slots,point,(cx,cy),radius,occupied,(width,bottom)):
            if spacing>=2.5:
                raise ValueError("Turnover event label exceeds available diagram space: "+s["event_id"])
            return draw_cycle(document,report,cycle,path,config,executed,spacing=spacing*1.35)
    keys=append_keys(root,all_folds)
    from .local_rules import append_rules
    append_rules(root,[vid for _,vid in stages],[events[eid] for eid,_ in stages if eid in events])
    write_svg(root,path)
    svg_to_png(path,path.with_suffix(".png"))
    result=dict(cycle_id=cycle.cycle_id,status=report["status"],closed=closed,svg=path.name,png=path.with_suffix(".png").name,
                initial_view=cycle.initial_view_id,final_view=cycle.final_view_id,intermediates=records,events=steps,
                activation=cycle.activation_event_ids,termination=cycle.termination_event_ids,product_atom_ids=sorted(products),
                orbit_geometry=dict(center=[cx,cy],radius=radius),
                carrier_mapping=cycle.carrier_mapping,material_flows=flows,abbreviation_key=keys,validation=report)
    write_json(path.with_suffix(".cycle.json"),result)
    return result


def draw_cycles(result,directory,config):
    cycles={c.cycle_id:c for c in getattr(result.document,"cycles",[])}
    rows=[]
    executed={r["event_id"] for r in result.events}
    for report in getattr(result,"cycle_reports",[]):
        cycle=cycles[report["cycle_id"]]
        key=hashlib.sha256(cycle.cycle_id.encode()).hexdigest()[:12]
        rows.append(draw_cycle(result.document,report,cycle,directory/f"cycle-{key}.svg",config,executed))
    return rows
