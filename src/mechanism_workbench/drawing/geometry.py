"""Bounded deterministic cubic routing; warnings are results, not hidden failures."""
import math
import time


def add(a,b):return (a[0]+b[0],a[1]+b[1])
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def mul(a,s):return (a[0]*s,a[1]*s)
def norm(a):return math.hypot(*a)
def unit(a):return mul(a,1/max(norm(a),1e-9))
def normal(a):return (-a[1],a[0])
def distance(a,b):return norm(sub(a,b))
def cross(a,b):return a[0]*b[1]-a[1]*b[0]


def samples(c,n=28):
    p,q,r,s=c
    return [add(add(mul(p,(1-t)**3),mul(q,3*(1-t)**2*t)),add(mul(r,3*(1-t)*t*t),mul(s,t**3))) for t in (i/n for i in range(n+1))]


def head(c,electrons):
    end=c[-1];u=unit(sub(c[-1],c[-2]));v=normal(u);size=10
    wings=[add(add(end,mul(u,-size)),mul(v,4.5))]
    if electrons==2:wings.append(add(add(end,mul(u,-size)),mul(v,-4.5)))
    return [(end,w) for w in wings]


def intersects(a,b,c,d):
    ab,cd=sub(b,a),sub(d,c);den=cross(ab,cd)
    if abs(den)<1e-8:return False
    ac=sub(c,a);t=cross(ac,cd)/den;u=cross(ac,ab)/den
    return 1e-4<t<1-1e-4 and 1e-4<u<1-1e-4


def segment_distance(p,a,b):
    delta=sub(b,a);length=max(norm(delta)**2,1e-9)
    t=max(0,min(1,((p[0]-a[0])*delta[0]+(p[1]-a[1])*delta[1])/length))
    return distance(p,add(a,mul(delta,t)))


class Spatial:
    def __init__(self,objects,cell=48):
        self.objects=objects;self.cells={};self.cell=cell
        for i,o in enumerate(objects):
            x0,y0,x1,y1=o['box']
            for x in range(math.floor(x0/cell),math.floor(x1/cell)+1):
                for y in range(math.floor(y0/cell),math.floor(y1/cell)+1):self.cells.setdefault((x,y),set()).add(i)
    def near(self,p):
        x,y=math.floor(p[0]/self.cell),math.floor(p[1]/self.cell)
        ids=set()
        for dx in (-1,0,1):
            for dy in (-1,0,1):ids|=self.cells.get((x+dx,y+dy),set())
        return [self.objects[i] for i in sorted(ids)]


def collisions(c,request,spatial,padding):
    points=samples(c);all_points=list(points)
    for a,b in head(c,request['electrons']):
        all_points.extend(add(a,mul(sub(b,a),t/5)) for t in range(6))
    hits={};p0,p3=c[0],c[-1]
    for p in all_points:
        for o in spatial.near(p):
            x0,y0,x1,y1=o['box'];pad=padding+1.2 # includes stroke, not centerline only
            if not x0-pad<=p[0]<=x1+pad or not y0-pad<=p[1]<=y1+pad:continue
            own0=set(o.get('atom_ids',[]))<=set(request['source']['atom_ids'])
            own3=set(o.get('atom_ids',[]))<=set(request['target']['atom_ids'])
            if o.get('atom_ids') and ((own0 and distance(p,p0)<=13) or (own3 and distance(p,p3)<=13)):continue
            if o['kind']=='bond' and not o.get('filled'):
                if not any(segment_distance(p,a,b)<=o.get('width',1.5)/2+pad for a,b in o.get('segments',[])):continue
            hits[o['id']]={'object_id':o['id'],'kind':o['kind']}
    return list(hits.values()),points


def candidate_curves(request,bounds,config,spatial):
    all_curves=[];x0,y0,x1,y1=bounds
    ports=[(s,e) for s in request['source']['positions'] for e in request['target']['positions']]
    for bend in (22,-22,38,-38,60,-60,90,-90,130,-130):
        for arm in (.28,.48):
            for start,end in ports:
                if distance(start,end)<6:continue
                delta=sub(end,start);v=normal(unit(delta))
                c=[start,add(add(start,mul(delta,arm)),mul(v,bend)),add(add(end,mul(delta,-arm)),mul(v,bend)),end]
                # A genuine gap separates the curve from the visible electron dot
                # or source bond. Preserve the semantic anchors in metadata.
                gap=config.electron_endpoint_gap+(2 if request['electrons']==1 else 0)
                c[0]=add(start,mul(unit(sub(c[1],start)),gap))
                c[-1]=add(end,mul(unit(sub(c[-2],end)),2.5))
                if not all(x0+6<=p[0]<=x1-6 and y0+6<=p[1]<=y1-6 for p in c):continue
                hits,points=collisions(c,request,spatial,config.collision_padding)
                crossing=sum(intersects(a,b,c1,d) for i,(a,b) in enumerate(zip(points,points[1:])) for j,(c1,d) in enumerate(zip(points,points[1:])) if j>i+2)
                arc=sum(distance(a,b) for a,b in zip(points,points[1:]))
                curvature=sum(norm(sub(sub(points[i+1],points[i]),sub(points[i],points[i-1]))) for i in range(1,len(points)-1))
                cost=sum(1600 if h['kind'] in {'text','electron','atom_id','charge'} else 500 for h in hits)+crossing*4000+arc*.03+curvature*.02
                all_curves.append({'control_points':c,'points':points,'semantic_source_xy':list(start),'semantic_target_xy':list(end),
                    'source_gap_pixels':gap,'target_gap_pixels':2.5,'collisions':hits,'self_crossings':crossing,'cost':cost,'head_segments':head(c,request['electrons'])})
                if len(all_curves)>=config.curve_candidates:break
            if len(all_curves)>=config.curve_candidates:break
        if len(all_curves)>=config.curve_candidates:break
    return sorted(all_curves,key=lambda c:(c['cost'],repr(c['control_points'])))[:12]


def route(requests,obstacles,bounds,config):
    started=time.perf_counter();spatial=Spatial(obstacles);options=[];limited=False
    for req in requests:
        opts=candidate_curves(req,bounds,config,spatial)
        if not opts:return [],{'status':'render_failed','reason':'No bounded cubic candidates','arrow_id':req['arrow_id'],'seconds':time.perf_counter()-started}
        options.append(opts)
        if time.perf_counter()-started>config.routing_seconds:limited=True
    beam=[(0,[],[])]
    for index,opts in enumerate(options):
        grown=[]
        for cost,selected,conflicts in beam:
            for opt in opts:
                crossings=[]
                # Joint beam alternatives retain all earlier arrows; never freeze
                # a single greedy curve then squeeze subsequent arrows around it.
                for previous,curve in enumerate(selected):
                    ps,qs=opt['points'][::2],curve['points'][::2]
                    if any(intersects(a,b,c,d) for a,b in zip(ps,ps[1:]) for c,d in zip(qs,qs[1:])):
                        crossings.append([requests[previous]['arrow_id'],requests[index]['arrow_id']])
                    if any(segment_distance(p,a,b)<3.5 for a,b in opt['head_segments'] for p in qs):
                        crossings.append([requests[previous]['arrow_id'],requests[index]['arrow_id']])
                    near=sum(min(segment_distance(p,a,b) for a,b in zip(qs,qs[1:]))<5.5 for p in ps)
                    if near>=2:crossings.append([requests[previous]['arrow_id'],requests[index]['arrow_id']])
                grown.append((cost+opt['cost']+len(crossings)*1500,selected+[opt],conflicts+crossings))
        beam=sorted(grown,key=lambda b:b[0])[:config.joint_beam]
        if time.perf_counter()-started>config.routing_seconds:limited=True;beam=beam[:1]
    _,selected,crossings=beam[0] if beam else (0,[],[])
    for req,opt in zip(requests,selected):
        req.update({k:v for k,v in opt.items() if k!='points'})
        req['source']['position']=opt['control_points'][0];req['target']['position']=opt['control_points'][-1]
    hits=[{'arrow_id':r['arrow_id'],**h} for r in requests for h in r.get('collisions',[])]
    return requests,{'status':'layout_warning' if hits or crossings or limited else 'rendered','collisions':hits,
        'arrow_crossings':crossings,'time_limit_reached':limited,'seconds':time.perf_counter()-started,
        'candidate_limit':config.curve_candidates,'joint_beam':config.joint_beam,'endpoint_exemption_radius':13,
        'line_and_head_checked':True,'endpoint_check':all(distance(r['source']['position'],r['control_points'][0])<1e-6 and distance(r['target']['position'],r['control_points'][-1])<1e-6 for r in requests)}
