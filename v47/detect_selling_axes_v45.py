#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import csv,json,math,hashlib
import cv2,numpy as np
B=Path(__file__).resolve().parent

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
# Hash-locked inputs.
G=B/'selling_fig1_gray_300_v37.png'; P=B/'Tafel-III-fig1.jpg'; PG=B/'selling_fig1_bsb_photo_gate_v38.json'
im=cv2.imread(str(G),0); full=cv2.imread(str(P),0)
assert im is not None and full is not None
assert sha(P)=='36fbabf24d710e4bfd570bf05b2dba8350862ffcd02d44efd4f755204653a2b0'
pg=json.loads(PG.read_text(encoding='utf-8'))
H=np.array(pg['registration_to_v37']['H_photo_to_v37_300'],float); Hinv=np.linalg.inv(H)
# v38 canonical photo crop.
photo=full[120:3000,1460:3550]
bg=cv2.GaussianBlur(photo,(0,0),sigmaX=18,sigmaY=18)
dark=np.clip(bg.astype(np.int16)-photo.astype(np.int16),0,255).astype(np.uint8)
photo_ink=cv2.dilate((dark>12).astype(np.uint8),np.ones((5,5),np.uint8),1)
# Source-style dashed-line detector: v37 LSD grouping, expanded only enough to expose
# the ambiguity that was hidden by v37's conservative threshold.
lsd=cv2.createLineSegmentDetector(cv2.LSD_REFINE_STD); L=lsd.detect(im)[0]
short=[]
for z in (L[:,0,:] if L is not None else []):
    x1,y1,x2,y2=map(float,z); le=math.hypot(x2-x1,y2-y1)
    if le<7: continue
    mx=(x1+x2)/2; my=(y1+y2)/2
    if mx<110 or my>2520 or (my<220 and mx<1300): continue
    ang=math.atan2(y2-y1,x2-x1); ang=ang+math.pi if ang<0 else ang
    dx,dy=math.cos(ang),math.sin(ang); nx,ny=-dy,dx
    rho=nx*mx+ny*my; t1=dx*x1+dy*y1; t2=dx*x2+dy*y2
    if t1>t2:t1,t2=t2,t1
    short.append((x1,y1,x2,y2,le,ang,rho,t1,t2))
bins=defaultdict(list)
for s in short: bins[(round(math.degrees(s[5])/2.5),round(s[6]/7))].append(s)
cands=[]
for key,arr in bins.items():
    if len(arr)<4: continue
    a=float(np.median([s[5] for s in arr])); dx,dy=math.cos(a),math.sin(a); nx,ny=-dy,dx
    rho=float(np.median([s[6] for s in arr])); ints=sorted((s[7],s[8]) for s in arr)
    merged=[]
    for lo,hi in ints:
        if not merged or lo>merged[-1][1]+3: merged.append([lo,hi])
        else: merged[-1][1]=max(merged[-1][1],hi)
    if len(merged)<4: continue
    tmin=min(x[0] for x in merged); tmax=max(x[1] for x in merged); span=tmax-tmin
    cov=sum(hi-lo for lo,hi in merged)/max(span,1)
    # Declared v45 exploratory window. Relaxation is diagnostic only and never a historical promotion.
    if span<250 or not(0.04<=cov<=0.30): continue
    x1=rho*nx+tmin*dx; y1=rho*ny+tmin*dy; x2=rho*nx+tmax*dx; y2=rho*ny+tmax*dy
    cands.append(dict(x1=x1,y1=y1,x2=x2,y2=y2,span=span,coverage=cov,pieces=len(merged),angle_deg=math.degrees(a),rho=rho))
uniq=[]
for c in sorted(cands,key=lambda q:-(q['span']*q['pieces'])):
    if all(not(min(abs(c['angle_deg']-d['angle_deg']),180-abs(c['angle_deg']-d['angle_deg']))<2 and abs(c['rho']-d['rho'])<12) for d in uniq):
        uniq.append(c)
    if len(uniq)>=18: break
# Dash-aware dual-source support: use only actual dark samples on v37 raster, project them into BSB photo.
for i,c in enumerate(uniq):
    p=np.array([c['x1'],c['y1']],float); q=np.array([c['x2'],c['y2']],float); Lq=float(np.linalg.norm(q-p)); n=max(2,int(Lq)+1)
    pts=np.linspace(p,q,n); srcdark=[]; photo_vals=[]
    for x,y in pts:
        xx=int(round(x)); yy=int(round(y)); v=False
        if 0<=xx<im.shape[1] and 0<=yy<im.shape[0]: v=bool(im[yy,xx]<200)
        srcdark.append(v)
        if not v: continue
        ph=Hinv@np.array([x,y,1.0]); xp,yp=ph[0]/ph[2],ph[1]/ph[2]; xxp=int(round(xp)); yyp=int(round(yp))
        if 0<=xxp<photo_ink.shape[1] and 0<=yyp<photo_ink.shape[0]: photo_vals.append(bool(photo_ink[yyp,xxp]))
    arr=np.array(srcdark,bool)
    c['source_dark_fraction']=float(arr.mean())
    c['source_dark_samples']=int(arr.sum())
    c['bsb_support_fraction_on_source_dark_samples']=float(np.mean(photo_vals)) if photo_vals else 0.0
    c['trace_id']=f'DAX{i:02d}'
# CSV
cols=['trace_id','angle_deg','rho','span','coverage','pieces','x1','y1','x2','y2','source_dark_fraction','source_dark_samples','bsb_support_fraction_on_source_dark_samples','status']
with (B/'selling_dotted_axis_traces_v45.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
    for c in uniq:
        r={k:(round(c[k],9) if isinstance(c.get(k),float) else c.get(k)) for k in cols if k!='status'}
        r['status']='DUAL_SOURCE_SUPPORTED_DOTTED_TRACE_DIAGNOSTIC_NOT_BOUND_TO_HAX'
        w.writerow(r)
# Overlay on canonical raster.
ov=cv2.cvtColor(im,cv2.COLOR_GRAY2BGR)
for i,c in enumerate(uniq):
    a=(int(round(c['x1'])),int(round(c['y1'])));b=(int(round(c['x2'])),int(round(c['y2'])))
    cv2.line(ov,a,b,(0,0,255),2,cv2.LINE_AA); cv2.putText(ov,c['trace_id'],a,cv2.FONT_HERSHEY_SIMPLEX,.45,(0,0,255),1,cv2.LINE_AA)
cv2.imwrite(str(B/'selling_dotted_axis_traces_overlay_v45.png'),ov)
# Check endpoint intersections among segment extents to expose non-uniqueness.
def seg_inter(c,d,tol=1e-9):
    p=np.array([c['x1'],c['y1']]); r=np.array([c['x2']-c['x1'],c['y2']-c['y1']]); q=np.array([d['x1'],d['y1']]); s=np.array([d['x2']-d['x1'],d['y2']-d['y1']])
    cr=lambda a,b:a[0]*b[1]-a[1]*b[0]; den=cr(r,s)
    if abs(den)<tol:return None
    t=cr(q-p,s)/den;u=cr(q-p,r)/den
    if -1e-6<=t<=1+1e-6 and -1e-6<=u<=1+1e-6:return [float(p[0]+t*r[0]),float(p[1]+t*r[1])]
    return None
ints=[]
for i in range(len(uniq)):
    for j in range(i+1,len(uniq)):
        x=seg_inter(uniq[i],uniq[j])
        if x is not None:ints.append({'a':uniq[i]['trace_id'],'b':uniq[j]['trace_id'],'point_px300':[round(x[0],3),round(x[1],3)]})
out={
 'version':'v45','source_gate':'Selling source semantics assert six external bipartite symmetry axes and six fourfold centres; source text also identifies dotted traces as symmetry-axis notation.',
 'detector':{'base':'v37 LSD collinear short-segment grouping','v45_diagnostic_changes':{'min_span_px300':250,'coverage_window':[0.04,0.30]},'historical_promotion':'FORBIDDEN_FROM_DETECTOR_ALONE'},
 'registration':{'witness':'selling_fig1_bsb_photo_gate_v38.json','p95_reprojection_error_px300':pg['registration_to_v37']['p95_reprojection_error_px_v37_300']},
 'trace_count':len(uniq),'trace_ids':[c['trace_id'] for c in uniq],
 'min_bsb_dash_support':min((c['bsb_support_fraction_on_source_dark_samples'] for c in uniq),default=0.0),
 'segment_extent_intersections':ints,
 'historical_external_axis_count':6,
 'historical_binding':{'HAX_bound':0,'HC_bound':0,'reason':'Seven robust dual-source-supported dotted traces are exposed, whereas the source requires six external axes; current unlabeled geometry does not canonically select the external six or their six-centre cyclic ordering.'},
 'C2_89':{'activated':False,'reason':'No HAX_i reflection is source-bound; node action cannot be certified.'},
 'decision':'PROVEN_SOURCE_DOTTED_TRACE_SET_7_DIAGNOSTIC; HISTORICAL_EXTERNAL_SIX_SUBSET_NONUNIQUE; HAX_HC_BINDING_NT_FAIL_CLOSED',
 'guard':'Seven trace candidates are not seven historical external axes. Cardinality/visual fit cannot select six without source labels/incidence.'
}
(B/'selling_bsb_axis_binding_attempt_v45.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'trace_count':len(uniq),'min_bsb_support':out['min_bsb_dash_support'],'intersections':len(ints),'decision':out['decision']},indent=2))
