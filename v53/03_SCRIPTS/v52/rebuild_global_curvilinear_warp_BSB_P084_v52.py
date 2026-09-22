#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict,Counter
import json,csv,ast,math
import numpy as np
from scipy.interpolate import RBFInterpolator
from shapely.geometry import LineString,Polygon,Point
ROOT=Path(__file__).resolve().parent;SRC=ROOT/'deps'
raster=json.load(open(SRC/'selling_fig1_raster_C6_HAX_closure_v52.json'))
geo=json.load(open(SRC/'selling_regular_62face_raster_transport_visualization_v52.geojson'))
ind=json.load(open(ROOT/'selling_curvilinear_edge_support_audit_v52.json'))
seg=list(csv.DictReader(open(SRC/'selling_fig1_segments_v37.csv',encoding='utf-8')));segmap={r['segment_id']:r for r in seg}
anchors={a['id']:np.array(a['xy_px300'],float) for a in raster['six_anchors']}; rhax={a['id']:a for a in raster['six_raster_axis_arcs']}
# historical raster path for each side
def segpts(sid): return np.array(ast.literal_eval(segmap[sid]['polyline_px300']),float)
def assemble(rid):
 a=rhax[rid];cur=anchors[a['from']].copy();out=[cur.copy()]
 for piece in a['path']:
  if piece.startswith('S'):
   P=segpts(piece)
   if np.linalg.norm(P[-1]-cur)<np.linalg.norm(P[0]-cur):P=P[::-1]
   if np.linalg.norm(P[0]-cur)>1e-6:out.append(P[0])
   out.extend(P[1:]);cur=P[-1]
 end=anchors[a['to']]
 if np.linalg.norm(np.array(out[-1])-end)>1e-6:out.append(end)
 return np.array(out,float)
paths={rid:assemble(rid) for rid in rhax}
def resample(P,n):
 ls=LineString(P);L=ls.length
 return np.array([ls.interpolate(s).coords[0] for s in np.linspace(0,L,n)],float)
# Current unique PL edge chords and endpoints
edges={}; endpoints={}
for f in geo['features']:
 if f['properties'].get('kind')=='EDGE_OCCURRENCE':
  eid=int(f['properties']['edge_id']);edges.setdefault(eid,np.array(f['geometry']['coordinates'],float)); endpoints[eid]=(int(f['properties']['v_start']),int(f['properties']['v_end']))
assert len(edges)==123
# Control pairs x(domain PL) -> y(curved target)
X=[];Y=[];labels=[]
# Boundary: each RHAX_i maps its chord RHC_i->RHC_{i+1} to actual historical arc.
for i in range(6):
 rid=f'RHAX{i}'; a=rhax[rid]; A=anchors[a['from']];B=anchors[a['to']]
 n=19; dom=np.linspace(0,1,n)[:,None]*B+(1-np.linspace(0,1,n))[:,None]*A
 tgt=resample(paths[rid],n)
 for k in range(n): X.append(dom[k]);Y.append(tgt[k]);labels.append(('boundary',rid,k))
# Internal soft controls only from independent edge candidates with strong support and modest offsets.
# Use interior points only; endpoints are globally shared and stay controlled by the smooth field.
strong=[]
curves={int(k):np.array(v,float) for k,v in ind['curves_px300'].items()}
for r in ind['edges']:
 pre=r.get('pre_crossing_geometric_status',r['geometric_status'])
 if pre=='STRONG-CURVILINEAR-STROKE-SUPPORT' and r['max_abs_normal_offset_px']<=18:
  strong.append(r['edge_id'])
for eid in strong:
 P=edges[eid]; A=P[0];B=P[-1]; Q=curves[eid]
 # 5 interior correspondences at equal arclength, but only if proposed displacement <=18 px.
 qrs=resample(Q,7)
 for k,t in enumerate(np.linspace(0,1,7)[1:-1],start=1):
  x=A*(1-t)+B*t; y=qrs[k]
  if np.linalg.norm(y-x)<=18:
   X.append(x);Y.append(y);labels.append(('internal-strong',eid,k))
# Additional direct BSB constraint for residual edge 96 from independently extracted photo candidate P084.
# Deterministic: project 13 interior prior-warp samples onto P084 and retain only distance <= 24 px.
# First reconstruct the previous global warp using controls accumulated so far.
_d0=defaultdict(list)
for x,y in zip(X,Y): _d0[tuple(np.round(x,5))].append(y)
_X0=np.array(list(_d0.keys()),float); _Y0=np.array([np.mean(_d0[k],axis=0) for k in _d0],float)
_r0=RBFInterpolator(_X0,_Y0-_X0,kernel='thin_plate_spline',degree=1,smoothing=500)
# photo candidate P084 in v37 frame
_pg0=json.load(open(SRC/'selling_fig1_bsb_photo_gate_v38.json')); _Hm=np.array(_pg0['registration_to_v37']['H_photo_to_v37_300'])
_ph=json.load(open(ROOT/'selling_fig1_photo_candidates_pre_v38.geojson'))
_ff=next(f for f in _ph['features'] if f['properties'].get('candidate_id')=='P084')
_PP=np.array(_ff['geometry']['coordinates'],float); _zz=(_Hm@np.c_[_PP,np.ones(len(_PP))].T).T; _zz=_zz[:,:2]/_zz[:,2,None]
_Lp=LineString(_zz)
_A,_B=edges[96][0],edges[96][-1]
for kk,t in enumerate(np.linspace(.08,.92,13)):
    xx=_A*(1-t)+_B*t; qq=xx+_r0(xx[None,:])[0]; yy=np.array(_Lp.interpolate(_Lp.project(Point(qq))).coords[0])
    if np.linalg.norm(yy-qq)<=24:
        X.append(xx);Y.append(yy);labels.append(('BSB-P084-edge96',96,kk))

# Add polygon centroid fixed-ish to suppress global drift; target=current.
poly=Polygon([anchors[f'RHC{i}'] for i in range(6)])
c=np.array(poly.centroid.coords[0]);X.append(c);Y.append(c);labels.append(('gauge-centroid',0,0))
# Deduplicate near-exact control domain points; average targets (boundary corners repeated consistently).
# Quantize 1e-5 px.
d=defaultdict(list);dl=defaultdict(list)
for x,y,l in zip(X,Y,labels):d[tuple(np.round(x,5))].append(y);dl[tuple(np.round(x,5))].append(l)
X=np.array(list(d.keys()),float);Y=np.array([np.mean(d[k],axis=0) for k in d],float)
print('controls',len(X),'strong edges',strong)
# Baseline crossing set
base_cross=set()
L0={i:LineString(edges[i]) for i in edges}
for i in range(123):
 for j in range(i+1,123):
  if set(endpoints[i])&set(endpoints[j]):continue
  inter=L0[i].intersection(L0[j])
  if not inter.is_empty:base_cross.add((i,j))
print('baseline crossings',base_cross)
# evaluate smoothing grid on displacement field, F(x)=x+RBF(x)
D=Y-X
# grid inside polygon
minx,miny,maxx,maxy=poly.bounds
grid=[]
for x in np.linspace(minx+2,maxx-2,52):
 for y in np.linspace(miny+2,maxy-2,62):
  if poly.contains(Point(x,y)):grid.append((x,y))
grid=np.array(grid,float)

def map_edge(rbf,P,n=80):
 A=P[0];B=P[-1]; t=np.linspace(0,1,max(3,int(np.linalg.norm(B-A)/2)+1)); Xq=A[None,:]*(1-t[:,None])+B[None,:]*t[:,None]
 return Xq+rbf(Xq)

def eval_sm(sm):
 rbf=RBFInterpolator(X,D,kernel='thin_plate_spline',degree=1,smoothing=sm)
 # control residual
 pred=X+rbf(X);cres=np.linalg.norm(pred-Y,axis=1)
 # jacobian F via finite difference
 h=.25
 Qxp=(grid+np.array([h,0]))+rbf(grid+np.array([h,0]));Qxm=(grid-np.array([h,0]))+rbf(grid-np.array([h,0]));
 Qyp=(grid+np.array([0,h]))+rbf(grid+np.array([0,h]));Qym=(grid-np.array([0,h]))+rbf(grid-np.array([0,h]));
 dx=(Qxp-Qxm)/(2*h);dy=(Qyp-Qym)/(2*h);det=dx[:,0]*dy[:,1]-dx[:,1]*dy[:,0]
 # mapped edge crossings
 L={i:LineString(map_edge(rbf,edges[i])) for i in range(123)}
 cr=set()
 for i in range(123):
  for j in range(i+1,123):
   if set(endpoints[i])&set(endpoints[j]):continue
   if not L[i].intersection(L[j]).is_empty:cr.add((i,j))
 return rbf,{'smoothing':sm,'control_residual_max':float(cres.max()),'control_residual_p95':float(np.quantile(cres,.95)),'jac_min':float(det.min()),'jac_p01':float(np.quantile(det,.01)),'jac_negative_or_zero':int(np.sum(det<=0)),'jac_samples':len(det),'crossing_count':len(cr),'new_crossings':sorted(cr-base_cross),'lost_baseline_crossings':sorted(base_cross-cr)}
fit_results=[]
for sm in [0,1e-4,1e-3,1e-2,.05,.1,.2,.5,1,2,5,10,20,50,100,200,500,1000]:
 try:
  r,m=eval_sm(sm);fit_results.append((r,m));print(m)
 except Exception as ex:print('FAIL',sm,ex)
# admissible: positive Jacobian everywhere and no new crossings. Prefer smallest boundary/control p95.
adm=[x for x in fit_results if x[1]['jac_negative_or_zero']==0 and not x[1]['new_crossings']]
if not adm:
 print('NO ADMISSIBLE GLOBAL WARP')
 sel=min(fit_results,key=lambda x:(x[1]['jac_negative_or_zero'],len(x[1]['new_crossings']),x[1]['control_residual_p95']))
else:sel=min(adm,key=lambda x:x[1]['control_residual_p95'])
rbf,metric=sel;print('SELECT',metric)
# serialize mapped curves and support against source photo geometry distance fields reused from independent audit by reconstructing masks
import cv2
from scipy.ndimage import distance_transform_edt
pg=json.load(open(SRC/'selling_fig1_bsb_photo_gate_v38.json'));Hmat=np.array(pg['registration_to_v37']['H_photo_to_v37_300']);photo=json.load(open(ROOT/'selling_fig1_photo_candidates_pre_v38.geojson'))
mask=np.zeros((2700,2225),np.uint8)
for f in photo['features']:
 if f['properties']['class'] not in ('PHOTO_MATCHED_TO_V37_GEOMETRY','PHOTO_NEW_LONG_GEOMETRY_CANDIDATE'):continue
 P=np.array(f['geometry']['coordinates'],float);q=(Hmat@np.c_[P,np.ones(len(P))].T).T;q=q[:,:2]/q[:,2,None]
 cv2.polylines(mask,[np.round(q).astype(np.int32).reshape(-1,1,2)],False,255,2)
dtg=distance_transform_edt(mask==0);gray=cv2.imread(str(ROOT/'bsb_warp_v37_300.png'),0);dtd=distance_transform_edt(~(gray<210));eps=float(pg['registration_to_v37']['p95_reprojection_error_px_v37_300'])
def vals(P,D):
 x=np.clip(np.rint(P[:,0]).astype(int),0,D.shape[1]-1);y=np.clip(np.rint(P[:,1]).astype(int),0,D.shape[0]-1);return D[y,x]
edgeout=[];curveout={}
for i in range(123):
 Q=map_edge(rbf,edges[i]);curveout[str(i)]=Q.tolist();vg=vals(Q,dtg);vd=vals(Q,dtd); sg=float(np.mean(vg<=eps));sd=float(np.mean(vd<=eps))
 if sg>=.75 and sd>=.75:st='STRONG-GLOBAL-CURVILINEAR-STROKE-SUPPORT'
 elif sg>=.35 and sd>=.50:st='PARTIAL-GLOBAL-CURVILINEAR-STROKE-SUPPORT'
 else:st='NT-GEOMETRIC-REALIZATION'
 edgeout.append({'edge_id':i,'endpoints':list(endpoints[i]),'geom_support_eps':sg,'dark_support_eps':sd,'geom_median_px':float(np.median(vg)),'geom_p95_px':float(np.quantile(vg,.95)),'status':st})
counts=Counter(e['status'] for e in edgeout);print('support',counts)
out={'version':'v52','phase':'next-immediate-gate','status':'PASS-GLOBAL-CURVILINEAR-HOMEOMORPHISM-CANDIDATE' if metric['jac_negative_or_zero']==0 and not metric['new_crossings'] else 'NT-GLOBAL-WARP-TOPOLOGY-NOT-CLOSED','domain':'current proven D6-selected PL raster realization','target_constraints':'six source-aware historical RHAX arcs + BSB-supported internal strong-edge soft controls','strong_internal_control_edges':strong,'control_count':len(X),'selected_metrics':metric,'all_fit_metrics':[m for _,m in fit_results],'baseline_nonincident_crossings':sorted(base_cross),'support_counts':dict(counts),'edges':edgeout,'curves_px300':curveout,'guard':'Only a single global positive-Jacobian/no-new-crossing warp may replace PL globally. Per-edge optimizers are evidence controls, not independent embeddings.'}
(ROOT/'selling_global_curvilinear_warp_BSB_P084_v52.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
