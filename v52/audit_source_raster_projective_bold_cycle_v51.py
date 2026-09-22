#!/usr/bin/env python3
import csv,json,itertools,numpy as np,pathlib,hashlib
B=pathlib.Path(__file__).resolve().parent
gate=json.load(open(B/'selling_raster_central_bold_cycle_gate_v51.json'));order=gate['cycle_order']
nodes={r['node_id']:(float(r['x_px300']),float(r['y_px300'])) for r in csv.DictReader(open(B/'selling_fig1_nodes_v40.csv'))};pts=np.array([nodes[n] for n in order],float)
src=np.array([[1,0],[1,.25],[2/3,1/3],[0,1/3],[-2/3,1/9],[-2/3,0]],float);scale=np.linalg.norm(pts.max(0)-pts.min(0))
def hom4(X,Y):
 A=[];b=[]
 for (x,y),(u,v) in zip(X,Y):
  A.append([x,y,1,0,0,0,-u*x,-u*y]);b.append(u)
  A.append([0,0,0,x,y,1,-v*x,-v*y]);b.append(v)
 try:h=np.linalg.solve(np.array(A),np.array(b))
 except np.linalg.LinAlgError:return None
 return np.r_[h,1].reshape(3,3)
def app(H,X):
 Z=np.c_[X,np.ones(len(X))]@H.T
 if np.any(abs(Z[:,2])<1e-10):return None
 return Z[:,:2]/Z[:,2,None]
res=[]
for inds in itertools.combinations(range(16),6):
 Y0=pts[list(inds)]
 for rev in [False,True]:
  Y=Y0[::-1] if rev else Y0
  for sh in range(6):
   X=np.roll(src,sh,axis=0);H=hom4(X[:4],Y[:4])
   if H is None:continue
   pred=app(H,X[4:])
   if pred is None or not np.all(np.isfinite(pred)):continue
   e=np.linalg.norm(pred-Y[4:],axis=1);res.append((float(np.sqrt(np.mean(e*e))),float(max(e)),inds,rev,sh))
res.sort();b=res[0]
out={'version':'v51','status':'PROVEN-DIAGNOSTIC-PROJECTIVE-BOLD-CYCLE-AUDIT','model':'necessary six-point global homography test; first four ordered source centres determine H, last two are held out','candidate_count':len(res),'source_centres':src.tolist(),'raster_cycle_order':order,'best':{'heldout_rms_px300':b[0],'heldout_max_px300':b[1],'normalized_rms_by_cycle_bbox_diagonal':b[0]/scale,'indices':list(b[2]),'nodes':[order[i] for i in b[2]],'reversed':b[3],'source_cyclic_shift':b[4]},'second_best_normalized_rms':res[1][0]/scale,'decision':'DIAGNOSTIC-ONLY/NO-HISTORICAL-PROMOTION','guard':'Selling explicitly permits a bent/curved representation of the symmetry axes. A global homography on the reinforced 16-cycle is not a source invariant; this test can neither certify nor refute the historical HC/HAX binding.'}
p=B/'selling_source_raster_projective_bold_cycle_audit_v51.json';p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,indent=2));print('sha256',hashlib.sha256(p.read_bytes()).hexdigest())
