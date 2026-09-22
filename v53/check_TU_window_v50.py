#!/usr/bin/env python3
import json,sys,hashlib
from pathlib import Path
import sympy as sp
HERE=Path('/mnt/data/v50_work');sys.path.insert(0,str(HERE))
from selling_global_wall_linearization_v50 import GT,GU
from selling_reconstruction_engine_v49 import RELABEL

def tup(A): return tuple(int(x) for x in list(A))
def mat(rows): return sp.Matrix([[sp.Integer(v) for v in row] for row in rows])
def powers(M,N):
 d={0:sp.eye(3)}
 for i in range(1,N+1): d[i]=d[i-1]*M
 Mi=M.inv()
 for i in range(1,N+1): d[-i]=d[-i+1]*Mi
 return d
D=json.load(open(HERE/'exact_TU_field_bfs_v50.json'))
reps=[mat(r['path_matrix']) for r in D['records']]
results=[]
for N in range(3,9):
 PT=powers(GT,N); PU=powers(GU,N)
 seen={}; collisions=[]
 for i,A in enumerate(reps):
  found=set()
  for a,Ta in PT.items():
   for b,Ub in PU.items():
    DA=Ta*Ub*A
    for bi,B in enumerate(RELABEL):
     for s in (1,-1):
      k=tup(s*DA*B)
      j=seen.get(k)
      if j is not None and j!=i: found.add(j)
      else: seen.setdefault(k,i)
  for j in sorted(found):
   if j<i: collisions.append((j,i))
 collisions=sorted(set(collisions))
 results.append({'bound':N,'cross_class_collisions':len(collisions),'pairs':collisions[:50],'key_count':len(seen)})
 print(N,len(collisions),collisions[:10],len(seen),flush=True)
out={'version':'v50','status':'PASS' if all(x['cross_class_collisions']==0 for x in results) else 'REFUTED','field_count':len(reps),'tested_bounds':[3,4,5,6,7,8],'results':results}
raw=json.dumps(out,indent=2).encode();out['sha256_pre_field']=hashlib.sha256(raw).hexdigest()
(HERE/'selling_TU_window_stability_v50.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
