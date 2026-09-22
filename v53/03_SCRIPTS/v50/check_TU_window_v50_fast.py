#!/usr/bin/env python3
import json,sys,hashlib
from pathlib import Path
HERE=Path('/mnt/data/v50_work');sys.path.insert(0,str(HERE))
import selling_global_wall_linearization_v50 as gl
from selling_reconstruction_engine_v49 import RELABEL

def tt(M): return tuple(tuple(int(M[i,j]) for j in range(3)) for i in range(3))
def mul(A,B):
 return tuple(tuple(A[i][0]*B[0][j]+A[i][1]*B[1][j]+A[i][2]*B[2][j] for j in range(3)) for i in range(3))
def neg(A): return tuple(tuple(-x for x in row) for row in A)
def inv3(A):
 import sympy as sp
 M=sp.Matrix(A).inv();return tt(M)
I=((1,0,0),(0,1,0),(0,0,1))
GT=tt(gl.GT);GU=tt(gl.GU);RBS=[tt(B) for B in RELABEL]
def powers(M,N):
 d={0:I}
 for i in range(1,N+1): d[i]=mul(d[i-1],M)
 Mi=inv3(M)
 for i in range(1,N+1): d[-i]=mul(d[-i+1],Mi)
 return d
D=json.load(open(HERE/'exact_TU_field_bfs_v50.json'))
reps=[tuple(tuple(int(x) for x in row) for row in r['path_matrix']) for r in D['records']]
results=[]
for N in range(3,9):
 PT=powers(GT,N); PU=powers(GU,N); decks=[mul(A,B) for A in PT.values() for B in PU.values()]
 seen={};pairs=set()
 for i,A in enumerate(reps):
  for De in decks:
   DA=mul(De,A)
   for B in RBS:
    X=mul(DA,B)
    for K in (X,neg(X)):
     j=seen.get(K)
     if j is None: seen[K]=i
     elif j!=i: pairs.add(tuple(sorted((i,j))))
 print('bound',N,'collisions',len(pairs),'keys',len(seen),flush=True)
 results.append({'bound':N,'cross_class_collisions':len(pairs),'pairs':sorted(pairs)[:200],'key_count':len(seen)})
out={'version':'v50','field_count':len(reps),'status':'PASS_WINDOW_3_TO_8' if all(r['cross_class_collisions']==0 for r in results) else 'REFUTED_BY_LARGER_WINDOW','results':results}
(HERE/'selling_TU_window_stability_v50.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
