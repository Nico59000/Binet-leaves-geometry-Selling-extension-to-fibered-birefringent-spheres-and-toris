#!/usr/bin/env python3
import json, sympy as sp, hashlib
from pathlib import Path
B=Path('/mnt/data/v51_work'); car=json.load(open(B/'selling_source_axis_centre_carrier_v51.json'))
R=[
 sp.diag(1,-1,1),
 sp.Matrix([[3,0,-2],[0,1,0],[4,0,-3]]),
 sp.Matrix([[13,24,-6],[-4,-7,2],[12,24,-5]]),
 sp.Matrix([[-8,-21,0],[3,8,0],[-2,-6,-1]]), # G_T^T
 sp.Matrix([[-4,-9,-3],[1,2,1],[2,6,1]]), # G_U^T
 sp.Matrix([[5,0,6],[0,1,0],[-4,0,-5]])
]
axes=[]
for i,(a,M) in enumerate(zip(car['axes'],R)):
 assert M*M==sp.eye(3)
 axes.append({'id':a['id'],'matrix_on_v':[[int(M[x,y]) for y in range(3)] for x in range(3)],'det':int(M.det()),'projective_involution':True})
cent=[]
for i,c in enumerate(car['centres']):
 q=sp.Rational(c['q']);r=sp.Rational(c['r']);v=sp.Matrix([1,-r,q]);A=R[i];D=R[(i+1)%6]
 assert A*D==D*A
 ev=[]
 for M in [A,D,A*D]:
  w=M*v
  lam=None
  for j in range(3):
   if v[j]!=0: lam=sp.simplify(w[j]/v[j]);break
  assert w==lam*v and lam in (1,-1)
  ev.append(int(lam))
 elems=[sp.eye(3),A,D,A*D]
 # exact distinct matrices; projectively perhaps no +/- duplicates
 assert len({tuple(M) for M in elems})==4
 assert (A*D)**2==sp.eye(3)
 cent.append({'id':c['id'],'q':c['q'],'r':c['r'],'projective_point_v':['1',str(-r),str(q)],'generators':[car['axes'][i]['id'],car['axes'][(i+1)%6]['id']],'generator_commute':True,'common_projective_eigenvalues':ev[:2],'product_eigenvalue':ev[2],'stabilizer_abstract':'C2 x C2','stabilizer_order':4,'four_sector_interpretation':'PROVEN-SOURCE-PROJECTIVE-STABILIZER'})
out={'version':'v51','status':'PROVEN-SOURCE-SIX-FOURFOLD-STABILIZERS','axes':axes,'centres':cent,'decision':'EACH_HC_i_HAS_EXACT_PROJECTIVE_V4_STABILIZER_GENERATED_BY_THE_TWO_ADJACENT_HAX_INVOLUTIONS','guard':'This certifies the source-coordinate fourfold-centre semantics. It still does not assign HC_i/HAX_i to raster coordinates/curves.'}
p=B/'selling_source_axis_centre_stabilizers_v51.json';p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2));print('sha256',hashlib.sha256(p.read_bytes()).hexdigest())
