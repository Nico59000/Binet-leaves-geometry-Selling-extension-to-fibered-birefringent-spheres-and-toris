import json,time,sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form
B='/mnt/data/v50_work'; D=json.load(open(B+'/selling_integral_orbifold_cochains_v50.json'))
def M(rows,cols,E): return sp.SparseMatrix(rows,cols,{(e['row'],e['col']):int(e['value']) for e in E})
for sys in ['Binet','Ad']:
 d=D['systems'][sys]; c=d['dimensions'];
 for nm in ['d0','d1']:
  A=M(c['C1'] if nm=='d0' else c['C2'],c['C0'] if nm=='d0' else c['C1'],d[nm+'_sparse'])
  print(sys,nm,A.shape,'start',flush=True); t=time.time(); S=smith_normal_form(A,domain=ZZ); dt=time.time()-t
  diag=[abs(int(S[i,i])) for i in range(min(S.rows,S.cols)) if S[i,i]!=0]
  print('done',dt,'rank',len(diag),'nonunit',[(x,diag.count(x)) for x in sorted(set(diag)) if x!=1],flush=True)
