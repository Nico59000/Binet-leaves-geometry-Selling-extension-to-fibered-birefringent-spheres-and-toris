#!/usr/bin/env python3
from __future__ import annotations
import json,sympy as sp,collections,time
from pathlib import Path
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
B=Path('/mnt/data/v51_work')
g3=json.load(open(B/'gamma3_level2_relative_3cells_v36.json'));cells=g3['cells'];gens=g3['generator_order'];qcyc=json.load(open(B/'gamma3_radius5_quotient_cycle_v51.json'))
STD={}
for i in range(1,4):F=sp.eye(3);F[i-1,i-1]=-1;STD[f'F{i}']=F
for i in range(1,4):
 for j in range(1,4):
  if i!=j:E=sp.eye(3);E[i-1,j-1]=2;STD[f'E{i}{j}']=E
def mt(M):return tuple(tuple(int(M[i,j]) for j in range(3)) for i in range(3))
def mm(A,C):return tuple(tuple(A[i][0]*C[0][j]+A[i][1]*C[1][j]+A[i][2]*C[2][j] for j in range(3)) for i in range(3))
def mi(A):return mt(sp.Matrix(A).inv())
def mod4(A):return tuple(tuple(x%4 for x in row) for row in A)
I=mt(sp.eye(3));MAT={g:mt(M) for g,M in STD.items()};MATI={g:mi(MAT[g]) for g in MAT};letters=[]
for g in gens:
 letters.append((g,1,MAT[g]));inv=MATI[g]
 if inv!=MAT[g]:letters.append((g,-1,inv))
ball={I};front=[I]
for d in range(1,6):
 nf=[]
 for a in front:
  for g,s,M in letters:
   h=mm(a,M)
   if h not in ball:ball.add(h);nf.append(h)
 front=nf
fibers=collections.defaultdict(list)
for A in ball:fibers[mod4(A)].append(A)
fox=[]
for c in cells:
 out={x:{} for x in gens};pref=I
 for tok in c['word']:
  inv=tok.endswith('^-1');b=tok[:-3] if inv else tok
  if inv:q=MATI[b];key=mm(pref,q);out[b][key]=out[b].get(key,0)-1;pref=mm(pref,q)
  else:key=pref;out[b][key]=out[b].get(key,0)+1;pref=mm(pref,MAT[b])
 assert pref==I
 fox.append({x:[(q,v,mi(q)) for q,v in dd.items() if v] for x,dd in out.items()})
rows=[];b=[]
for t in qcyc['actual_boundary_residual']:
 rows.append((t['generator'],tuple(tuple(int(x) for x in rr) for rr in t['group_matrix'])));b.append(-int(t['coefficient']))
ri={r:i for i,r in enumerate(rows)}
fset=set()
for x,G in rows:
 for ci in range(43):
  for q,c,qinv in fox[ci][x]:
   A=mm(G,qinv)
   if A in ball:fset.add((ci,mod4(A)))
def proj(ci,A):
 d={}
 for x in gens:
  for q,c,qinv in fox[ci][x]:
   j=ri.get((x,mm(A,q)))
   if j is not None:d[j]=d.get(j,0)+c
 return {j:v for j,v in d.items() if v}
cols=[];meta=[]
for ci,k in sorted(fset,key=lambda z:(z[0],repr(z[1]))):
 vals=fibers[k];r=proj(ci,vals[0])
 for A in vals[1:]:
  q=proj(ci,A);d=dict(q)
  for j,v in r.items():d[j]=d.get(j,0)-v
  d={j:v for j,v in d.items() if v}
  if d:cols.append(d);meta.append((ci,k,vals[0],A))
print('matrix',96,len(cols),flush=True)
M=sp.zeros(96,len(cols))
for j,d in enumerate(cols):
 for i,v in d.items():M[i,j]=v
print('rankQ',M.rank(), 'aug',M.row_join(sp.Matrix(b)).rank(),flush=True)
t=time.time();S=smith_normal_form(M,domain=ZZ);print('snf sec',time.time()-t,flush=True)
diag=[abs(int(S[i,i])) for i in range(min(S.rows,S.cols)) if S[i,i]!=0]
print('diag nontrivial',[x for x in diag if x!=1][:30],'rank',len(diag),flush=True)
# Solve rational, use gauss_jordan_solve maybe
sol=sp.linsolve((M,sp.Matrix(b)))
vec=next(iter(sol)); # substitute all parameters 0
params=set().union(*(e.free_symbols for e in vec)); subs={s:0 for s in params};v0=[sp.simplify(e.subs(subs)) for e in vec]
den=sp.ilcm(*[sp.denom(x) for x in v0]);nz=[(i,x) for i,x in enumerate(v0) if x]
print('solution den',den,'nz',len(nz),'maxnum',max(abs(int(x*den)) for i,x in nz) if nz else 0)
out={'version':'v51','status':'PROVEN-EXACT-RESIDUAL-PROJECTION-LATTICE','matrix_shape':[96,len(cols)],'rank_Q':M.rank(),'rank_aug_Q':M.row_join(sp.Matrix(b)).rank(),'smith_nonzero_count':len(diag),'smith_nonunit_factors':[x for x in diag if x!=1],'residual_projection_has_integral_solution':all(x==1 for x in diag), 'one_rational_solution_common_denominator':int(den),'one_rational_solution_nonzero_directions':len(nz),'guard':'This solves only the 96 original residual coordinates. The selected correction generally creates side residuals outside this coordinate set and is not yet a full group-ring cycle.'}
(B/'gamma3_radius5_residual_projection_lattice_v51.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,indent=2))
