#!/usr/bin/env python3
from __future__ import annotations
import json,sympy as sp,collections,time,hashlib
from pathlib import Path
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
for k in fibers:fibers[k].sort(key=repr)
fox=[]
for c in cells:
 out={x:{} for x in gens};pref=I
 for tok in c['word']:
  inv=tok.endswith('^-1');b=tok[:-3] if inv else tok
  if inv:q=MATI[b];key=mm(pref,q);out[b][key]=out[b].get(key,0)-1;pref=mm(pref,q)
  else:key=pref;out[b][key]=out[b].get(key,0)+1;pref=mm(pref,MAT[b])
 assert pref==I
 fox.append({x:[(q,v,mi(q)) for q,v in dd.items() if v] for x,dd in out.items()})

def boundary_col(ci,A):
 d={}
 for x in gens:
  for q,c,qinv in fox[ci][x]:
   r=(x,mm(A,q));d[r]=d.get(r,0)+c
 return {r:v for r,v in d.items() if v}
def addto(dst,src,a=1):
 for r,v in src.items():
  x=dst.get(r,0)+a*v
  if x:dst[r]=x
  elif r in dst:del dst[r]
# initial boundary of multiplied-by-4 witness from stored residual
R={}
for t in qcyc['actual_boundary_residual']:
 r=(t['generator'],tuple(tuple(int(x) for x in rr) for rr in t['group_matrix']));R[r]=int(t['coefficient'])
chain=collections.Counter() # correction only: (ci,A)->integer coeff, quotient kernel by construction
history=[]
seen_residual={}
for step in range(1,17):
 rows=sorted(R,key=repr);
 # canonical residual signature for exact periodicity detection
 sig=hashlib.sha256(repr([(r,R[r]) for r in rows]).encode()).hexdigest()
 if sig in seen_residual:
  history.append({'step':step,'status':'EXACT-RESIDUAL-PERIOD-DETECTED','repeats_step':seen_residual[sig],'residual_count':len(rows),'residual_sha256':sig})
  break
 seen_residual[sig]=step
 ri={r:i for i,r in enumerate(rows)};b=[-R[r] for r in rows]
 # find fibers touching current residual rows
 fset=set()
 for x,G in rows:
  for ci in range(43):
   for q,c,qinv in fox[ci][x]:
    A=mm(G,qinv)
    if A in ball:fset.add((ci,mod4(A)))
 # collect distinct nonzero projected directions, but keep metadata first occurrence
 cols=[];meta=[];seen={}
 def proj(ci,A):
  d={}
  for x in gens:
   for q,c,qinv in fox[ci][x]:
    j=ri.get((x,mm(A,q)))
    if j is not None:d[j]=d.get(j,0)+c
  return {j:v for j,v in d.items() if v}
 for ci,k in sorted(fset,key=lambda z:(z[0],repr(z[1]))):
  vals=fibers[k];ref=vals[0];pr=proj(ci,ref)
  for A in vals[1:]:
   pa=proj(ci,A);d=dict(pa)
   for j,v in pr.items():d[j]=d.get(j,0)-v
   d={j:v for j,v in d.items() if v}
   if not d:continue
   key=tuple(sorted(d.items()))
   if key in seen:continue
   seen[key]=len(cols);cols.append(d);meta.append((ci,k,ref,A))
 print('step',step,'res',len(rows),'fibers',len(fset),'projdirs',len(cols),flush=True)
 if not cols:
  history.append({'step':step,'residual_count_before':len(rows),'status':'STUCK_NO_PROJECTED_DIRECTIONS'});break
 M=sp.zeros(len(rows),len(cols))
 for j,d in enumerate(cols):
  for i,v in d.items():M[i,j]=v
 rb=M.rank();ra=M.row_join(sp.Matrix(b)).rank()
 print('ranks',rb,ra,flush=True)
 if rb!=ra:
  history.append({'step':step,'residual_count_before':len(rows),'fiber_count':len(fset),'projected_direction_count':len(cols),'rank_Q':rb,'rank_aug_Q':ra,'status':'STUCK_RATIONAL_INCONSISTENT'});break
 sol=next(iter(sp.linsolve((M,sp.Matrix(b)))))
 params=set().union(*(e.free_symbols for e in sol));sub={s:0 for s in params};v=[sp.simplify(e.subs(sub)) for e in sol]
 den=sp.ilcm(*[sp.denom(x) for x in v]);nz=[(j,x) for j,x in enumerate(v) if x]
 print('sol den',den,'nz',len(nz),flush=True)
 if den!=1:
  history.append({'step':step,'residual_count_before':len(rows),'fiber_count':len(fset),'projected_direction_count':len(cols),'rank_Q':rb,'rank_aug_Q':ra,'solution_denominator':int(den),'solution_nz':len(nz),'status':'STUCK_NONINTEGRAL_CANONICAL_SOLUTION'});break
 # apply exact correction directions: x*(col(A)-col(ref)) to full boundary and chain
 step_terms=[]
 for j,x in nz:
  coeff=int(x);ci,k,ref,A=meta[j]
  chain[(ci,A)]+=coeff;chain[(ci,ref)]-=coeff
  addto(R,boundary_col(ci,A),coeff);addto(R,boundary_col(ci,ref),-coeff)
  step_terms.append({'relator':cells[ci]['name'],'fiber_mod4':[list(r) for r in k],'A':[list(r) for r in A],'ref':[list(r) for r in ref],'coefficient':coeff})
 # remove zero chain entries
 chain += collections.Counter()
 hist={'step':step,'residual_sha256_before':sig,'residual_count_before':len(rows),'fiber_count':len(fset),'projected_direction_count':len(cols),'rank_Q':rb,'rank_aug_Q':ra,'solution_denominator':1,'solution_nz':len(nz),'residual_count_after':len(R),'residual_l1_after':sum(abs(v) for v in R.values()),'correction_terms':step_terms,'status':'APPLIED'};history.append(hist)
 print('after',len(R),'L1',sum(abs(v) for v in R.values()),'chain nz',sum(1 for v in chain.values() if v),flush=True)
 if not R:break
out={'version':'v51','status':'PROVEN-EXACT-ITERATED-KERNEL-CORRECTION-'+('CLOSED' if not R else 'CHECKPOINT'),'fixed_quotient_witness':'gamma3_radius5_quotient_cycle_v51.json multiplied by common denominator 4','steps':history,'final_residual_count':len(R),'final_residual_l1':sum(abs(v) for v in R.values()),'kernel_correction_chain_nonzero_terms':sum(1 for v in chain.values() if v),'integral_lift_found':not R,'residual_period_detected':any(h.get('status')=='EXACT-RESIDUAL-PERIOD-DETECTED' for h in history),'guard':'Every correction is a within-(relator,mod4-class) difference, so the finite quotient 3-chain is unchanged. Failure/stopping of this deterministic elimination is not a global no-go unless an explicit inconsistency certificate is recorded.'}
# sample final residual, chain
out['final_residual_sample']=[{'generator':x,'group_matrix':[list(rr) for rr in A],'coefficient':c} for (x,A),c in list(sorted(R.items(),key=lambda z:(z[0][0],repr(z[0][1]))))[:200]]
out['kernel_chain_sample']=[{'relator':cells[ci]['name'],'group_matrix':[list(rr) for rr in A],'coefficient':c} for (ci,A),c in list(sorted(((k,v) for k,v in chain.items() if v),key=lambda z:(z[0][0],repr(z[0][1]))))[:200]]
p=B/'gamma3_radius5_iterated_kernel_correction_v51.json';p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({k:out[k] for k in ['status','final_residual_count','final_residual_l1','kernel_correction_chain_nonzero_terms','integral_lift_found']},indent=2));print('sha256',hashlib.sha256(p.read_bytes()).hexdigest())
