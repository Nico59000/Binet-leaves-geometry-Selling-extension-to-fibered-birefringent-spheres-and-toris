#!/usr/bin/env python3
import json,hashlib,collections
from pathlib import Path
import sympy as sp
B=Path(__file__).resolve().parent
g3=json.load(open(B/'gamma3_level2_relative_3cells_v36.json')); cells=g3['cells']; gens=g3['generator_order']; name_to_ci={c['name']:i for i,c in enumerate(cells)}
qcyc=json.load(open(B/'gamma3_radius5_quotient_cycle_v51.json')); cert=json.load(open(B/'gamma3_radius5_iterated_kernel_correction_v51.json'))
STD={}
for i in range(1,4):
 F=sp.eye(3);F[i-1,i-1]=-1;STD[f'F{i}']=F
for i in range(1,4):
 for j in range(1,4):
  if i!=j:
   E=sp.eye(3);E[i-1,j-1]=2;STD[f'E{i}{j}']=E
def mt(M):return tuple(tuple(int(M[i,j]) for j in range(3)) for i in range(3))
def mm(A,C):return tuple(tuple(A[i][0]*C[0][j]+A[i][1]*C[1][j]+A[i][2]*C[2][j] for j in range(3)) for i in range(3))
def mi(A):return mt(sp.Matrix(A).inv())
I=mt(sp.eye(3)); MAT={g:mt(M) for g,M in STD.items()}; MATI={g:mi(MAT[g]) for g in MAT}
fox=[]
for c in cells:
 out={x:{} for x in gens};pref=I
 for tok in c['word']:
  inv=tok.endswith('^-1');b=tok[:-3] if inv else tok
  if inv:
   q=MATI[b];key=mm(pref,q);out[b][key]=out[b].get(key,0)-1;pref=mm(pref,q)
  else:
   key=pref;out[b][key]=out[b].get(key,0)+1;pref=mm(pref,MAT[b])
 assert pref==I
 fox.append({x:[(q,v) for q,v in dd.items() if v] for x,dd in out.items()})
def boundary_col(ci,A):
 d={}
 for x in gens:
  for q,c in fox[ci][x]:
   k=(x,mm(A,q));d[k]=d.get(k,0)+c
 return {k:v for k,v in d.items() if v}
def add(dst,src,a=1):
 for k,v in src.items():
  z=dst.get(k,0)+a*v
  if z:dst[k]=z
  else:dst.pop(k,None)
def sig(R):
 rows=sorted(R,key=repr);return hashlib.sha256(repr([(r,R[r]) for r in rows]).encode()).hexdigest()
# reconstruct multiplied-by-4 canonical witness boundary
R={}
for t in qcyc['terms']:
 ci=int(t['relator_index']);A=tuple(tuple(int(x) for x in row) for row in t['integral_matrix_representative']);coef=int(t['coefficient_numerator']) # denom 4; multiplied by 4
 assert int(t['coefficient_denominator'])==4
 add(R,boundary_col(ci,A),coef)
stored={(t['generator'],tuple(tuple(int(x) for x in row) for row in t['group_matrix'])):int(t['coefficient']) for t in qcyc['actual_boundary_residual']}
assert R==stored and len(R)==96
hashes=[]; chain=collections.Counter()
for h in cert['steps']:
 if h['status']=='EXACT-RESIDUAL-PERIOD-DETECTED':
  assert sig(R)==h['residual_sha256']
  assert h['repeats_step']==2
  hashes.append(sig(R));break
 assert h['status']=='APPLIED'
 assert len(R)==h['residual_count_before'] and sig(R)==h['residual_sha256_before']
 for t in h['correction_terms']:
  ci=name_to_ci[t['relator']];A=tuple(tuple(int(x) for x in row) for row in t['A']);ref=tuple(tuple(int(x) for x in row) for row in t['ref']);coef=int(t['coefficient'])
  # exact within-fiber condition
  assert tuple(tuple(x%4 for x in row) for row in A)==tuple(tuple(x%4 for x in row) for row in ref)==tuple(tuple(int(x) for x in row) for row in t['fiber_mod4'])
  add(R,boundary_col(ci,A),coef);add(R,boundary_col(ci,ref),-coef)
  chain[(ci,A)]+=coef;chain[(ci,ref)]-=coef
 chain += collections.Counter()
 assert len(R)==h['residual_count_after'] and sum(abs(v) for v in R.values())==h['residual_l1_after']
 hashes.append(sig(R))
chain += collections.Counter()
assert len(R)==cert['final_residual_count']==96
assert sum(abs(v) for v in R.values())==cert['final_residual_l1']==96
assert sum(1 for v in chain.values() if v)==cert['kernel_correction_chain_nonzero_terms']==8
assert cert['integral_lift_found'] is False and cert['residual_period_detected'] is True
out={'version':'v51','status':'PASS-EXACT-LIGHTWEIGHT-ITERATED-CERTIFICATE-REPLAY','initial_residual_count':96,'step_or_period_hashes':hashes,'final_residual_count':len(R),'final_residual_l1':sum(abs(v) for v in R.values()),'kernel_chain_nonzero_terms':sum(1 for v in chain.values() if v),'period_repeats_step':2,'guard':'This replays the stored exact correction terms and their group-ring boundaries without re-solving the large local search problems.'}
p=B/'gamma3_radius5_iterated_kernel_replay_v51.json';p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2));print('sha256',hashlib.sha256(p.read_bytes()).hexdigest())
