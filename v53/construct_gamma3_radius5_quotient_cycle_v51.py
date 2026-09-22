#!/usr/bin/env python3
# Deterministic exact radius-5 quotient witness: lexicographically minimal 4-term R3b_123 cycle.
from __future__ import annotations
import json, sympy as sp, hashlib
from pathlib import Path
from fractions import Fraction
B=Path(__file__).resolve().parent
g3=json.load(open(B/'gamma3_level2_relative_3cells_v36.json'));cells=g3['cells'];gens=g3['generator_order'];gi={g:i for i,g in enumerate(gens)}
r3=next(i for i,c in enumerate(cells) if c['name']=='R3b_123')
cert4=json.load(open(B/'gamma3_radius4_mod4_dual_certificate_v50.json'));bits=tuple(cert4['critical_character']['bits']);chi_g={g:(-1 if bits[gi[g]] else 1) for g in gens}
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
def mod4(A):return tuple(tuple(x%4 for x in row) for row in A)
def mm4(A,C):return tuple(tuple((A[i][0]*C[0][j]+A[i][1]*C[1][j]+A[i][2]*C[2][j])%4 for j in range(3)) for i in range(3))
I=mt(sp.eye(3));I4=mod4(I);MAT={g:mt(M) for g,M in STD.items()};MATI={g:mi(MAT[g]) for g in MAT};MAT4={g:mod4(MAT[g]) for g in MAT};MATI4={g:mod4(MATI[g]) for g in MAT}
letters=[]
for g in gens:
 letters.append((g,1,MAT[g]))
 if MATI[g]!=MAT[g]:letters.append((g,-1,MATI[g]))
ball={I:(0,1,())};front=[I]
for d in range(1,6):
 nf=[]
 for a in front:
  _,cha,wa=ball[a]
  for g,s,M in letters:
   h=mm(a,M);ch=cha*chi_g[g]
   if h not in ball:ball[h]=(d,ch,wa+((g,s),));nf.append(h)
   else:assert ball[h][1]==ch
 front=nf
b5={}
for A,(d,ch,w) in ball.items():
 k=mod4(A)
 if k not in b5 or (d,repr(w),repr(A))<(b5[k][0],repr(b5[k][2]),repr(b5[k][3])):b5[k]=(d,ch,w,A)
 else:assert b5[k][1]==ch
assert len(b5)==382
def fox4(w):
 out={x:{} for x in gens};pref=I4
 for t in w:
  inv=t.endswith('^-1');b=t[:-3] if inv else t
  if inv:q=MATI4[b];key=mm4(pref,q);out[b][key]=out[b].get(key,0)-1;pref=mm4(pref,q)
  else:out[b][pref]=out[b].get(pref,0)+1;pref=mm4(pref,MAT4[b])
 assert pref==I4;return out
F4=fox4(cells[r3]['word'])
def col(h):
 d={}
 for x in gens:
  for q,c in F4[x].items():
   key=(x,mm4(h,q));d[key]=d.get(key,0)+c
 return {k:v for k,v in d.items() if v}
def add(a,b):
 d=dict(a)
 for k,v in b.items():
  d[k]=d.get(k,0)+v
  if d[k]==0:del d[k]
 return d
def key(d):return tuple(sorted(((repr(k),v) for k,v in d.items())))
def negkey(k):return tuple((r,-v) for r,v in k)
items=sorted([(h,meta) for h,meta in b5.items() if meta[1]==-1],key=lambda z:repr(z[0]))
cols=[col(h) for h,m in items]
# Store all pairs by exact sparse sum; enumerate lexicographically and select minimal disjoint quadruple.
pairs={}
for i in range(len(items)):
 for j in range(i+1,len(items)):
  K=key(add(cols[i],cols[j]));pairs.setdefault(K,[]).append((i,j))
best=None
for K,ps in pairs.items():
 qs=pairs.get(negkey(K),[])
 if not qs:continue
 for a in ps:
  for b in qs:
   inds=tuple(sorted(a+b))
   if len(set(inds))<4:continue
   if best is None or inds<best:best=inds
if best is None:raise SystemExit('no 4-term exact quotient cycle')
assert not add(add(cols[best[0]],cols[best[1]]),add(cols[best[2]],cols[best[3]]))
# integral Fox boundary for canonical representatives
def fox_int(word):
 out={x:{} for x in gens};pref=I
 for tok in word:
  inv=tok.endswith('^-1');b=tok[:-3] if inv else tok
  if inv:q=MATI[b];kk=mm(pref,q);out[b][kk]=out[b].get(kk,0)-1;pref=mm(pref,q)
  else:out[b][pref]=out[b].get(pref,0)+1;pref=mm(pref,MAT[b])
 assert pref==I;return out
FI=fox_int(cells[r3]['word']);actual={};terms=[]
for idx in best:
 h,(d,ch,w,A)=items[idx]
 terms.append({'relator_index':r3,'relator_name':'R3b_123','group_mod4':[list(row) for row in h],'depth':d,'character':ch,'word_representative':[[a,s] for a,s in w],'coefficient_numerator':-1,'coefficient_denominator':4,'integral_matrix_representative':[list(row) for row in A]})
 # multiplied-by-4 chain has coefficient -1
 for x in gens:
  for q,c in FI[x].items():
   kk=(x,mm(A,q));actual[kk]=actual.get(kk,0)-c
actual={k:v for k,v in actual.items() if v}
resid=[{'generator':x,'group_matrix':[list(row) for row in A],'coefficient':c} for (x,A),c in sorted(actual.items(),key=lambda z:(z[0][0],repr(z[0][1])))]
out={'version':'v51','status':'PROVEN-EXACT-DETERMINISTIC-FINITE-QUOTIENT-CYCLE','support_radius':5,'selection_rule':'lexicographically minimal four distinct chi=-1 R3b_123 translation classes whose exact mod4 Fox columns sum to zero','candidate_chi_minus_classes':len(items),'quotient_cycle_term_count':4,'common_denominator':4,'target_projection':'1','terms':terms,'chosen_representatives_actual_group_ring_boundary_nonzero_terms':len(actual),'chosen_representatives_form_integral_lift':len(actual)==0,'actual_boundary_residual':resid,'decision':'FINITE_QUOTIENT_OBSTRUCTION_DISAPPEARS_AT_RADIUS5; DETERMINISTIC QUOTIENT WITNESS FIXED; CANONICAL-SECTION INTEGRAL LIFT '+('PROVEN' if not actual else 'REFUTED-TYPED'),'guard':'The exact quotient witness is not by itself an integral group-ring lift. Nonzero actual boundary for its canonical representatives does not exclude kernel corrections in other mod4 fibers.'}
p=B/'gamma3_radius5_quotient_cycle_v51.json';p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:out[k] for k in ['status','selection_rule','candidate_chi_minus_classes','quotient_cycle_term_count','target_projection','chosen_representatives_actual_group_ring_boundary_nonzero_terms','chosen_representatives_form_integral_lift','decision']},indent=2));print('terms',[(t['depth'],t['word_representative']) for t in terms]);print('sha256',hashlib.sha256(p.read_bytes()).hexdigest())
