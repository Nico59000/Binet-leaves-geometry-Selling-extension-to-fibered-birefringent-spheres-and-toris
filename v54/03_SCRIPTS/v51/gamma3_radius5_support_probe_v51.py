#!/usr/bin/env python3
import json, sympy as sp, hashlib
from pathlib import Path
B=Path('/mnt/data/v51_work')
g3=json.load(open(B/'gamma3_level2_relative_3cells_v36.json'))
cells=g3['cells']; gens=g3['generator_order']; gi={g:i for i,g in enumerate(gens)}
STD={}
for i in range(1,4):
 F=sp.eye(3); F[i-1,i-1]=-1; STD[f'F{i}']=F
for i in range(1,4):
 for j in range(1,4):
  if i!=j:
   E=sp.eye(3); E[i-1,j-1]=2; STD[f'E{i}{j}']=E
def mt(M):return tuple(tuple(int(M[i,j]) for j in range(3)) for i in range(3))
def mm(A,C):return tuple(tuple(A[i][0]*C[0][j]+A[i][1]*C[1][j]+A[i][2]*C[2][j] for j in range(3)) for i in range(3))
def mi(A):return mt(sp.Matrix(A).inv())
def mod4(A):return tuple(tuple(x%4 for x in row) for row in A)
def mm4(A,C):return tuple(tuple((A[i][0]*C[0][j]+A[i][1]*C[1][j]+A[i][2]*C[2][j])%4 for j in range(3)) for i in range(3))
I=mt(sp.eye(3)); MAT={g:mt(M) for g,M in STD.items()}; MATI={g:mi(MAT[g]) for g in MAT}; MAT4={g:mod4(MAT[g]) for g in MAT}; MATI4={g:mod4(MATI[g]) for g in MAT}
cert=json.load(open(B/'gamma3_radius4_mod4_dual_certificate_v50.json'))
bits=tuple(cert['critical_character']['bits']); chi_g={g:(-1 if bits[gi[g]] else 1) for g in gens}
letters=[]
for g in gens:
 letters.append((g,1,MAT[g]))
 if MATI[g]!=MAT[g]:letters.append((g,-1,MATI[g]))
ball={I:(0,1,())};front=[I]; prof={0:1}; imgprof={0:1}; seenmod={mod4(I)}
for d in range(1,6):
 nf=[]
 for a in front:
  _,cha,wa=ball[a]
  for g,s,M in letters:
   h=mm(a,M); ch=cha*chi_g[g]
   if h not in ball:
    ball[h]=(d,ch,wa+((g,s),));nf.append(h)
   else:
    assert ball[h][1]==ch
 front=nf;prof[d]=len(nf)
 newmods=0
 for a in nf:
  k=mod4(a)
  if k not in seenmod:seenmod.add(k);newmods+=1
 imgprof[d]=newmods
 print('depth',d,'new',len(nf),'ball',len(ball),'newmod',newmods,'modtotal',len(seenmod),flush=True)
# select minimal reps mod4 and check character consistency on all B5 collisions
b5={}; collisions=0
for A,(d,ch,w) in ball.items():
 k=mod4(A)
 if k not in b5 or d<b5[k][0]: b5[k]=(d,ch,w)
 else:
  if b5[k][1]!=ch:
   collisions+=1
out={'version':'v51','status':'PASS-EXACT-RADIUS5-SUPPORT-PROBE' if collisions==0 else 'FAIL-CHARACTER-COLLISION','ball_profile':prof,'ball_size':len(ball),'mod4_new_profile':imgprof,'mod4_image_size':len(b5),'character_collision_count':collisions,'source_radius4_image_size':cert['radius4_image_mod4_size'],'full_mod4_quotient_size':512,'decision':'RADIUS5_FINITE_QUOTIENT_SUPPORT_READY' if collisions==0 else 'STOP'}
p=B/'gamma3_radius5_support_gate_v51.json';p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2));print('sha256',hashlib.sha256(p.read_bytes()).hexdigest())
