#!/usr/bin/env python3
import json,sys
from pathlib import Path
import sympy as sp
B=Path(__file__).resolve().parent
cert=json.loads((B/'gamma3_radius3_mod4_dual_certificate_v48.json').read_text())
g3=json.loads((B/'gamma3_level2_relative_3cells_v36.json').read_text()); cells=g3['cells']; gens=g3['generator_order']; gi={g:i for i,g in enumerate(gens)}
r3=next(i for i,c in enumerate(cells) if c['name']=='R3b_123')
STD={}
for i in range(1,4):
 F=sp.eye(3);F[i-1,i-1]=-1;STD[f'F{i}']=F
for i in range(1,4):
 for j in range(1,4):
  if i!=j:
   E=sp.eye(3);E[i-1,j-1]=2;STD[f'E{i}{j}']=E
def mt(M):return tuple(tuple(int(M[i,j]) for j in range(3)) for i in range(3))
def mm(A,C):return tuple(tuple(sum(A[i][k]*C[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def mi(A):return mt(sp.Matrix(A).inv())
def mod4(A):return tuple(tuple(x%4 for x in row) for row in A)
def mm4(A,C):return tuple(tuple(sum(A[i][k]*C[k][j] for k in range(3))%4 for j in range(3)) for i in range(3))
I=mt(sp.eye(3)); I4=mod4(I); MAT={g:mt(M) for g,M in STD.items()}; MATI={g:mi(MAT[g]) for g in MAT}; MAT4={g:mod4(MAT[g]) for g in MAT}; MATI4={g:mod4(MATI[g]) for g in MAT}
# full quotient group
letters4=[]
for g in gens:
 letters4.append(MAT4[g])
 if MATI4[g]!=MAT4[g]: letters4.append(MATI4[g])
allG={I4};front=[I4]
while front:
 nf=[]
 for a in front:
  for M in letters4:
   h=mm4(a,M)
   if h not in allG:allG.add(h);nf.append(h)
 front=nf
assert len(allG)==512
# exact B3 and collapsed image; character F_i -> -1, E_ij -> +1
bits=(0,0,0,0,0,0,1,1,1); chi_g={g:(-1 if bits[gi[g]] else 1) for g in gens}
letters=[]
for g in gens:
 letters.append((g,1,MAT[g]))
 if MATI[g]!=MAT[g]: letters.append((g,-1,MATI[g]))
ball={I:(0,1,())};front=[I]
prof={0:1}
for d in (1,2,3):
 nf=[]
 for a in front:
  for g,s,M in letters:
   h=mm(a,M); ch=ball[a][1]*chi_g[g]
   if h not in ball: ball[h]=(d,ch,ball[a][2]+((g,s),)); nf.append(h)
 front=nf;prof[d]=len(nf)
assert len(ball)==1376 and prof=={0:1,1:15,2:147,3:1213}
b4={}
for A,(d,ch,w) in ball.items():
 k=mod4(A)
 if k not in b4 or d<b4[k][0]: b4[k]=(d,ch,w)
 else: assert b4[k][1]==ch
assert len(b4)==130
# Fox boundary in quotient
def fox4(w):
 out={x:{} for x in gens};pref=I4
 for t in w:
  inv=t.endswith('^-1');b=t[:-3] if inv else t
  if inv:
   q=MATI4[b]; key=mm4(pref,q); out[b][key]=out[b].get(key,0)-1; pref=mm4(pref,q)
  else:
   out[b][pref]=out[b].get(pref,0)+1; pref=mm4(pref,MAT4[b])
 assert pref==I4
 return out
FOX=[fox4(c['word']) for c in cells]
cols=[];rowset=set();colmaps=[]
for k in range(43):
 for h,(d,ch,w) in b4.items():
  dd={}
  for x in gens:
   for q,c in FOX[k][x].items():
    key=(x,mm4(h,q));dd[key]=dd.get(key,0)+c
  dd={a:c for a,c in dd.items() if c};rowset.update(dd);cols.append((k,h,d,ch,w));colmaps.append(dd)
rows=sorted(rowset,key=repr);ri={r:i for i,r in enumerate(rows)}
assert [len(rows),len(cols)]==cert['matrix_shape']==[3624,5590]
# Reconstruct the dual functional by semantic row key, not stored row index.
ykey={}
for term in cert['dual_terms']:
    key=(term['generator'],tuple(tuple(int(x) for x in row) for row in term['group_mod4']))
    ykey[key]=int(term['coefficient'])
den=int(cert['denominator'])
assert len(ykey)==127 and den==1
# Verify y^T D = l exactly, where l is character amplitude on R3b_123.
for j,(k,h,d,ch,w) in enumerate(cols):
    got=sum(ykey.get(key,0)*c for key,c in colmaps[j].items())
    want=den*(ch if k==r3 else 0)
    if got!=want:
        print('FAIL column',j,got,want,file=sys.stderr);sys.exit(1)
print(json.dumps({'status':'PASS','group_mod4':len(allG),'B3':len(ball),'B3_image_mod4':len(b4),'shape':[len(rows),len(cols)],'dual_terms':len(ykey),'identity':'EXACT'}))
