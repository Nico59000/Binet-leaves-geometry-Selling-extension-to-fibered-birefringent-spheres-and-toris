#!/usr/bin/env python3
from __future__ import annotations
import json,sys
from pathlib import Path
import sympy as sp
B=Path(__file__).resolve().parent
cert=json.loads((B/'gamma3_radius4_mod4_dual_certificate_v50.json').read_text())
g3=json.loads((B/'gamma3_level2_relative_3cells_v36.json').read_text())
cells=g3['cells']; gens=g3['generator_order']; gi={g:i for i,g in enumerate(gens)}
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
# Full quotient size, independently of the bounded support image.
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
# Exact radius-four ball and the critical character.
bits=tuple(cert['critical_character']['bits']); assert cert['critical_character']['generator_order']==gens
chi_g={g:(-1 if bits[gi[g]] else 1) for g in gens}
letters=[]
for g in gens:
 letters.append((g,1,MAT[g]))
 if MATI[g]!=MAT[g]: letters.append((g,-1,MATI[g]))
ball={I:(0,1,())};front=[I]; prof={0:1}
for d in (1,2,3,4):
 nf=[]
 for a in front:
  for g,s,M in letters:
   h=mm(a,M); ch=ball[a][1]*chi_g[g]
   if h not in ball:
    ball[h]=(d,ch,ball[a][2]+((g,s),)); nf.append(h)
   else:
    assert ball[h][1]==ch
 front=nf; prof[d]=len(nf)
assert len(ball)==10640 and prof=={0:1,1:15,2:147,3:1213,4:9264}
b4={}
for A,(d,ch,w) in ball.items():
 k=mod4(A)
 if k not in b4 or d<b4[k][0]: b4[k]=(d,ch,w)
 else: assert b4[k][1]==ch
assert len(b4)==256
# Fox boundary in quotient.
def fox4(w):
 out={x:{} for x in gens};pref=I4
 for t in w:
  inv=t.endswith('^-1'); b=t[:-3] if inv else t
  if inv:
   q=MATI4[b]; key=mm4(pref,q); out[b][key]=out[b].get(key,0)-1; pref=mm4(pref,q)
  else:
   out[b][pref]=out[b].get(pref,0)+1; pref=mm4(pref,MAT4[b])
 assert pref==I4
 return out
FOX=[fox4(c['word']) for c in cells]
items=sorted(b4.items(),key=lambda z:repr(z[0]))
cols=[];rowset=set();colmaps=[]
for k in range(43):
 for h,(d,ch,w) in items:
  dd={}
  for x in gens:
   for q,c in FOX[k][x].items():
    key=(x,mm4(h,q)); dd[key]=dd.get(key,0)+c
  dd={a:c for a,c in dd.items() if c}; rowset.update(dd); cols.append((k,h,d,ch,w)); colmaps.append(dd)
rows=sorted(rowset,key=repr)
assert [len(rows),len(cols)]==cert['matrix_shape']==[4424,11008]
assert sum(len(x) for x in colmaps)==cert['matrix_nnz']==50688
# Semantic certificate, common denominator 2.
ykey={}
for term in cert['dual_terms']:
 key=(term['generator'],tuple(tuple(int(x) for x in row) for row in term['group_mod4']))
 assert key not in ykey
 ykey[key]=int(term['coefficient_numerator_over_2'])
assert len(ykey)==cert['dual_term_count']==319
assert int(cert['dual_denominator'])==2
assert set(ykey)<=set(rows)
# Exact integral identity: (2y)^T D = 2l.
for j,(k,h,d,ch,w) in enumerate(cols):
 got=sum(ykey.get(key,0)*c for key,c in colmaps[j].items())
 want=2*(ch if k==r3 else 0)
 if got!=want:
  print(json.dumps({'status':'FAIL','column':j,'got':got,'want':want}),file=sys.stderr);sys.exit(1)
print(json.dumps({
 'status':'PASS','full_mod4_quotient_size':len(allG),'B4':len(ball),
 'B4_image_mod4':len(b4),'shape':[len(rows),len(cols)],'nnz':sum(len(x) for x in colmaps),
 'dual_terms':len(ykey),'dual_denominator':2,'identity':'EXACT_INTEGER_COMMON_DENOMINATOR_2',
 'decision':'CHIRES14_NO_LIFT_RADIUS_LE_4_MIN_IF_ANY_GE_5'
},indent=2))
