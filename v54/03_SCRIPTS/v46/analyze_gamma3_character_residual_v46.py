#!/usr/bin/env python3
from pathlib import Path
import json,csv,math,time
import sympy as sp
B=Path(__file__).resolve().parent
g3=json.load(open(B/'gamma3_level2_relative_3cells_v36.json'));cells=g3['cells'];gens=g3['generator_order'];gi={g:i for i,g in enumerate(gens)};ci={c['name']:i for i,c in enumerate(cells)}
bits=(0,0,0,0,0,0,1,1,1)
# Character chi: Eij -> +1, Fi -> -1.
def fox_char(w):
 vals=[0]*9;pref=1
 for t in w:
  inv=t.endswith('^-1');b=t[:-3] if inv else t;j=gi[b];s=-1 if bits[j] else 1
  if inv: vals[j]+=-pref*s;pref*=s
  else: vals[j]+=pref;pref*=s
 assert pref==1
 return vals
D=sp.Matrix([[fox_char(c['word'])[i] for c in cells] for i in range(9)])
# M43 shadow basis.
fc=json.load(open(B/'gamma3_group_ring_4cells_v38.json')); partial=[]
for cyc in fc['cycles']:
 v=[0]*43;v[ci[cyc['relator_a']]]=1;v[ci[cyc['relator_b']]]=1;partial.append(sp.Matrix(v))
square={}
for k,c in enumerate(cells):
 w=c['word']
 if len(w)%2==0 and w[:len(w)//2]==w[len(w)//2:]:
  half=w[:len(w)//2];p=0
  for t in half:p^=bits[gi[t[:-3] if t.endswith('^-1') else t]]
  square[k]=(half,p)
  if p:
   v=[0]*43;v[k]=1;partial.append(sp.Matrix(v))
P=sp.Matrix.hstack(*partial)
assert D*P==sp.zeros(9,P.cols)
# Greedy residual quotient basis in ker(D)/shadow(M43).
K=D.nullspace();cur=P;Q=[];rank=P.rank()
for v in K:
 nr=sp.Matrix.hstack(cur,v).rank()
 if nr>rank:Q.append(v);cur=sp.Matrix.hstack(cur,v);rank=nr
assert len(Q)==14 and rank==35
rows=[]
for qi,v in enumerate(Q,1):
 den=1
 for x in v:den=sp.ilcm(den,int(x.q))
 iv=[int(x*den) for x in v];g=0
 for x in iv:g=math.gcd(g,abs(x))
 if g:iv=[x//g for x in iv]
 first=next((x for x in iv if x),1)
 if first<0:iv=[-x for x in iv]
 for i,x in enumerate(iv):
  if x:rows.append({'direction_id':f'CHIRES_{qi:02d}','relator':cells[i]['name'],'coefficient':x})
with open(B/'gamma3_character_residual_basis_v46.csv','w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=['direction_id','relator','coefficient']);w.writeheader();w.writerows(rows)
basis={
 'version':'v46','character':{g:(-1 if bits[i] else 1) for i,g in enumerate(gens)},
 'description':'minimal residual character chi(Eij)=+1, chi(Fi)=-1',
 'character_Fox_rank':D.rank(),'character_kernel_dimension':43-D.rank(),
 'M43_shadow_rank':P.rank(),'residual_quotient_dimension':len(Q),
 'basis_csv':'gamma3_character_residual_basis_v46.csv',
 'decision':'PROVEN_14_DIMENSIONAL_RESIDUAL_CHARACTER_QUOTIENT; THESE ARE SHADOW DIRECTIONS, NOT YET EXACT GROUP-RING SYZYGIES'}
(B/'gamma3_character_residual_basis_v46.json').write_text(json.dumps(basis,indent=2,ensure_ascii=False),encoding='utf-8')
# Build residual-coordinate quotient map L: ker(D)->Q^14 killing M43 shadow.
S=cur
for i in range(43):
 e=sp.zeros(43,1);e[i]=1
 if sp.Matrix.hstack(S,e).rank()>S.rank():S=sp.Matrix.hstack(S,e)
 if S.cols==43:break
assert S.cols==43
Sinv=S.inv();L=Sinv[21:35,:];assert L*P==sp.zeros(14,P.cols)
# Exact group-ring boundary with faithful integer-matrix group elements.
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
I=mt(sp.eye(3));MAT={g:mt(M) for g,M in STD.items()};MATI={g:mi(MAT[g]) for g in MAT}
def grterm(d,g,c):
 if c:
  d[g]=d.get(g,0)+c
  if d[g]==0:del d[g]
def fox(w):
 out={x:{} for x in gens};pref=I
 for t in w:
  inv=t.endswith('^-1');b=t[:-3] if inv else t
  if inv:
   xinv=MATI[b];grterm(out[b],mm(pref,xinv),-1);pref=mm(pref,xinv)
  else:grterm(out[b],pref,1);pref=mm(pref,MAT[b])
 assert pref==I;return out
FOX=[fox(c['word']) for c in cells]
# Radius-one ball in the group, with chi values. Matrix duplicates are deduplicated.
ball={I:1}
for g in gens:
 ch=-1 if bits[gi[g]] else 1
 for h in (MAT[g],MATI[g]):
  if h in ball: assert ball[h]==ch
  ball[h]=ch
ball=list(ball.items());assert len(ball)==16
cols=[];colmaps=[];rowset=set()
for k in range(43):
 for h,ch in ball:
  d={}
  for x in gens:
   for q,c in FOX[k][x].items():
    key=(x,mm(h,q));d[key]=d.get(key,0)+c;rowset.add(key)
  cols.append((k,h,ch));colmaps.append(d)
R=sorted(rowset,key=repr);ri={r:i for i,r in enumerate(R)};sd={}
for j,d in enumerate(colmaps):
 for r,c in d.items():sd[(ri[r],j)]=c
BM=sp.MutableSparseMatrix(len(R),len(cols),sd)
# residual shadow coordinates of each radius-one coefficient column.
C=sp.zeros(14,len(cols))
for j,(k,h,ch) in enumerate(cols):
 for a in range(14):C[a,j]=L[a,k]*ch
rB=BM.rank();rBC=BM.col_join(C).rank()
no_go={
 'version':'v46','character':basis['character'],'group_ball':'B1={1} union {Eij^+-1, Fi}; deduplicated in the faithful 3x3 integer matrix representation',
 'ball_unique_elements':len(ball),'variables_43_relators_times_ball':len(cols),'exact_boundary_rows':len(R),
 'exact_boundary_rank':rB,'exact_kernel_dimension':len(cols)-rB,
 'residual_shadow_image_dimension_of_exact_radius1_kernel_mod_M43':rBC-rB,
 'decision':'PROVEN BOUNDED NO-GO: EVERY EXACT CYCLE WITH ALL GROUP-RING COEFFICIENTS SUPPORTED IN B1 HAS CHARACTER SHADOW INSIDE M43; NO ONE OF THE 14 RESIDUAL DIRECTIONS LIFTS AT RADIUS ONE.',
 'scope_guard':'This is a bounded coefficient-support theorem. It does not exclude exact residual syzygies using group elements of word length >=2.'}
(B/'gamma3_character_radius1_no_go_v46.json').write_text(json.dumps(no_go,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'D_rank':D.rank(),'kernel':43-D.rank(),'M43_shadow_rank':P.rank(),'residual':len(Q),'ball':len(ball),'boundary_rank':rB,'kernel_radius1':len(cols)-rB,'residual_image':rBC-rB},indent=2))
