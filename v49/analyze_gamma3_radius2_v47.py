#!/usr/bin/env python3
from pathlib import Path
import json,csv,math,collections
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ
B=Path(__file__).resolve().parent
g3=json.load(open(B/'gamma3_level2_relative_3cells_v36.json',encoding='utf-8'))
cells=g3['cells'];gens=g3['generator_order'];gi={g:i for i,g in enumerate(gens)};ci={c['name']:i for i,c in enumerate(cells)}
bits=(0,0,0,0,0,0,1,1,1)
# Character chi(Eij)=+1, chi(Fi)=-1 and residual quotient from v46 construction.
def fox_char(w):
 vals=[0]*9;pref=1
 for t in w:
  inv=t.endswith('^-1');b=t[:-3] if inv else t;j=gi[b];s=-1 if bits[j] else 1
  if inv: vals[j]+=-pref*s;pref*=s
  else: vals[j]+=pref;pref*=s
 assert pref==1
 return vals
D=sp.Matrix([[fox_char(c['word'])[i] for c in cells] for i in range(9)])
fc=json.load(open(B/'gamma3_group_ring_4cells_v38.json',encoding='utf-8'));partial=[]
for cyc in fc['cycles']:
 v=[0]*43;v[ci[cyc['relator_a']]]=1;v[ci[cyc['relator_b']]]=1;partial.append(sp.Matrix(v))
for k,c in enumerate(cells):
 w=c['word']
 if len(w)%2==0 and w[:len(w)//2]==w[len(w)//2:]:
  p=0
  for t in w[:len(w)//2]:p^=bits[gi[t[:-3] if t.endswith('^-1') else t]]
  if p:
   v=[0]*43;v[k]=1;partial.append(sp.Matrix(v))
P=sp.Matrix.hstack(*partial)
assert D*P==sp.zeros(9,P.cols) and P.rank()==21
K=D.nullspace();cur=P;Q=[];rk=P.rank()
for v in K:
 nr=sp.Matrix.hstack(cur,v).rank()
 if nr>rk:Q.append(v);cur=sp.Matrix.hstack(cur,v);rk=nr
assert len(Q)==14 and rk==35
S=cur
for i in range(43):
 e=sp.zeros(43,1);e[i]=1
 if sp.Matrix.hstack(S,e).rank()>S.rank():S=sp.Matrix.hstack(S,e)
 if S.cols==43:break
L=S.inv()[21:35,:]
assert L*P==sp.zeros(14,P.cols)
# Faithful integer-matrix carrier for Gamma_3(2).
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
 assert pref==I
 return out
FOX=[fox(c['word']) for c in cells]
# Exact radius <=2 group ball with shortest word representatives.
letters=[]
for g in gens:
 letters.append((g,1,MAT[g]))
 if MATI[g]!=MAT[g]:letters.append((g,-1,MATI[g]))
chi_g={g:(-1 if bits[gi[g]] else 1) for g in gens}
ball={I:{'depth':0,'chi':1,'word':()}};front=[I]
for dep in (1,2):
 new=[]
 for a in front:
  for g,sg,M in letters:
   h=mm(a,M);ch=ball[a]['chi']*chi_g[g]
   if h not in ball:
    ball[h]={'depth':dep,'chi':ch,'word':ball[a]['word']+((g,sg),)};new.append(h)
   else:assert ball[h]['chi']==ch
 front=new
assert len(ball)==163
ballitems=sorted(ball.items(),key=lambda z:repr(z[0]))
# Exact Fox boundary matrix for all 43 relator coefficients supported in B2.
cols=[];colmaps=[];rowset=set()
for k in range(43):
 for h,meta in ballitems:
  dd={}
  for x in gens:
   for q,c in FOX[k][x].items():
    key=(x,mm(h,q));dd[key]=dd.get(key,0)+c;rowset.add(key)
  cols.append((k,h,meta['chi']));colmaps.append(dd)
rows=sorted(rowset,key=repr);ri={r:i for i,r in enumerate(rows)};dok={}
for j,dd in enumerate(colmaps):
 for r,c in dd.items():
  if c:dok[(ri[r],j)]=QQ.convert(c)
DM=DomainMatrix.from_dok(dok,(len(rows),len(cols)),QQ)
rB=DM.rank();N=DM.nullspace();assert N.shape[1]==len(cols)
# Residual shadow image of the exact B2 kernel.
nd=N.to_dok();nrows=[{} for _ in range(N.shape[0])]
for (i,j),v in nd.items(): nrows[i][j]=sp.Rational(int(v.numerator),int(v.denominator))
Y=sp.zeros(14,N.shape[0])
for i,row in enumerate(nrows):
 for j,x in row.items():
  k,h,ch=cols[j]
  for a in range(14):
   if L[a,k]:Y[a,i]+=L[a,k]*ch*x
rY=Y.rank();assert rY==13
left=Y.T.nullspace();assert len(left)==1
lv=left[0];den=1
for x in lv:den=sp.ilcm(den,int(x.q))
rel=[int(x*den) for x in lv];gg=0
for x in rel:gg=math.gcd(gg,abs(x))
if gg:rel=[x//gg for x in rel]
if next(x for x in rel if x)<0:rel=[-x for x in rel]
assert rel==[0]*13+[1]
# Construct pure exact lifts for CHIRES_01..CHIRES_13.
_,piv=Y.rref();piv=list(piv);assert len(piv)==13
A=Y[:13,piv];assert A.det()!=0
outrows=[];summaries=[]
for target in range(13):
 e=sp.zeros(13,1);e[target]=1;alpha=A.inv()*e
 vec={}
 for coeff,ni in zip(alpha,piv):
  if coeff==0:continue
  for j,x in nrows[ni].items():
   vec[j]=vec.get(j,sp.Rational(0))+coeff*x
   if vec[j]==0:del vec[j]
 d=1
 for x in vec.values():d=sp.ilcm(d,int(x.q))
 iv={j:int(x*d) for j,x in vec.items()};g=0
 for x in iv.values():g=math.gcd(g,abs(x))
 if g:iv={j:x//g for j,x in iv.items()}
 # normalize sign to make residual amplitude positive.
 yy=sp.zeros(14,1)
 for j,x in iv.items():
  k,h,ch=cols[j]
  for a in range(14):yy[a]+=L[a,k]*ch*x
 nz=[(i,z) for i,z in enumerate(yy) if z]
 assert len(nz)==1 and nz[0][0]==target
 if nz[0][1]<0:
  iv={j:-x for j,x in iv.items()};yy=-yy;nz=[(target,-nz[0][1])]
 amp=int(nz[0][1]) if nz[0][1].q==1 else str(nz[0][1])
 assert amp==2
 lift=f'R2PURE_{target+1:02d}'
 maxdep=max(ball[cols[j][1]]['depth'] for j in iv)
 summaries.append({'lift_id':lift,'target':f'CHIRES_{target+1:02d}','term_count':len(iv),'max_group_word_radius':maxdep,'residual_amplitude':amp,'exact_Fox_boundary':'ZERO'})
 for j,x in sorted(iv.items()):
  k,h,ch=cols[j];w=ball[h]['word'];ws='1' if not w else '*'.join(g if sg==1 else g+'^-1' for g,sg in w)
  outrows.append({'lift_id':lift,'target':f'CHIRES_{target+1:02d}','relator':cells[k]['name'],'group_word':ws,'group_radius':ball[h]['depth'],'coefficient':x})
with (B/'gamma3_radius2_pure_lifts_v47.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=['lift_id','target','relator','group_word','group_radius','coefficient']);w.writeheader();w.writerows(outrows)
cert={
 'version':'v47','character':{g:(-1 if bits[i] else 1) for i,g in enumerate(gens)},
 'radius2_ball':{'unique_group_elements':len(ball),'depth_profile':dict(collections.Counter(m['depth'] for m in ball.values()))},
 'exact_boundary':{'rows':len(rows),'columns':len(cols),'rank_Q':rB,'kernel_dimension':len(cols)-rB,'kernel_basis_rows':N.shape[0]},
 'v46_residual_quotient_dimension':14,'residual_image_dimension_of_exact_radius2_kernel_mod_M43':rY,
 'residual_annihilator_relation':{'coordinates':[f'CHIRES_{i:02d}' for i in range(1,15)],'primitive_relation':rel,'interpretation':'CHIRES_14 coordinate is zero on every exact radius<=2 cycle modulo M43 shadow'},
 'pure_exact_lifts':summaries,
 'lifted_directions':[f'CHIRES_{i:02d}' for i in range(1,14)],
 'unlifted_direction_radius2':'CHIRES_14','unlifted_direction_support':'R3b_123',
 'minimal_radius_conclusion':'Because v46 proved zero residual image at radius 1, the thirteen lifted directions have exact minimal coefficient-support radius 2. CHIRES_14 has no exact lift at radius <=2 and therefore requires radius >=3 if it lifts at all.',
 'decision':'PROVEN: EXACT RADIUS-2 GROUP-RING CYCLES LIFT 13 OF THE 14 CHARACTER-RESIDUAL DIRECTIONS; CHIRES_14 REMAINS A BOUNDED NO-GO THROUGH RADIUS 2.',
 'scope_guard':'Character-shadow independence proves the thirteen cycles are new modulo M43 after chi-specialization. It does not by itself prove a free rank-13 R-submodule or completeness of the full identity module.'}
(B/'gamma3_character_radius2_lifts_v47.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False),encoding='utf-8')
partial={
 'version':'v47','base_module':'M43 = R^6 direct-sum 18 order-two cyclic sectors (certified partial resolution)',
 'new_exact_cycles':[x['lift_id'] for x in summaries],
 'new_cycle_count':13,'shadow_independence_mod_M43':13,
 'character_quotient_before':14,'character_quotient_after_adjoining_new_cycles':1,
 'remaining_character_generator':'CHIRES_14 = R3b_123',
 'full_R_module_intersections_among_new_generators':'NT','complete_resolution':'NT',
 'status':'PROVEN-EXTENSION-AT-CHARACTER-SHADOW / FULL-R-MODULE-STRUCTURE-NT'}
(B/'gamma3_partial_resolution_v47.json').write_text(json.dumps(partial,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'ball':len(ball),'rows':len(rows),'cols':len(cols),'rank':rB,'kernel':len(cols)-rB,'residual_image':rY,'relation':rel,'pure_lifts':len(summaries)},indent=2))
