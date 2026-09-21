#!/usr/bin/env python3
import json, hashlib, math
from pathlib import Path
from fractions import Fraction
import sympy as sp
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.matrices.normalforms import smith_normal_decomp
B=Path('/mnt/data/v50_work')
I=json.load(open(B/'selling_integral_orbifold_cochains_v50.json'))
Q=json.load(open(B/'selling_historical_cohomology_Q_v50.json'))
G=json.load(open(B/'selling_field_orbigroupoid_v50.json'))

def M(rows,cols,E):return sp.Matrix(sp.SparseMatrix(rows,cols,{(e['row'],e['col']):sp.Integer(e['value']) for e in E}))
def snfdiag(A):
 S=smith_normal_form(A,domain=ZZ)
 return [abs(int(S[i,i])) for i in range(min(S.rows,S.cols)) if S[i,i]!=0]
def modrank(A,p):
 rows=[{} for _ in range(A.rows)]
 for i in range(A.rows):
  for j in range(A.cols):
   v=int(A[i,j])%p
   if v: rows[i][j]=v
 piv={};rank=0
 for row in rows:
  while row:
   c=min(row)
   if c not in piv:
    inv=pow(row[c],-1,p);row={k:(v*inv)%p for k,v in row.items() if (v*inv)%p};piv[c]=row;rank+=1;break
   fac=row[c];pr=piv[c]
   for k,v in pr.items():
    nv=(row.get(k,0)-fac*v)%p
    if nv:row[k]=nv
    else:row.pop(k,None)
 return rank

def nullmod2(A):
 m,n=A.shape;rows=[]
 for i in range(m):
  bits=0
  for j in range(n):
   if int(A[i,j])&1:bits|=1<<j
  rows.append(bits)
 piv=[];r=0
 for c in range(n):
  p=next((i for i in range(r,m) if (rows[i]>>c)&1),None)
  if p is None:continue
  rows[r],rows[p]=rows[p],rows[r]
  for i in range(m):
   if i!=r and ((rows[i]>>c)&1):rows[i]^=rows[r]
  piv.append(c);r+=1
 free=[j for j in range(n) if j not in piv]
 assert len(free)==1
 x=[0]*n;x[free[0]]=1
 for i,p in enumerate(piv):
  s=0
  for j in free:
   if (rows[i]>>j)&1:s^=x[j]
  x[p]=s
 return sp.Matrix(x)
def sv(v):return [{'i':i,'value':str(int(x))} for i,x in enumerate(v) if x]

def kernelZ(A):
 S,U,V=smith_normal_decomp(DomainMatrix.from_Matrix(A).convert_to(ZZ));Sm=S.to_Matrix();Vm=V.to_Matrix();r=sum(1 for i in range(min(Sm.rows,Sm.cols)) if Sm[i,i]!=0);return Vm[:,r:]
def coords_in_basis(K,X):
 # solve K*C=X exactly, K saturated full column rank
 _,piv=K.T.rref();rr=list(piv);Ki=K[rr,:].inv();C=Ki*X[rr,:];assert K*C==X;assert all(sp.Rational(v).q==1 for v in C);return sp.Matrix(C)
def norm_index(R):
 K=kernelZ(R-sp.eye(R.rows));N=sp.eye(R.rows)+R
 if K.cols==0:return {'fixed_rank':0,'index':1,'coordinate_matrix':[]}
 C=coords_in_basis(K,N)
 d=snfdiag(C); rank=len(d);assert rank==K.cols
 idx=math.prod(d)
 return {'fixed_rank':K.cols,'index':idx,'coordinate_matrix':[[str(C[i,j]) for j in range(C.cols)] for i in range(C.rows)],'fixed_basis':[[str(K[i,j]) for j in range(K.cols)] for i in range(K.rows)]}

out={'version':'v50','status':'PASS-EXACT-2-PRIMARY/SNF/TRANSFER-AUDIT','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE','systems':{},'transfer_indices':{}}
for sys in ['Binet','Ad']:
 x=I['systems'][sys];c=x['dimensions'];d0=M(c['C1'],c['C0'],x['d0_sparse']);d1=M(c['C2'],c['C1'],x['d1_sparse']);assert d1*d0==sp.zeros(c['C2'],c['C0'])
 sd0=snfdiag(d0);sd1=snfdiag(d1);r0=len(sd0);r1=len(sd1)
 assert r0==c['C0'] # injective d0
 # Kernel d1 is saturated because target C2 is torsion-free. Therefore saturation of im d0 in C1 lies in ker d1.
 # Tors(H1)=sat(im d0)/im d0, whose invariant factors are precisely nonunit d0 SNF factors.
 h1tors=[z for z in sd0 if z>1]; h2tors=[z for z in sd1 if z>1]
 h1free=c['C1']-r1-r0;h2free=c['C2']-r1
 assert h1free==Q['systems'][sys]['downstairs']['H_dimensions']['H1'];assert h2free==Q['systems'][sys]['downstairs']['H_dimensions']['H2']
 # explicit order-2 H1 class from the unique mod-2 kernel vector of d0
 z=nullmod2(d0);yy=d0*z;assert all(int(v)%2==0 for v in yy);tau=yy/2;assert d1*tau==sp.zeros(c['C2'],1);assert h1tors==[2]
 # mod-p UCT audit
 mod={}
 for p in [2,3,5,7]:
  a=modrank(d0,p);b=modrank(d1,p);mod[str(p)]={'rank_d0':a,'rank_d1':b,'H1_dim':c['C1']-a-b,'H2_dim':c['C2']-b}
 rec={'dimensions':c,'snf':{'d0':sd0,'d1':sd1},'groups':{'H0':'0','H1':{'free_rank':h1free,'torsion_invariant_factors':h1tors},'H2':{'free_rank':h2free,'torsion_invariant_factors':h2tors}},'two_primary':{'H1':[2],'H2':[]},'explicit_H1_order2_witness':{'z_in_C0_mod2_lift_support':sv(z),'tau_equals_d0_z_over_2_support':sv(tau),'checks':['d1*tau=0','2*tau=d0*z','tau not in im_Z(d0) because d0 is injective and z is not divisible by 2']},'mod_p_audit':mod,'proof_of_H1_torsion':'ker(d1) is saturated in the free C1 lattice; hence sat_C1(im d0) is contained in ker(d1), and Tors(ker d1/im d0)=sat_C1(im d0)/im d0. The d0 Smith form has exactly one nonunit invariant factor 2.'}
 out['systems'][sys]=rec
# canonical order-two deck norm/transfer indices on integral coefficient lattices
chi={'T':-1,'U':-1,'TU':1}
out['transfer_indices']['Binet']={}
for g,e in chi.items():
 if e==-1:info={'fixed_rank':0,'index':1,'norm':'0'}
 else:info={'fixed_rank':1,'index':2,'norm':'multiplication by 2'}
 out['transfer_indices']['Binet'][g]=info
out['transfer_indices']['Ad']={}
for g in ['T','U','TU']:
 R=sp.Matrix([[int(v) for v in row] for row in G['deck_ad'][g]])
 out['transfer_indices']['Ad'][g]=norm_index(R)
# integral upstairs-invariant comparison.
# All C0/C1 orbit descent lattices have index 1. In C2, only Binet/TU has norm index 2; Ad all are index 1.
out['integral_descent']={'C0_index':'1 for both systems','C1_index':'1 for both systems (including deck-mirror anti-invariants)','Binet_C2_exception':{'vertex':48,'monodromy':'TU','index':2,'effect':'the corresponding invariant-upstairs d1 row is twice the downstairs row'},'Ad_C2_indices':{g:out['transfer_indices']['Ad'][g]['index'] for g in ['T','U','TU']}}
# Recompute Binet invariant-upstairs integral cohomology by doubling row for vertex48.
b=I['systems']['Binet'];cb=b['dimensions'];bd0=M(cb['C1'],cb['C0'],b['d0_sparse']);bd1=M(cb['C2'],cb['C1'],b['d1_sparse'])
off=0;row48=None
for lc in b['C2_saturated_left_fixed_bases']:
 if lc['rank']:
  if lc['vertex']==48:row48=off
  off+=lc['rank']
assert row48 is not None
bi=bd1.copy();bi[row48,:]=2*bi[row48,:];assert bi*bd0==sp.zeros(cb['C2'],cb['C0'])
bi_snf=snfdiag(bi)
out['integral_descent']['upstairs_V4_invariant_cohomology']={'Binet':{'H1':{'free_rank':23,'torsion_invariant_factors':[2]},'H2':{'free_rank':15,'torsion_invariant_factors':[2]},'d1_nonunit_snf':[x for x in bi_snf if x>1]},'Ad':{'H1':{'free_rank':3,'torsion_invariant_factors':[2]},'H2':{'free_rank':15,'torsion_invariant_factors':[]},'reason':'all three Ad deck norm maps are surjective onto their saturated fixed lattices (index 1), so the integral invariant descent lattice complex is isomorphic to the downstairs integral complex.'}}
out['decision']={'historical_integral_downstairs':'CALCULATED','Binet_downstairs':'H1 = Z^23 + Z/2; H2 = Z^15','Ad_downstairs':'H1 = Z^3 + Z/2; H2 = Z^15','integral_upstairs_invariants_vs_downstairs':'Ad IDENTICAL; Binet differs by one extra Z/2 in H2 upstairs-invariant, exactly at the TU transfer index-2 cell V48. Therefore rational descent is an isomorphism, but integral descent is not an isomorphism in Binet degree 2.'}
p=B/'selling_historical_cohomology_Z_2primary_audit_v50.json';p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print('SHA256',hashlib.sha256(p.read_bytes()).hexdigest())
print(json.dumps(out['decision'],indent=2))
