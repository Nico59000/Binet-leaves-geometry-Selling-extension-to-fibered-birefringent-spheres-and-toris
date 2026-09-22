#!/usr/bin/env python3
import json,time,hashlib
from pathlib import Path
import sympy as sp
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp, smith_normal_form
B=Path('/mnt/data/v50_work')
D=json.load(open(B/'selling_integral_orbifold_cochains_v50.json'))
def M(rows,cols,E):return sp.SparseMatrix(rows,cols,{(e['row'],e['col']):sp.Integer(e['value']) for e in E})
def dmZ(A):return DomainMatrix.from_Matrix(A).convert_to(ZZ)
def diag_nonzero(S):
 m=S.to_Matrix() if hasattr(S,'to_Matrix') else S
 return [abs(int(m[i,i])) for i in range(min(m.rows,m.cols)) if m[i,i]!=0]
def snf_summary(A):
 S=smith_normal_form(dmZ(A));diag=diag_nonzero(S);return diag
out={'version':'v50','status':'PASS-EXACT-INTEGRAL-HISTORICAL-COHOMOLOGY','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE','systems':{}}
for sys in ['Binet','Ad']:
 print('SYS',sys,flush=True)
 x=D['systems'][sys];c=x['dimensions'];d0=M(c['C1'],c['C0'],x['d0_sparse']);d1=M(c['C2'],c['C1'],x['d1_sparse'])
 assert d1*d0==sp.zeros(c['C2'],c['C0'])
 print('smith decomp d1',d1.shape,flush=True);t=time.time(); S,U,V=smith_normal_decomp(dmZ(d1));print('sec',time.time()-t,flush=True)
 Sm=S.to_Matrix();Vm=V.to_Matrix();Um=U.to_Matrix(); diag1=diag_nonzero(S);r1=len(diag1)
 assert Um*d1*Vm==Sm
 assert all(q==1 for q in diag1),diag1
 # V unimodular, columns r1: give saturated integral kernel of d1
 K=Vm[:,r1:]; assert d1*K==sp.zeros(c['C2'],K.cols)
 # V^{-1} d0: first r1 coordinates must vanish; bottom are coordinates in K basis
 Vinv=Vm.inv(); assert all(z.q==1 for z in Vinv)
 C=Vinv*d0
 assert C[:r1,:]==sp.zeros(r1,c['C0'])
 Z=C[r1:,:]
 assert K*Z==d0
 print('kernel',K.shape,'coord',Z.shape,'snf H1 relation',flush=True)
 sz=smith_normal_form(dmZ(Z)); dz=diag_nonzero(sz); rz=len(dz)
 # H1 = Z^(ker rank)/im Z ; torsion dz>1 + free ker-rz
 h1free=K.cols-rz; h1tors=[q for q in dz if q>1]
 # H2 coker d1, since diag1 all 1 => free C2-r1
 h2free=c['C2']-r1; h2tors=[q for q in diag1 if q>1]
 # H0 = ker d0. rational ranks say zero; exact d0 rank columns full.
 d0diag=diag_nonzero(smith_normal_form(dmZ(d0))); h0free=c['C0']-len(d0diag)
 # 2-primary summary
 def v2(n):
  a=0
  while n%2==0 and n:a+=1;n//=2
  return a
 rec={'dimensions':c,'snf_d0_nonunit':[q for q in d0diag if q>1],'snf_d1_nonunit':[q for q in diag1 if q>1],
      'kernel_d1_rank':K.cols,'H0_Z':{'free_rank':h0free,'torsion':[]},
      'H1_Z':{'free_rank':h1free,'torsion_invariant_factors':h1tors},
      'H2_Z':{'free_rank':h2free,'torsion_invariant_factors':h2tors},
      'H1_relation_snf_nonunit':[q for q in dz if q>1],
      'two_primary':{'H1_2_primary':[2**v2(q) for q in h1tors if v2(q)>0],'H2_2_primary':[2**v2(q) for q in h2tors if v2(q)>0]},
      'checks':{'d1d0_zero':True,'d1_image_primitive':all(q==1 for q in diag1),'kernel_basis_saturated_from_unimodular_V':True,'d0_coordinates_integral_in_kernel_basis':True}}
 out['systems'][sys]=rec
 print('RESULT',rec['H1_Z'],rec['H2_Z'],flush=True)
p=B/'selling_historical_cohomology_Z_v50.json';p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print('sha',hashlib.sha256(p.read_bytes()).hexdigest())
