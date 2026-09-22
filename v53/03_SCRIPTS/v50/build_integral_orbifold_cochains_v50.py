#!/usr/bin/env python3
import json, math, hashlib
from pathlib import Path
from collections import defaultdict, Counter
import sympy as sp
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
B=Path('/mnt/data/v50_work')
G=json.load(open(B/'selling_field_orbigroupoid_v50.json'))
Links=json.load(open(B/'selling_orbifold_links_v50.json'))['links']
ar={a['id']:a for a in G['arrows']}
def M(x): return sp.Matrix([[sp.Integer(v) for v in row] for row in x])

def kernel_Z(A):
    # Saturated Z-kernel basis as columns, via Smith decomposition U*A*V=S.
    A=sp.Matrix(A)
    dm=DomainMatrix.from_Matrix(A).convert_to(ZZ)
    S,U,V=smith_normal_decomp(dm)
    Sm=S.to_Matrix(); Vm=V.to_Matrix()
    rank=sum(1 for i in range(min(Sm.rows,Sm.cols)) if Sm[i,i]!=0)
    K=Vm[:,rank:]
    assert A*K==sp.zeros(A.rows,K.cols)
    # columns are part of unimodular basis -> saturated
    return K

def coord_map_Z(Bas):
    # Exact rational left inverse; for vectors known to lie in saturated image coordinates must be integral.
    if Bas.cols==0:return sp.zeros(0,Bas.rows)
    # choose independent rows to invert square minor
    _,piv=Bas.T.rref(); rows=list(piv); assert len(rows)==Bas.cols
    Q=Bas[rows,:]
    Qi=Q.inv()
    L=sp.zeros(Bas.cols,Bas.rows)
    for j,r in enumerate(rows): L[:,r]=Qi[:,j]
    assert L*Bas==sp.eye(Bas.cols)
    return L

def clean(d):return {k:sp.simplify(v) for k,v in d.items() if sp.simplify(v)!=0}
def add_block(acc,off,C):
    for i in range(C.rows):
      for j in range(C.cols):
       if C[i,j]:acc[i][off+j]+=C[i,j]
def rows_to_entries(rows):
    return [{'row':i,'col':j,'value':str(int(v))} for i,r in enumerate(rows) for j,v in sorted(r.items()) if v]

def build(kind):
    n=1 if kind=='Binet' else 3
    edgeB={}; left={}; off={}; c1=0; basis_cert=[]
    # C1 lattices
    for e in range(123):
      a=ar[e]
      R=sp.Matrix([[sp.Integer(a['binet_transport'])]]) if kind=='Binet' else M(a['ad_transport'])
      if a['kind'].endswith('ISOTROPY'):
        Bas=kernel_Z(R+sp.eye(n))
      else:
        Bas=sp.eye(n)
      L=coord_map_Z(Bas)
      edgeB[e]=Bas;left[e]=L;off[e]=c1;c1+=Bas.cols
      basis_cert.append({'edge':e,'kind':a['kind'],'rank':Bas.cols,'basis':[[str(Bas[i,j]) for j in range(Bas.cols)] for i in range(n)]})
    # d0
    D0=[defaultdict(lambda:sp.Integer(0)) for _ in range(c1)]
    for e in range(123):
      a=ar[e];s=a['source'];t=a['target']
      R=sp.Matrix([[sp.Integer(a['binet_transport'])]]) if kind=='Binet' else M(a['ad_transport'])
      if edgeB[e].cols:
        Ct=left[e];Cs=left[e]*(-R)
        for i in range(edgeB[e].cols):
          for j in range(n):
            if Ct[i,j]:D0[off[e]+i][n*t+j]+=Ct[i,j]
            if Cs[i,j]:D0[off[e]+i][n*s+j]+=Cs[i,j]
    D0=[clean(x) for x in D0]
    for ri,row in enumerate(D0):
      for c,x in row.items(): assert sp.Rational(x).q==1,(kind,'nonintegral final d0',ri,c,x)
    # d1 using saturated integer left-fixed lattices
    D1=[]; linkcert=[]
    for lk in Links:
      word=lk['chosen']['word']
      Rs=[];Ss=[];Es=[]
      for z in word:
        e=z['edge'];a=ar[e]
        R=sp.Matrix([[sp.Integer(a['binet_transport'])]]) if kind=='Binet' else M(a['ad_transport'])
        d=z['direction']
        if d==1:Rp=R;Sg=sp.eye(n)
        else:Rp=R.inv();Sg=-R.inv()
        # all should remain integral/unimodular
        assert all(sp.Rational(x).q==1 for x in list(Rp)+list(Sg))
        Rs.append(Rp);Ss.append(Sg);Es.append(e)
      tail=sp.eye(n);raw=[defaultdict(lambda:sp.Integer(0)) for _ in range(n)]
      for i in range(len(word)-1,-1,-1):
        e=Es[i];Cblk=tail*Ss[i]*edgeB[e]
        for rr in range(n):
          for j in range(edgeB[e].cols):
            if Cblk[rr,j]:raw[rr][off[e]+j]+=Cblk[rr,j]
        tail=tail*Rs[i]
      H=tail
      K=kernel_Z(H.T-sp.eye(n)) # columns fixed covectors transposed
      P=K.T
      out=[]
      for rr in range(P.rows):
        d=defaultdict(lambda:sp.Integer(0))
        for i in range(n):
          for c,v in raw[i].items():d[c]+=P[rr,i]*v
        d=clean(d)
        assert all(sp.Rational(v).q==1 for v in d.values())
        out.append(d)
      # local d1d0 exact
      for row in out:
        z=defaultdict(lambda:sp.Integer(0))
        for i,a0 in row.items():
          for c,b0 in D0[i].items():z[c]+=a0*b0
        assert not clean(z),(kind,lk['vertex'],clean(z))
      D1.extend(out)
      linkcert.append({'vertex':lk['vertex'],'holonomy':[[str(H[i,j]) for j in range(n)] for i in range(n)],'rank':P.rows,'basis_rows':[[str(P[i,j]) for j in range(n)] for i in range(P.rows)]})
    # dimensions should match rational model
    exp=json.load(open(B/'selling_field_orbifold_cochains_v50.json'))['dimensions'][kind]
    assert c1==exp['C1'] and len(D1)==exp['C2']
    # global product
    for row in D1:
      z=defaultdict(lambda:sp.Integer(0))
      for i,a0 in row.items():
        for c,b0 in D0[i].items():z[c]+=a0*b0
      assert not clean(z)
    return {'dimensions':{'C0':62*n,'C1':c1,'C2':len(D1)},'C1_saturated_bases':basis_cert,'C2_saturated_left_fixed_bases':linkcert,'d0_sparse':rows_to_entries(D0),'d1_sparse':rows_to_entries(D1),'d1d0_zero':True}

out={'version':'v50','status':'PASS-EXACT-INTEGRAL-SATURATED-ORBIFOLD-COCHAINS','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE','lattice_rule':'C1 isotropy lattices and C2 left-fixed lattices are saturated Z-kernels obtained from unimodular Smith decompositions before forming d0,d1.','systems':{}}
for sys in ['Binet','Ad']:
 print('BUILD',sys,flush=True);out['systems'][sys]=build(sys);print(out['systems'][sys]['dimensions'],len(out['systems'][sys]['d0_sparse']),len(out['systems'][sys]['d1_sparse']),flush=True)
p=B/'selling_integral_orbifold_cochains_v50.json';p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('SHA256',hashlib.sha256(p.read_bytes()).hexdigest())
