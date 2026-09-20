#!/usr/bin/env python3
import json, sympy as sp
from pathlib import Path
from collections import Counter,defaultdict
B=Path('/mnt/data/v50_work')
g=json.load(open(B/'selling_field_orbigroupoid_v50.json'))
links=json.load(open(B/'selling_orbifold_links_v50.json'))['links']
ar={a['id']:a for a in g['arrows']}
def M(x): return sp.Matrix([[sp.Rational(v) for v in row] for row in x])
def left_inverse(Bas):
    if Bas.cols==0:return sp.zeros(0,Bas.rows)
    return (Bas.T*Bas).inv()*Bas.T
# C1 orientation modules
edgeB={};edgeA={};leftB={};leftA={};offB={};offA={};c1b=c1a=0
for e in range(123):
    a=ar[e]; eps=sp.Integer(a['binet_transport']); R=M(a['ad_transport'])
    if a['kind'].endswith('ISOTROPY'):
        bb=[] if eps!=-1 else [sp.Matrix([1])]
        aa=(R+sp.eye(3)).nullspace()
    else:
        bb=[sp.Matrix([1])]; aa=[sp.eye(3)[:,i] for i in range(3)]
    Bb=sp.Matrix.hstack(*bb) if bb else sp.zeros(1,0)
    Ba=sp.Matrix.hstack(*aa) if aa else sp.zeros(3,0)
    edgeB[e]=Bb;edgeA[e]=Ba;leftB[e]=left_inverse(Bb);leftA[e]=left_inverse(Ba)
    offB[e]=c1b;offA[e]=c1a;c1b+=Bb.cols;c1a+=Ba.cols
# sparse rows d0: list dict col->val
D0B=[defaultdict(lambda:sp.Integer(0)) for _ in range(c1b)]
D0A=[defaultdict(lambda:sp.Integer(0)) for _ in range(c1a)]
for e in range(123):
    a=ar[e];s=a['source'];t=a['target'];eps=sp.Integer(a['binet_transport']);R=M(a['ad_transport'])
    if edgeB[e].cols:
        Ct=leftB[e]; Cs=leftB[e]*(-sp.Matrix([[eps]]))
        for i in range(edgeB[e].cols):
            D0B[offB[e]+i][t]+=Ct[i,0];D0B[offB[e]+i][s]+=Cs[i,0]
    if edgeA[e].cols:
        Ct=leftA[e]; Cs=leftA[e]*(-R)
        for i in range(edgeA[e].cols):
            for j in range(3):
                if Ct[i,j]:D0A[offA[e]+i][3*t+j]+=Ct[i,j]
                if Cs[i,j]:D0A[offA[e]+i][3*s+j]+=Cs[i,j]
# clean zeros
def clean(d): return {k:sp.simplify(v) for k,v in d.items() if sp.simplify(v)!=0}
D0B=[clean(x) for x in D0B];D0A=[clean(x) for x in D0A]
def left_fixed(H):
    ns=(H.T-sp.eye(H.rows)).nullspace()
    return sp.Matrix.vstack(*[v.T for v in ns]) if ns else sp.zeros(0,H.rows)
def add_block(acc, baseoff, C):
    # C d x k; acc list of d dicts
    for i in range(C.rows):
        for j in range(C.cols):
            if C[i,j]: acc[i][baseoff+j]+=C[i,j]
# d1 rows and certs
D1B=[];D1A=[];cert=[];rowoffB=rowoffA=0
for lk in links:
    word=lk['chosen']['word']
    # compute suffix transports first
    Rb=[];Ra=[];Sb=[];Sa=[];Es=[]
    for z in word:
        e=z['edge'];a=ar[e];eps=sp.Integer(a['binet_transport']);R=M(a['ad_transport']);d=z['direction']
        if d==1: rb=eps;ra=R;sb=sp.Integer(1);sa=sp.eye(3)
        else: rb=1/eps;ra=R.inv();sb=-1/eps;sa=-R.inv()
        Rb.append(rb);Ra.append(ra);Sb.append(sb);Sa.append(sa);Es.append(e)
    tailb=sp.Integer(1);taila=sp.eye(3)
    rawB=[defaultdict(lambda:sp.Integer(0))]
    rawA=[defaultdict(lambda:sp.Integer(0)) for _ in range(3)]
    for i in range(len(word)-1,-1,-1):
        e=Es[i]
        if edgeB[e].cols:
            C=sp.Matrix([[tailb*Sb[i]]])*edgeB[e]
            add_block(rawB,offB[e],C)
        if edgeA[e].cols:
            C=taila*Sa[i]*edgeA[e]
            add_block(rawA,offA[e],C)
        tailb=sp.simplify(tailb*Rb[i]);taila=taila*Ra[i]
    HB=sp.Matrix([[tailb]]);HA=taila
    assert str(tailb)==lk['chosen']['binet_holonomy']; assert HA==M(lk['chosen']['ad_holonomy'])
    PB=left_fixed(HB);PA=left_fixed(HA)
    # project sparse raw maps
    outB=[]
    for rr in range(PB.rows):
        d=defaultdict(lambda:sp.Integer(0))
        for i in range(1):
            for c,vv in rawB[i].items(): d[c]+=PB[rr,i]*vv
        outB.append(clean(d))
    outA=[]
    for rr in range(PA.rows):
        d=defaultdict(lambda:sp.Integer(0))
        for i in range(3):
            for c,vv in rawA[i].items(): d[c]+=PA[rr,i]*vv
        outA.append(clean(d))
    # local exact chain check via sparse multiplication
    def check_rows(rows,D0):
        for row in rows:
            z=defaultdict(lambda:sp.Integer(0))
            for i,a0 in row.items():
                for c,b0 in D0[i].items(): z[c]+=a0*b0
            z=clean(z)
            assert not z,z
    check_rows(outB,D0B);check_rows(outA,D0A)
    D1B.extend(outB);D1A.extend(outA)
    cert.append({'vertex':lk['vertex'],'word':[{'edge':z['edge'],'direction':z['direction'],'from':z['from_face'],'to':z['to_face']} for z in word],
                 'cycle_count':lk['cycle_count'],'cycle_length':len(word),'Binet_holonomy':str(tailb),'Binet_C2_dim':PB.rows,
                 'Ad_holonomy':[[str(HA[i,j]) for j in range(3)] for i in range(3)],'Ad_C2_dim':PA.rows,
                 'd1d0_local_Binet':True,'d1d0_local_Ad':True})
# global exact sparse check (local already sufficient, but replay globally)
def global_product_zero(D1,D0):
    for row in D1:
        z=defaultdict(lambda:sp.Integer(0))
        for i,a0 in row.items():
            for c,b0 in D0[i].items():z[c]+=a0*b0
        if clean(z): return False
    return True
assert global_product_zero(D1B,D0B);assert global_product_zero(D1A,D0A)
def entries(rows):
    return [{'row':i,'col':j,'value':str(v)} for i,r in enumerate(rows) for j,v in sorted(r.items())]
out={'version':'v50','status':'PASS-EXACT-ORBIFOLD-COCHAIN-COMPLEX','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
     'carrier':'dual/orbigroupoid of X_Sell^cell regular quotient',
     'dimensions':{'Binet':{'C0':62,'C1':c1b,'C2':len(D1B)},'Ad':{'C0':186,'C1':c1a,'C2':len(D1A)}},
     'identities':{'Binet_d1d0_zero':True,'Ad_d1d0_zero':True},
     'C1_orientation_module_rule':'ordinary crossing: full coefficient fiber; MULTIWALL/DECK mirror isotropy: anti-invariant subspace ker(rho+I)',
     'C2_link_module_rule':'for each exact orbifold-link word with coefficient holonomy H, use left fixed covectors ker(H^T-I); d1 is the transported path sum projected to this space',
     'link_dimension_histograms':{'Binet':dict(Counter(x['Binet_C2_dim'] for x in cert)),'Ad':dict(Counter(x['Ad_C2_dim'] for x in cert))},
     'links':cert,'d0_Binet_sparse':entries(D0B),'d1_Binet_sparse':entries(D1B),'d0_Ad_sparse':entries(D0A),'d1_Ad_sparse':entries(D1A),
     'cohomology_gate':'LOCKED: no H1/H2 until the complete F_Sell^field -> F29_can object/arrow/link functor is certified.'}
json.dump(out,open(B/'selling_field_orbifold_cochains_v50.json','w'),indent=2)
print('PASS dims',out['dimensions'])
print('hist',out['link_dimension_histograms'])
print('nnz',len(out['d0_Binet_sparse']),len(out['d1_Binet_sparse']),len(out['d0_Ad_sparse']),len(out['d1_Ad_sparse']))
