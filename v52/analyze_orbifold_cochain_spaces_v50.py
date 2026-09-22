#!/usr/bin/env python3
import json, sympy as sp
from pathlib import Path
from collections import Counter
B=Path('/mnt/data/v50_work')
g=json.load(open(B/'selling_field_orbigroupoid_v50.json'))
reg=json.load(open(B/'selling_regular_TU_quotient_v50.json'))

def M(x): return sp.Matrix([[sp.Rational(v) for v in row] for row in x])
# edge bases for C1
edgeB_B={}; edgeB_A={}; edgeoffB={}; edgeoffA={}; c1b=c1a=0
for a in g['arrows']:
    e=a['id']; kind=a['kind']; eps=sp.Integer(a['binet_transport']); R=M(a['ad_transport'])
    if kind.endswith('ISOTROPY'):
        # orientation-reversing edge stabilizer: rho v = -v
        bb=[] if eps!= -1 else [sp.Matrix([1])]
        aa=(R+sp.eye(3)).nullspace()
    else:
        bb=[sp.Matrix([1])]
        aa=[sp.eye(3)[:,i] for i in range(3)]
    Bb=sp.Matrix.hstack(*bb) if bb else sp.zeros(1,0)
    Ba=sp.Matrix.hstack(*aa) if aa else sp.zeros(3,0)
    edgeB_B[e]=Bb; edgeB_A[e]=Ba
    edgeoffB[e]=c1b; edgeoffA[e]=c1a
    c1b+=Bb.cols; c1a+=Ba.cols
# helper coordinates in basis with exact check
def coords(Bas,Y):
    if Bas.cols==0:
        assert Y==sp.zeros(Y.rows,Y.cols)
        return sp.zeros(0,Y.cols)
    # solve Bas*C=Y columnwise
    C=sp.zeros(Bas.cols,Y.cols)
    for j in range(Y.cols):
        sol=sp.linsolve((Bas,Y[:,j])); vv=sp.Matrix(list(next(iter(sol))))
        assert Bas*vv==Y[:,j]
        C[:,j]=vv
    return C
D0B=sp.zeros(c1b,62)
D0A=sp.zeros(c1a,62*3)
for a in g['arrows']:
    e=a['id']; s=a['source']; t=a['target']; eps=sp.Integer(a['binet_transport']); R=M(a['ad_transport'])
    Yb=sp.zeros(1,62); Yb[0,t]+=1; Yb[0,s]-=eps
    Cb=coords(edgeB_B[e],Yb)
    if edgeB_B[e].cols: D0B[edgeoffB[e]:edgeoffB[e]+edgeB_B[e].cols,:]=Cb
    Ya=sp.zeros(3,62*3); Ya[:,3*t:3*t+3]+=sp.eye(3); Ya[:,3*s:3*s+3]-=R
    Ca=coords(edgeB_A[e],Ya)
    if edgeB_A[e].cols: D0A[edgeoffA[e]:edgeoffA[e]+edgeB_A[e].cols,:]=Ca
# local relations
local=[]
for v in range(reg['C0']):
    es=[e['id'] for e in reg['edges'] if v in e['endpoints']]
    fs=sorted(set(sum(([a['source'],a['target']] for a in g['arrows'] if a['id'] in es),[])))
    rb=[]
    for e in es: rb += list(range(edgeoffB[e],edgeoffB[e]+edgeB_B[e].cols))
    ra=[]
    for e in es: ra += list(range(edgeoffA[e],edgeoffA[e]+edgeB_A[e].cols))
    cb=fs
    ca=[]
    for f in fs: ca+=list(range(3*f,3*f+3))
    Mb=D0B.extract(rb,cb) if rb else sp.zeros(0,len(cb))
    Ma=D0A.extract(ra,ca) if ra else sp.zeros(0,len(ca))
    nb=Mb.T.nullspace(); na=Ma.T.nullspace()
    local.append({'vertex':v,'edges':es,'faces':fs,'C1B_local_dim':len(rb),'C1A_local_dim':len(ra),'Binet_relation_dim':len(nb),'Ad_relation_dim':len(na)})
# serialize bases and summary
out={'version':'v50','status':'PASS-EXACT-ORBICHAIN-SPACE-DIAGNOSTIC','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
     'C0_dims':{'Binet':62,'Ad':186},'C1_dims':{'Binet':c1b,'Ad':c1a},
     'C1_edge_dims':{str(e):{'Binet':edgeB_B[e].cols,'Ad':edgeB_A[e].cols} for e in range(123)},
     'd0_ranks':{'Binet':D0B.rank(),'Ad':D0A.rank()},
     'local_relation_dims':local,
     'relation_dim_histograms':{'Binet':dict(Counter(x['Binet_relation_dim'] for x in local)),'Ad':dict(Counter(x['Ad_relation_dim'] for x in local))},
     'guard':'Local left-kernel dimensions are used to discover the orbifold link coefficient spaces; d1 is not promoted until each basis is matched to a certified cyclic link relation.'}
json.dump(out,open(B/'selling_orbichain_space_diagnostic_v50.json','w'),indent=2)
print('C1 dims',c1b,c1a,'ranks',D0B.rank(),D0A.rank())
print('hist B',Counter(x['Binet_relation_dim'] for x in local))
print('hist A',Counter(x['Ad_relation_dim'] for x in local))
for x in local:
    if x['Binet_relation_dim'] not in (0,1) or x['Ad_relation_dim'] not in (0,1,2,3): print('odd',x)
