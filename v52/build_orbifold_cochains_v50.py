#!/usr/bin/env python3
import json, sympy as sp
from pathlib import Path
from collections import Counter
B=Path('/mnt/data/v50_work')
g=json.load(open(B/'selling_field_orbigroupoid_v50.json'))
links=json.load(open(B/'selling_orbifold_links_v50.json'))['links']
reg=json.load(open(B/'selling_regular_TU_quotient_v50.json'))
ar={a['id']:a for a in g['arrows']}
def M(x): return sp.Matrix([[sp.Rational(v) for v in row] for row in x])
# C1 orientation modules: full on ordinary dual edges; anti-invariant under reflection/isotropy on orbifold mirror edges.
edgeB={}; edgeA={}; offB={}; offA={}; c1b=c1a=0
for e in range(123):
 a=ar[e]; eps=sp.Integer(a['binet_transport']); R=M(a['ad_transport'])
 if a['kind'].endswith('ISOTROPY'):
  bb=[] if eps!= -1 else [sp.Matrix([1])]
  aa=(R+sp.eye(3)).nullspace()
 else:
  bb=[sp.Matrix([1])]
  aa=[sp.eye(3)[:,i] for i in range(3)]
 Bb=sp.Matrix.hstack(*bb) if bb else sp.zeros(1,0)
 Ba=sp.Matrix.hstack(*aa) if aa else sp.zeros(3,0)
 edgeB[e]=Bb; edgeA[e]=Ba; offB[e]=c1b; offA[e]=c1a
 c1b+=Bb.cols; c1a+=Ba.cols
# exact coordinate solve via precomputed left inverses
def left_inverse(Bas):
 if Bas.cols==0: return sp.zeros(0,Bas.rows)
 return sp.simplify((Bas.T*Bas).inv()*Bas.T)
leftB={e:left_inverse(edgeB[e]) for e in range(123)}
leftA={e:left_inverse(edgeA[e]) for e in range(123)}
def coords_with(L,Bas,Y):
 if Bas.cols==0:
  assert Y==sp.zeros(Y.rows,Y.cols); return sp.zeros(0,Y.cols)
 C=sp.simplify(L*Y); assert Bas*C==Y; return C
# d0
D0B=sp.zeros(c1b,62); D0A=sp.zeros(c1a,186)
for e in range(123):
 a=ar[e]; s=a['source']; t=a['target']; eps=sp.Integer(a['binet_transport']); R=M(a['ad_transport'])
 Y=sp.zeros(1,62); Y[0,t]+=1; Y[0,s]-=eps
 C=coords_with(leftB[e],edgeB[e],Y)
 if C.rows: D0B[offB[e]:offB[e]+C.rows,:]=C
 YA=sp.zeros(3,186); YA[:,3*t:3*t+3]+=sp.eye(3); YA[:,3*s:3*s+3]-=R
 CA=coords_with(leftA[e],edgeA[e],YA)
 if CA.rows: D0A[offA[e]:offA[e]+CA.rows,:]=CA
# link d1; raw path integral then project to left fixed space of holonomy.
def left_fixed_basis(H):
 # rows lambda with lambda H=lambda => columns of (H.T-I) nullspace, transposed
 ns=(H.T-sp.eye(H.rows)).nullspace()
 return sp.Matrix.vstack(*[v.T for v in ns]) if ns else sp.zeros(0,H.rows)
rowsB=[]; rowsA=[]; linkcert=[]
for lk in links:
 v=lk['vertex']; word=lk['chosen']['word']
 # Build ordered R_i and edge contribution maps y_i from C1 global coordinates to fiber at step target.
 Rs=[]; YBs=[]; YAs=[]; objects=[]
 for z in word:
  e=z['edge']; a=ar[e]; R=M(a['ad_transport']); eps=sp.Integer(a['binet_transport']); d=z['direction']
  if d==1:
   rb=eps; ra=R; sb=sp.Integer(1); sa=sp.eye(3)
  else:
   rb=1/eps; ra=R.inv(); sb=-1/eps; sa=-R.inv()
  # global coefficient selection map into current target fiber
  Yb=sp.zeros(1,c1b)
  if edgeB[e].cols:
   Yb[:,offB[e]:offB[e]+edgeB[e].cols]=sb*edgeB[e]
  Ya=sp.zeros(3,c1a)
  if edgeA[e].cols:
   Ya[:,offA[e]:offA[e]+edgeA[e].cols]=sa*edgeA[e]
  Rs.append((rb,ra)); YBs.append(Yb); YAs.append(Ya)
 # forward-to-final accumulation, process from end backwards
 tailb=sp.Integer(1); taila=sp.eye(3); rawB=sp.zeros(1,c1b); rawA=sp.zeros(3,c1a)
 for i in range(len(word)-1,-1,-1):
  rawB += tailb*YBs[i]; rawA += taila*YAs[i]
  tailb = sp.simplify(tailb*Rs[i][0]); taila = sp.simplify(taila*Rs[i][1])
 HB=sp.Matrix([[tailb]]); HA=taila
 # compare stored holonomy
 assert str(sp.simplify(tailb))==lk['chosen']['binet_holonomy']
 assert HA==M(lk['chosen']['ad_holonomy'])
 PB=left_fixed_basis(HB); PA=left_fixed_basis(HA)
 DB=sp.simplify(PB*rawB); DA=sp.simplify(PA*rawA)
 # local chain identity follows from path telescoping and left-fixed projection; verified globally below.
 assert PB*(sp.eye(1)-HB)==sp.zeros(PB.rows,1)
 assert PA*(sp.eye(3)-HA)==sp.zeros(PA.rows,3)
 rowsB += [DB.row(i) for i in range(DB.rows)]
 rowsA += [DA.row(i) for i in range(DA.rows)]
 linkcert.append({'vertex':v,'word':[{'edge':z['edge'],'direction':z['direction'],'from':z['from_face'],'to':z['to_face']} for z in word],
                  'cycle_count':lk['cycle_count'],'cycle_length':len(word),
                  'Binet_holonomy':str(tailb),'Binet_C2_dim':PB.rows,
                  'Ad_holonomy':[[str(HA[i,j]) for j in range(3)] for i in range(3)],'Ad_C2_dim':PA.rows})
D1B=sp.Matrix.vstack(*rowsB) if rowsB else sp.zeros(0,c1b)
D1A=sp.Matrix.vstack(*rowsA) if rowsA else sp.zeros(0,c1a)
assert D1B*D0B==sp.zeros(D1B.rows,D0B.cols)
assert D1A*D0A==sp.zeros(D1A.rows,D0A.cols)
# save sparse matrices
def sparse_entries(X):
 return [{'row':i,'col':j,'value':str(X[i,j])} for i,j in zip(*X.todok().keys())] if False else [ {'row':i,'col':j,'value':str(v)} for (i,j),v in X.todok().items() ]
out={'version':'v50','status':'PASS-EXACT-ORBIFOLD-COCHAIN-COMPLEX','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
 'carrier':'dual/orbigroupoid of the PROVEN regular primal Selling complex; 62 field objects, 123 dual arrows, 62 primal-event orbifold links',
 'dimensions':{'Binet':{'C0':62,'C1':c1b,'C2':D1B.rows},'Ad':{'C0':186,'C1':c1a,'C2':D1A.rows}},
 'ranks_for_verification_only':'DEFERRED-NOT-NEEDED-FOR-GATE',
 'identities':{'Binet_d1d0_zero':D1B*D0B==sp.zeros(D1B.rows,62),'Ad_d1d0_zero':D1A*D0A==sp.zeros(D1A.rows,186)},
 'C1_edge_dims':{str(e):{'Binet':edgeB[e].cols,'Ad':edgeA[e].cols,'kind':ar[e]['kind']} for e in range(123)},
 'links':linkcert,
 'd0_Binet_sparse':sparse_entries(D0B),'d1_Binet_sparse':sparse_entries(D1B),
 'd0_Ad_sparse':sparse_entries(D0A),'d1_Ad_sparse':sparse_entries(D1A),
 'cohomology_gate':'LOCKED: H1/H2 are not computed or published before F_Sell^field -> F29_can closes object/arrow/link compatibility.'}
json.dump(out,open(B/'selling_field_orbifold_cochains_v50.json','w'),indent=2)
print('dimensions',out['dimensions'])
print('identities',out['identities'])
print('Binet C2 link dims',Counter(x['Binet_C2_dim'] for x in linkcert))
print('Ad C2 link dims',Counter(x['Ad_C2_dim'] for x in linkcert))
