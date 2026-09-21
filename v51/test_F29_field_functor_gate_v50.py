#!/usr/bin/env python3
import json,itertools
from collections import deque,Counter
from pathlib import Path
import numpy as np, sympy as sp, sys
B=Path('/mnt/data/v50_work')
sys.path.insert(0,str(B))
from selling_reconstruction_engine_v50 import S,RELABEL
Gd=json.load(open(B/'selling_field_orbigroupoid_v50.json'))
links=json.load(open(B/'selling_orbifold_links_v50.json'))['links']
bfs=json.load(open(B/'exact_TU_field_bfs_v50.json'))
glob=json.load(open(B/'selling_global_wall_linearization_v50.json'))
ar={a['id']:a for a in Gd['arrows']}; rec={x['id']:x for x in bfs['records']}
first={'g':9,'h':2,'k':1,'l':3,'m':5,'n':4}
def n2(M):return np.array([[int(x)%2 for x in row] for row in M],dtype=np.uint8)
def sp2(M):return np.array(M.tolist(),dtype=int)%2
def key(A):return tuple(int(x) for x in np.array(A,dtype=np.uint8).reshape(-1))
def rank2(A):
 A=np.array(A,dtype=np.uint8).copy()%2;r=0
 for c in range(A.shape[1]):
  p=next((i for i in range(r,A.shape[0]) if A[i,c]),None)
  if p is not None:
   A[[r,p]]=A[[p,r]]
   for i in range(A.shape[0]):
    if i!=r and A[i,c]:A[i]^=A[r]
   r+=1
 return r
# GL(3,2)
G=[]
for bits in itertools.product([0,1],repeat=9):
 A=np.array(bits,dtype=np.uint8).reshape(3,3)
 if rank2(A)==3:G.append(A)
G.sort(key=key);gid={key(A):i for i,A in enumerate(G)};I=np.eye(3,dtype=np.uint8)
assert len(G)==168
S2=[sp2(x) for x in S]
# shortest words from identity for every group element (right Cayley)
root=gid[key(I)]; word={root:[]};dq=deque([root])
while dq:
 u=dq.popleft()
 for j,s in enumerate(S2,1):
  v=gid[key(G[u]@s%2)]
  if v not in word:word[v]=word[u]+[j];dq.append(v)
assert len(word)==168
# deck V4 in Q2 chart
GT=n2(glob['deck_in_Q2_chart']['G_T']);GU=n2(glob['deck_in_Q2_chart']['G_U'])
V4=[I,GT,GU,GT@GU%2]
assert len({key(x) for x in V4})==4
# left action free and 42 orbits
assert all(not np.array_equal(d@A%2,A) for d in V4[1:] for A in G)
orb={};orbits=[]
for A in G:
 k=min(key(d@A%2) for d in V4)
 if k not in orb:orb[k]=len(orbits);orbits.append(k)
assert len(orbits)==42
# source object charts A_f mod2
Af={r['id']:n2(r['path_matrix']) for r in bfs['records']}
# arrow mod2 matrices in field coordinates and exact endpoint deck witness
arr=[];strict=0;quot=0
for a in Gd['arrows']:
 e=a['id'];s=a['source'];t=a['target']
 if a['kind']=='ORDINARY-CROSSING':
  M=S[first[a['wall']]-1]*RELABEL[a.get('relabel_index',0)].inv(); mm=sp2(M)
 elif a['kind'] in ('MULTIWALL-ISOTROPY','DECK-MIRROR-ISOTROPY'):
  mm=n2(a['field_matrix'])
 else:raise RuntimeError(a['kind'])
 start=Af[s];end=start@mm%2;target=Af[t]
 ds=[i for i,d in enumerate(V4) if np.array_equal(end,d@target%2)]
 okq=len(ds)>0; oks=np.array_equal(end,target)
 assert okq,(e,s,t,end,target)
 quot+=okq;strict+=oks
 wid=gid[key(mm)]; w=word[wid]
 arr.append({'edge':e,'kind':a['kind'],'source':s,'target':t,'mod2_matrix':mm.astype(int).tolist(),'shortest_S_word':['S'+str(j) for j in w],
             'strict_F29_endpoint':bool(oks),'V4_endpoint_witness':['I','T','U','TU'][ds[0]],'quotient42_endpoint':True})
# link closures under arrow matrices
linkout=[]; strictlinks=0
for lk in links:
 w=lk['chosen']['word']; startf=w[0]['from_face']; cur=Af[startf].copy()
 for z in w:
  aa=arr[z['edge']]; mm=np.array(aa['mod2_matrix'],dtype=np.uint8)
  if z['direction']==-1:
   # inverse by search; matrices here generally involutive? use group inverse exact by lookup brute
   inv=next(X for X in G if np.array_equal(mm@X%2,I));mm=inv
  cur=cur@mm%2
 strictclose=np.array_equal(cur,Af[startf])
 ds=[i for i,d in enumerate(V4) if np.array_equal(cur,d@Af[startf]%2)]
 assert ds
 strictlinks+=strictclose
 linkout.append({'vertex':lk['vertex'],'strict_F29_closed':bool(strictclose),'V4_closure':['I','T','U','TU'][ds[0]],'quotient42_closed':True,
                 'F29_fillable_by_841_cycle_basis':bool(strictclose),'reason_if_not_strict':None if strictclose else 'endpoint differs by nontrivial left V4 deck action'})
# explicit no-go witness edge 118
witness=next(x for x in arr if x['edge']==118)
assert not witness['strict_F29_endpoint'] and witness['V4_endpoint_witness']!='I'
out={'version':'v50','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
 'strict_field_functor_to_F29_can':{'status':'REFUTED-TYPED','reason':'F29_can is the free right GL(3,2) Cayley/congruence chart. Source orbifold loops with nonidentity mod-2 transport cannot map to endomorphism paths at a single target object while preserving reduction mod 2.',
   'object_count_target':168,'arrow_endpoint_matches_strict':strict,'arrow_total':123,'link_loops_strictly_closed':strictlinks,'link_total':62,'no_go_witness_edge_118':witness,
   'free_right_action_exhaustive_check_168x3_nonidentity_V4':True},
 'canonical_repair':{'target':'left-V4 quotient/orbifold of F29_can','status':'PROVEN-OBJECT/ARROW/LINK-COMPATIBLE-MOD2','V4_order':4,'object_orbits':42,
   'all_123_arrow_endpoints_close_in_quotient':quot==123,'all_62_orbifold_links_close_in_quotient':all(x['quotient42_closed'] for x in linkout),
   'interpretation':'Repeated-field T/U identifications are retained as orbifold/deck data instead of being forced to disappear inside the free 168-object chart.'},
 'arrows':arr,'links':linkout,
 'cohomology_gate':'LOCKED: the requested strict F_Sell^field -> F29_can gate fails. Do not compute historical H1/H2. A successor may retarget to the certified left-V4 quotient/orbifold or lift Selling back to its prequotient development before mapping to F29_can.'}
json.dump(out,open(B/'selling_F29_field_functor_gate_v50.json','w'),indent=2)
print('strict arrows',strict,'/123; strict links',strictlinks,'/62; quotient arrows',quot,'/123; quotient links',sum(x['quotient42_closed'] for x in linkout),'/62')
print('witness',witness)
print('link V4 closures',Counter(x['V4_closure'] for x in linkout))
