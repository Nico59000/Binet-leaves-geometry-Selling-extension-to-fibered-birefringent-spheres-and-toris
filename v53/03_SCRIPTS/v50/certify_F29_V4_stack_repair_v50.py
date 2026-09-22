#!/usr/bin/env python3
import json,itertools
from collections import deque,Counter
from pathlib import Path
import numpy as np,sys
B=Path('/mnt/data/v50_work')
F=json.load(open(B/'selling_F29_field_functor_gate_v50.json'))
L=json.load(open(B/'selling_orbifold_links_v50.json'))['links']
BFS=json.load(open(B/'exact_TU_field_bfs_v50.json'))
GLOB=json.load(open(B/'selling_global_wall_linearization_v50.json'))
sys.path.insert(0,str(B))
from selling_reconstruction_engine_v50 import S

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
G=[]
for bits in itertools.product([0,1],repeat=9):
 A=np.array(bits,dtype=np.uint8).reshape(3,3)
 if rank2(A)==3:G.append(A)
G.sort(key=key);gid={key(A):i for i,A in enumerate(G)};I=np.eye(3,dtype=np.uint8)
S2=[np.array(M.tolist(),dtype=int)%2 for M in S]
root=gid[key(I)]
# BFS tree and shortest generator words
parent=[None]*168; pgen=[None]*168; paths=[None]*168
parent[root]=root;paths[root]=[];dq=deque([root])
while dq:
 u=dq.popleft()
 for j,s in enumerate(S2,1):
  v=gid[key(G[u]@s%2)]
  if parent[v] is None:parent[v]=u;pgen[v]=j;paths[v]=paths[u]+[j];dq.append(v)
assert all(x is not None for x in parent)
# all undirected Cayley edges and tree edges
edges=[]
for u,A in enumerate(G):
 for j,s in enumerate(S2,1):
  v=gid[key(A@s%2)]
  if u<v:edges.append((u,v,j))
tree=set()
for v in range(168):
 if v==root:continue
 u=parent[v];j=pgen[v];tree.add((min(u,v),max(u,v),j))
non=[e for e in edges if e not in tree];cell_index={e:i for i,e in enumerate(non)}
assert len(non)==841
# shortest word for arbitrary matrix uses paths[gid]
def inv2(A):
 return next(X for X in G if np.array_equal(A@X%2,I))
# unique tree vertex path u->v as list steps (a,b,j)
def root_chain(v):
 arr=[]
 while v!=root:
  u=parent[v];j=pgen[v];arr.append((u,v,j));v=u
 arr.reverse();return arr # root->orig
def tree_path(u,v):
 cu=root_chain(u);cv=root_chain(v)
 # vertices root paths; find common prefix edge count
 k=0
 while k<min(len(cu),len(cv)) and cu[k]==cv[k]:k+=1
 out=[]
 # u back to LCA: reverse cu[k:]
 for a,b,j in reversed(cu[k:]):out.append((b,a,j))
 for a,b,j in cv[k:]:out.append((a,b,j))
 return out
# reduce a tree walk by cancelling immediate inverse edge traversals
def reduce_tree_walk(steps):
 st=[]
 for a,b,j in steps:
  if st and st[-1][0]==b and st[-1][1]==a and st[-1][2]==j:st.pop()
  else:st.append((a,b,j))
 return st
# arrow map edge->word ints
aw={a['edge']:[int(x[1:]) for x in a['shortest_S_word']] for a in F['arrows']}
Af={r['id']:np.array(r['path_matrix'],dtype=int)%2 for r in BFS['records']}
# V4 for labels
GT=np.array(GLOB['deck_in_Q2_chart']['G_T'],dtype=int)%2;GU=np.array(GLOB['deck_in_Q2_chart']['G_U'],dtype=int)%2
V4={'I':I,'T':GT,'U':GU,'TU':GT@GU%2}
cert=[]
for lk,fg in zip(L,F['links']):
 assert lk['vertex']==fg['vertex']
 word=lk['chosen']['word'];startf=word[0]['from_face'];start=Af[startf];cur=start.copy();steps=[]
 for z in word:
  ww=aw[z['edge']] if z['direction']==1 else list(reversed(aw[z['edge']]))
  for j in ww:
   u=gid[key(cur)];cur=cur@S2[j-1]%2;v=gid[key(cur)];steps.append((u,v,j))
 # check closure label
 dlab=fg['V4_closure'];assert np.array_equal(cur,V4[dlab]@start%2)
 # append canonical connector to raw start if nontrivial deck closure
 conn=inv2(cur)@start%2
 cw=paths[gid[key(conn)]]
 connector=[]
 for j in cw:
  u=gid[key(cur)];cur=cur@S2[j-1]%2;v=gid[key(cur)];connector.append((u,v,j));steps.append((u,v,j))
 assert np.array_equal(cur,start)
 # fundamental cell fill: replace every non-tree chord by tree path
 fill=[];rewritten=[]
 for u,v,j in steps:
  e=(min(u,v),max(u,v),j)
  if e in tree:rewritten.append((u,v,j))
  else:
   idx=cell_index[e];sign=1 if (u,v,j)==e else -1
   fill.append({'fundamental_cell':idx,'sign':sign,'chord':[u,v,j]})
   rewritten.extend(tree_path(u,v))
 red=reduce_tree_walk(rewritten)
 assert not red,(lk['vertex'],red)
 cert.append({'vertex':lk['vertex'],'V4_closure':dlab,'source_link_step_count':len(steps)-len(connector),'connector_word':['S'+str(j) for j in cw],
              'closed_F29_step_count':len(steps),'fundamental_cell_occurrences':fill,'fundamental_cell_count':len(fill),'tree_rewrite_reduces_to_empty':True})
out={'version':'v50','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE','status':'PASS-EXACT-V4-STACK-REPAIR-2CELL-COMPATIBILITY',
 'target':'[F29_can / V4_left] action-groupoid/stack completion; not the strict free 168-object F29_can chart',
 'V4_closure_histogram':dict(Counter(x['V4_closure'] for x in cert)),'links':cert,
 'checks':{'all_62_links_close_after_canonical_V4_connector':True,'all_62_closed_lifts_fill_by_841_fundamental_cells':True,'F29_fundamental_cell_count':841},
 'interpretation':'For 59 links the connector is constant. For the three deck links, the connector records T, U or TU isotropy. After adding that stack isotropy, the lifted loop is filled constructively by the existing 841 F29 fundamental 2-cells.',
 'strict_functor_status':'REFUTED-TYPED remains unchanged; this is a typed target repair, not a promotion of F_Sell^field -> F29_can.',
 'cohomology_gate':'LOCKED pending a decision to accept the repaired stack target as the intended comparison object and a naturality audit for Binet/Ad transports.'}
json.dump(out,open(B/'selling_F29_V4_stack_repair_v50.json','w'),indent=2)
print('PASS',out['V4_closure_histogram'],'fill cells total',sum(x['fundamental_cell_count'] for x in cert),'max',max(x['fundamental_cell_count'] for x in cert))
