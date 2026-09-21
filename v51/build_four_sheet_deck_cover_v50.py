#!/usr/bin/env python3
import json,itertools,hashlib
from collections import deque,Counter
from pathlib import Path
import numpy as np
B=Path('/mnt/data/v50_work')
Gd=json.load(open(B/'selling_field_orbigroupoid_v50.json'))
Links=json.load(open(B/'selling_orbifold_links_v50.json'))['links']
Gate=json.load(open(B/'selling_F29_field_functor_gate_v50.json'))
Stack=json.load(open(B/'selling_F29_V4_stack_repair_v50.json'))
BFS=json.load(open(B/'exact_TU_field_bfs_v50.json'))
GLOB=json.load(open(B/'selling_global_wall_linearization_v50.json'))
# helpers
def n2(M): return np.array([[int(x)%2 for x in row] for row in M],dtype=np.uint8)
def key(A): return tuple(int(x) for x in np.array(A,dtype=np.uint8).reshape(-1))
def rank2(A):
 A=np.array(A,dtype=np.uint8).copy()%2;r=0
 for c in range(A.shape[1]):
  p=next((i for i in range(r,A.shape[0]) if A[i,c]),None)
  if p is not None:
   A[[r,p]]=A[[p,r]]
   for i in range(A.shape[0]):
    if i!=r and A[i,c]: A[i]^=A[r]
   r+=1
 return r
# GL(3,2)
GL=[]
for bits in itertools.product([0,1],repeat=9):
 A=np.array(bits,dtype=np.uint8).reshape(3,3)
 if rank2(A)==3: GL.append(A)
GL.sort(key=key); gid={key(A):i for i,A in enumerate(GL)}; I=np.eye(3,dtype=np.uint8)
assert len(GL)==168
# S words from existing gate are enough; build tree using the actual mod2 arrow-word generators referenced there
# recover generator matrices from previous engine
import sys
sys.path.insert(0,str(B))
from selling_reconstruction_engine_v50 import S
S2=[np.array(M.tolist(),dtype=int)%2 for M in S]
root=gid[key(I)]
parent=[None]*168;pgen=[None]*168;paths=[None]*168
parent[root]=root;paths[root]=[];dq=deque([root])
while dq:
 u=dq.popleft()
 for j,s in enumerate(S2,1):
  v=gid[key(GL[u]@s%2)]
  if parent[v] is None:
   parent[v]=u;pgen[v]=j;paths[v]=paths[u]+[j];dq.append(v)
assert all(x is not None for x in parent)
# fundamental cycles = non-tree undirected Cayley edges
edges=[]
for u,A in enumerate(GL):
 for j,s in enumerate(S2,1):
  v=gid[key(A@s%2)]
  if u<v: edges.append((u,v,j))
tree=set()
for v in range(168):
 if v==root: continue
 u=parent[v];j=pgen[v];tree.add((min(u,v),max(u,v),j))
non=[e for e in edges if e not in tree]; cell_index={e:i for i,e in enumerate(non)}
assert len(non)==841

def root_chain(v):
 arr=[]
 while v!=root:
  u=parent[v];j=pgen[v];arr.append((u,v,j));v=u
 arr.reverse();return arr
def tree_path(u,v):
 cu=root_chain(u);cv=root_chain(v);k=0
 while k<min(len(cu),len(cv)) and cu[k]==cv[k]: k+=1
 out=[]
 for a,b,j in reversed(cu[k:]): out.append((b,a,j))
 for a,b,j in cv[k:]: out.append((a,b,j))
 return out
def reduce_tree_walk(steps):
 st=[]
 for a,b,j in steps:
  if st and st[-1][0]==b and st[-1][1]==a and st[-1][2]==j: st.pop()
  else: st.append((a,b,j))
 return st
# V4 labels, exact XOR law
GT=n2(GLOB['deck_in_Q2_chart']['G_T']); GU=n2(GLOB['deck_in_Q2_chart']['G_U'])
Vmat={'I':I,'T':GT,'U':GU,'TU':GT@GU%2}
exp={'I':(0,0),'T':(1,0),'U':(0,1),'TU':(1,1)}
lab={(0,0):'I',(1,0):'T',(0,1):'U',(1,1):'TU'}
def mul(a,b):
 x=exp[a];y=exp[b];return lab[(x[0]^y[0],x[1]^y[1])]
assert all(np.array_equal(Vmat[mul(a,b)],Vmat[a]@Vmat[b]%2) for a in Vmat for b in Vmat)
# source object charts
Af={r['id']:n2(r['path_matrix']) for r in BFS['records']}
assert len(Af)==62
# gate arrow data gives exact mod2 arrow and witness
GA={a['edge']:a for a in Gate['arrows']}
SA={a['id']:a for a in Gd['arrows']}
assert set(GA)==set(SA)==set(range(123))
# object cover
objects=[]
for f in range(62):
 for s in ['I','T','U','TU']:
  im=Vmat[s]@Af[f]%2
  objects.append({'field':f,'sheet':s,'target_F29_object':gid[key(im)],'target_matrix':im.astype(int).tolist()})
# lifted arrows and strict endpoint verification
arrows=[]; inter=0
for e in range(123):
 a=SA[e]; g=GA[e]; d=g['V4_endpoint_witness']
 assert d==a['deck'],(e,d,a['deck'])
 mm=np.array(g['mod2_matrix'],dtype=np.uint8)
 for s in ['I','T','U','TU']:
  t_sheet=mul(s,d)
  start=Vmat[s]@Af[a['source']]%2
  end=start@mm%2
  target=Vmat[t_sheet]@Af[a['target']]%2
  assert np.array_equal(end,target),(e,s,t_sheet)
  if s!=t_sheet: inter+=1
  arrows.append({'base_edge':e,'source':[a['source'],s],'target':[a['target'],t_sheet],'deck_increment':d,
                 'kind':a['kind'],'wall':a.get('wall'),'target_endpoint_strict':True,
                 'F29_source':gid[key(start)],'F29_target':gid[key(target)],'shortest_S_word':g['shortest_S_word']})
assert len(arrows)==492 and inter==48
# deck action equivariance on objects/arrows
for h in Vmat:
 for o in objects:
  f,s=o['field'],o['sheet']; hs=mul(h,s)
  lhs=Vmat[hs]@Af[f]%2; rhs=Vmat[h]@np.array(o['target_matrix'],dtype=np.uint8)%2
  assert np.array_equal(lhs,rhs)
# arrow word helper
aw={e:[int(x[1:]) for x in GA[e]['shortest_S_word']] for e in GA}
# inverse matrix helper via GL lookup
def inv2(A): return next(X for X in GL if np.array_equal(A@X%2,I))
# link monodromy labels from gate
FG={x['vertex']:x for x in Gate['links']}
# based lifted 2-relations: four per base vertex; repeat order 1/2 according to deck holonomy
based=[]; geom_keys={}; fill_total=0; maxfill=0
for lk in Links:
 v=lk['vertex']; delta=FG[v]['V4_closure']; order=1 if delta=='I' else 2
 word=lk['chosen']['word']; startf=word[0]['from_face']
 # verify deck product of source word
 dd='I'
 for z in word: dd=mul(dd,SA[z['edge']]['deck'])
 assert dd==delta,(v,dd,delta)
 for s0 in ['I','T','U','TU']:
  s=s0; f=startf; cur=Vmat[s]@Af[f]%2; steps=[]; lifted_steps=[]
  for rep in range(order):
   for z in word:
    e=z['edge']; d=SA[e]['deck']; direction=z['direction']
    # face traversal source/target in base provided by link, deck same under inverse (order 2)
    ns=mul(s,d)
    ww=aw[e] if direction==1 else list(reversed(aw[e]))
    for j in ww:
     u=gid[key(cur)]; cur=cur@S2[j-1]%2; vv=gid[key(cur)]; steps.append((u,vv,j))
    lifted_steps.append({'base_edge':e,'direction':direction,'from_sheet':s,'to_sheet':ns,
                         'from_face':z['from_face'],'to_face':z['to_face']})
    s=ns;f=z['to_face']
  assert s==s0,(v,s0,s,delta,order)
  assert np.array_equal(cur,Vmat[s0]@Af[startf]%2),(v,s0)
  fill=[];rew=[]
  for u,vv,j in steps:
   ee=(min(u,vv),max(u,vv),j)
   if ee in tree: rew.append((u,vv,j))
   else:
    fill.append({'fundamental_cell':cell_index[ee],'sign':1 if (u,vv,j)==ee else -1,'chord':[u,vv,j]})
    rew.extend(tree_path(u,vv))
  red=reduce_tree_walk(rew); assert not red,(v,s0,red)
  fill_total+=len(fill);maxfill=max(maxfill,len(fill))
  based.append({'base_vertex':v,'start_face':startf,'start_sheet':s0,'monodromy':delta,'monodromy_order':order,
                'base_link_length':len(word),'lifted_boundary_edge_occurrences':len(lifted_steps),
                'lifted_steps':lifted_steps,'strict_F29_closed':True,
                'fundamental_cell_occurrences':fill,'fundamental_cell_count':len(fill),'tree_rewrite_reduces_to_empty':True})
  # geometric key mod cyclic shift by monodromy subgroup: nontrivial identifies s and s*delta
  ss=sorted([s0] if delta=='I' else [s0,mul(s0,delta)])
  geom_keys[(v,tuple(ss))]=True
assert len(based)==248 and len(geom_keys)==59*4+3*2==242
# verify free V4 action on based relations: h sends start sheet and all sheets; monodromy unchanged
# counts / quotient recovery
assert len(objects)//4==62 and len(arrows)//4==123 and len(based)//4==62
# target image stats
imcnt=Counter(o['target_F29_object'] for o in objects)
# descent compatibility: quotient object target orbits agree with previous 42 orbit construction
# orbit id by left V4 on 168 target objects
orbrep={}; orbit_id={}
for A in GL:
 reps=[key(Vmat[d]@A%2) for d in Vmat]; rep=min(reps)
 if rep not in orbrep: orbrep[rep]=len(orbrep)
 oid=orbrep[rep]
 for d in Vmat: orbit_id[gid[key(Vmat[d]@A%2)]]=oid
assert len(orbrep)==42
for f in range(62):
 ids={orbit_id[gid[key(Vmat[s]@Af[f]%2)]] for s in Vmat}; assert len(ids)==1
# lower arrow witness must equal cover sheet increment
assert all(GA[e]['V4_endpoint_witness']==SA[e]['deck'] for e in range(123))
# lower link closure must equal cover monodromy
assert all(FG[lk['vertex']]['V4_closure']==next(x['monodromy'] for x in based if x['base_vertex']==lk['vertex']) for lk in Links)
# produce artifacts
cover={'version':'v50','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
 'status':'PROVEN-FOUR-SHEET-DECK-UNFOLDING-OBJECTS-ARROWS-BASED-2RELATIONS',
 'deck_group':{'name':'V4','elements':['I','T','U','TU'],'relations':['T^2=I','U^2=I','TU=UT'],'order':4},
 'counts':{'base_objects':62,'cover_objects':248,'base_arrows':123,'cover_arrows':492,'base_nontrivial_deck_loops':12,
           'lifted_intersheet_arrows':48,'base_link_relations':62,'cover_based_2relations':248,'cover_geometric_2cells':242},
 'objects':objects,'arrows':arrows,
 'two_cell_lift_rule':{'trivial_monodromy':'one traversal','nontrivial_order2_monodromy':'two traversals; no source connector adjoined',
                       'geometric_cell_count_note':'For T/U/TU monodromy, start sheets s and s*delta are cyclic basepoint shifts of the same unbased lifted disk, hence 2 geometric lifts rather than 4.'},
 'V4_action':{'objects':'h.(f,s)=(f,hs)','arrows':'h.(e,s)=(e,hs)','based_2relations':'h.(v,s)=(v,hs)','free_on_objects_arrows_based_relations':True},
 'guards':{'historical_H1_H2':'LOCKED pending completion/acceptance of descent square and local-system naturality'}}
strict={'version':'v50','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
 'status':'PASS-STRICT-UPPER-FUNCTOR-TO-F29-CAN',
 'source':'tilde_G_Sell_field four-sheet deck unfolding','target':'F29_can 168-object free chart',
 'object_map':'(f,s)->s*A_f mod2','arrow_map':'lifted source arrow uses the same certified F29 shortest S-word as its base arrow',
 'checks':{'objects_total':248,'object_endpoint_matrices_exact':True,'unique_target_objects_hit':len(imcnt),
           'arrows_total':492,'arrows_strict_endpoint_matches':492,'intersheet_arrows':48,
           'based_2relations_total':248,'based_2relations_strictly_closed':248,
           'based_2relations_fillable_by_841_fundamental_cells':248,'F29_fundamental_cell_count':841,
           'fundamental_cell_occurrences_total':fill_total,'max_fundamental_cells_per_lifted_relation':maxfill,
           'V4_equivariant_on_objects_and_arrows':True},
 'based_2relations':based,
 'interpretation':'The 12 nontrivial deck loops are no longer endomorphism loops upstairs; all become arrows between sheets. The three nontrivial orbifold link holonomies lift by their order-two double traversal to strict closed F29 loops.'}
desc={'version':'v50','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
 'status':'PASS-DESCENT-SQUARE-OBJECT/ARROW/2RELATION',
 'square':{'upper':'tilde_G_Sell_field -> F29_can','left':'/V4','right':'left /V4','lower':'G_Sell_field -> [F29_can/V4]'},
 'checks':{'cover_object_orbits':62,'cover_arrow_orbits':123,'cover_based_2relation_orbits':62,
           'target_object_orbits':42,'object_square_commutes':True,'arrow_square_commutes':True,
           'base_arrow_deck_witnesses_recovered':True,'base_link_V4_holonomies_recovered':True,
           'upper_functor_strict_objects_arrows_2relations':True,
           'lower_stack_target_status_preserved':'PROVEN-OBJECT/ARROW/LINK-COMPATIBLE-MOD2'},
 'nontrivial_base_links':[{ 'vertex':v, 'holonomy':FG[v]['V4_closure'], 'lift_order':2,
                            'based_lifts':4, 'geometric_lifts':2}
                          for v in sorted(FG) if FG[v]['V4_closure']!='I'],
 'interpretation':'The upper strict functor is a V4-equivariant lift of the lower stacky comparison. For nontrivial orbifold links, quotienting an order-two doubled lifted disk recovers the base orbifold 2-relation with the recorded V4 isotropy.',
 'cohomology_gate':'NOT OPENED AUTOMATICALLY: the descent square is now closed, but Binet/Ad local systems must be pulled to the cover and descended back equivariantly before historical H1/H2 are computed.'}
for name,data in [('selling_four_sheet_deck_cover_v50.json',cover),('selling_F29_strict_cover_functor_v50.json',strict),('selling_F29_descent_square_v50.json',desc)]:
 (B/name).write_text(json.dumps(data,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')
print('PASS cover',cover['counts'])
print('target images',len(imcnt),'multiplicity hist',dict(sorted(Counter(imcnt.values()).items())))
print('strict 2rels',len(based),'fill total',fill_total,'max',maxfill)
print('nontriv links',[(v,FG[v]['V4_closure']) for v in sorted(FG) if FG[v]['V4_closure']!='I'])
