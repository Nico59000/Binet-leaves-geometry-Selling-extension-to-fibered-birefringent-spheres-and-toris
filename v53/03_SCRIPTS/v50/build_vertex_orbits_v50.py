#!/usr/bin/env python3
import json,sys,math
from pathlib import Path
import sympy as sp
HERE=Path('/mnt/data/v50_work');sys.path.insert(0,str(HERE))
from selling_global_wall_linearization_v50 import GT,GU
BFS=json.load(open(HERE/'exact_TU_field_bfs_v50.json'));recs={r['id']:r for r in BFS['records']}
F=json.load(open(HERE/'selling_face_cycles_v50.json'));faces={x['field']:x for x in F['faces']}
E=json.load(open(HERE/'selling_edge_orbits_v50.json'))

def phi(D,p):
 q,r=p;v=sp.Matrix([1,-sp.Float(r,40),sp.Float(q,40)]);u=D.T*v
 return (float(sp.N(u[2]/u[0],30)),float(sp.N(-u[1]/u[0],30)))
def ds(a,b):return (a[0]-b[0])**2+(a[1]-b[1])**2

def cycle_info(fid,ci=0):
 f=faces[fid];cy=f['cycles'][ci];V=cy['vertices'];out=[]
 for k,s in enumerate(cy['segments']):
  co=f['components'][s['component']]
  v0=V[k];v1=V[(k+1)%len(V)]
  p0=(f['vertices'][v0]['q'],f['vertices'][v0]['r']);p1=(f['vertices'][v1]['q'],f['vertices'][v1]['r'])
  out.append({'field':fid,'cycle':ci,'si':k,'kind':co['kind'],'walls':co['walls'],'v0':v0,'v1':v1,'p0':p0,'p1':p1})
 return out
# canonical vertex occurrences cycle0
occ=[]
for i in range(62):
 for v in faces[i]['cycles'][0]['vertices']:
  k=(i,v)
  if k not in occ:occ.append(k)
parent={k:k for k in occ}
def find(x):
 while parent[x]!=x:
  parent[x]=parent[parent[x]];x=parent[x]
 return x
def union(a,b):
 if a not in parent or b not in parent:raise KeyError((a,b))
 a,b=find(a),find(b)
 if a!=b:parent[max(a,b)]=min(a,b)
# map lobeB vertex -> lobeA vertex under TU
Averts=faces[55]['cycles'][0]['vertices'];Bverts=faces[55]['cycles'][1]['vertices']
mapBv={}
for vb in Bverts:
 p=(faces[55]['vertices'][vb]['q'],faces[55]['vertices'][vb]['r']);pp=phi(GT*GU,p)
 va=min(Averts,key=lambda x:ds(pp,(faces[55]['vertices'][x]['q'],faces[55]['vertices'][x]['r'])))
 dd=ds(pp,(faces[55]['vertices'][va]['q'],faces[55]['vertices'][va]['r']))
 assert dd<1e-10,(vb,va,dd)
 mapBv[vb]=va
print('55 Bvertex->A',mapBv)
# helper target dart raw info and canonicalize vertex if lobeB
def raw_dart(fid,ci,si):return cycle_info(fid,ci)[si]
def canon_target_vertex(fid,ci,v):
 if fid==55 and ci==1:return mapBv[v]
 return v
# ordinary links recompute endpoint mapping
endpoint_links=[];fails=[]
for i in range(62):
 for d in cycle_info(i,0):
  if d['kind']!='ORDINARY':continue
  w=d['walls'][0];meta=recs[i]['neighbors'][w];j=meta['target']
  if meta.get('new'):D=sp.eye(3)
  else:
   a,b=meta['deck_exponents'];D=GT**a*GU**b
  P=[phi(D,d['p0']),phi(D,d['p1'])]
  # target nearest ordinary dart, allowing both lobes for 55
  cand=[]
  for ci in ([0,1] if j==55 else [0]):
   for td in cycle_info(j,ci):
    if td['kind']!='ORDINARY':continue
    # best endpoint matching
    same=max(ds(P[0],td['p0']),ds(P[1],td['p1']));rev=max(ds(P[0],td['p1']),ds(P[1],td['p0']))
    cand.append((min(same,rev),same<=rev,ci,td))
  dist,same,ci,td=min(cand,key=lambda x:x[0])
  if dist>1e-8:fails.append((i,d['si'],j,dist));continue
  tv0=canon_target_vertex(j,ci,td['v0']);tv1=canon_target_vertex(j,ci,td['v1'])
  if same:
   pairs=[((i,d['v0']),(j,tv0)),((i,d['v1']),(j,tv1))]
  else:pairs=[((i,d['v0']),(j,tv1)),((i,d['v1']),(j,tv0))]
  for a,b in pairs:union(a,b)
  endpoint_links.append({'from':[i,d['si']],'to_raw':[j,ci,td['si']],'same_orientation':same,'distance2':dist,'vertex_pairs':[[list(a),list(b)] for a,b in pairs]})
assert not fails,fails
# direct deck identifications, including multiwalls
for i in range(62):
 srcs=cycle_info(i,0)
 for d in srcs:
  for nm,D in [('T',GT),('U',GU),('TU',GT*GU)]:
   P=[phi(D,d['p0']),phi(D,d['p1'])];cand=[]
   for ci in ([0,1] if i==55 else [0]):
    for td in cycle_info(i,ci):
     if td['kind']!=d['kind']:continue
     same=max(ds(P[0],td['p0']),ds(P[1],td['p1']));rev=max(ds(P[0],td['p1']),ds(P[1],td['p0']))
     cand.append((min(same,rev),same<=rev,ci,td))
   dist,same,ci,td=min(cand,key=lambda x:x[0])
   if dist<1e-10:
    tv0=canon_target_vertex(i,ci,td['v0']);tv1=canon_target_vertex(i,ci,td['v1'])
    if same:pairs=[((i,d['v0']),(i,tv0)),((i,d['v1']),(i,tv1))]
    else:pairs=[((i,d['v0']),(i,tv1)),((i,d['v1']),(i,tv0))]
    for a,b in pairs:union(a,b)
# collect orbits
orbd={}
for k in parent:orbd.setdefault(find(k),[]).append(k)
orbits=sorted([sorted(v) for v in orbd.values()],key=lambda x:x[0])
idx={k:i for i,O in enumerate(orbits) for k in O}
from collections import Counter
hist=Counter(len(O) for O in orbits)
print('vertex occurrences',len(occ),'orbits',len(orbits),'hist',dict(hist))
# Edge endpoints at quotient using each edge orbit representative
edge_end=[];edge_loops=0
for ei,O in enumerate(E['edge_orbits']):
 k=tuple(O[0]);d=cycle_info(k[0],0)[k[1]];a=idx[(k[0],d['v0'])];b=idx[(k[0],d['v1'])]
 if a==b:edge_loops+=1
 edge_end.append({'edge':ei,'v0':a,'v1':b,'loop':a==b,'representative':list(k)})
print('edge loops',edge_loops)
out={'version':'v50','status':'PASS','vertex_occurrences':len(occ),'vertex_orbit_count':len(orbits),'orbit_size_histogram':dict(hist),'vertex_orbits':orbits,'field55_lobeB_vertex_to_A_under_TU':mapBv,'ordinary_endpoint_links':endpoint_links,'edge_endpoints':edge_end,'edge_loop_count':edge_loops}
(HERE/'selling_vertex_orbits_v50.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
