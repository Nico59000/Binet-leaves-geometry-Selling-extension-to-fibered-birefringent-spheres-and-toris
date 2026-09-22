#!/usr/bin/env python3
import json,sys,math
from pathlib import Path
import sympy as sp
HERE=Path('/mnt/data/v50_work');sys.path.insert(0,str(HERE))
from selling_global_wall_linearization_v50 import GT,GU
BFS=json.load(open(HERE/'exact_TU_field_bfs_v50.json'));recs={r['id']:r for r in BFS['records']}
F=json.load(open(HERE/'selling_face_cycles_v50.json'));faces={x['field']:x for x in F['faces']}

def phi(D,p):
 q,r=p;v=sp.Matrix([1,-sp.Float(r,40),sp.Float(q,40)]);u=D.T*v
 return (float(sp.N(u[2]/u[0],30)),float(sp.N(-u[1]/u[0],30)))
def edist(E1,E2):
 (a,b),(c,d)=E1,E2
 ds=lambda x,y:(x[0]-y[0])**2+(x[1]-y[1])**2
 return min(max(ds(a,c),ds(b,d)),max(ds(a,d),ds(b,c)))
def cycle_darts(fid,ci=0):
 f=faces[fid];cy=f['cycles'][ci];out=[]
 for si,s in enumerate(cy['segments']):
  co=f['components'][s['component']];pts=[(f['vertices'][v]['q'],f['vertices'][v]['r']) for v in s['v']]
  out.append({'key':(fid,si),'field':fid,'cycle':ci,'si':si,'segment':s,'component':s['component'],'kind':co['kind'],'walls':co['walls'],'pts':pts})
 return out
# lobe B->A canonicalization under TU
lobeA=cycle_darts(55,0);lobeB=cycle_darts(55,1)
mapBtoA={}
for db in lobeB:
 E=[phi(GT*GU,p) for p in db['pts']]
 dist,da=min((edist(E,da['pts']),da) for da in lobeA)
 assert dist<1e-10,(db,dist)
 mapBtoA[db['si']]=da['si']
# canonical face-edge occurrences: cycle0 only
all_darts=[];bykey={}
for i in range(62):
 for d in cycle_darts(i,0):all_darts.append(d);bykey[d['key']]=d
ordinary=[d for d in all_darts if d['kind']=='ORDINARY']
# DSU
parent={d['key']:d['key'] for d in all_darts}
def find(x):
 while parent[x]!=x:
  parent[x]=parent[parent[x]];x=parent[x]
 return x
def union(a,b):
 a,b=find(a),find(b)
 if a!=b: parent[max(a,b)]=min(a,b)
# pair ordinary each to nearest target occurrence under exact deck metadata
links=[];fails=[]
for d in ordinary:
 i=d['field'];w=d['walls'][0];meta=recs[i]['neighbors'][w];j=meta['target']
 if meta.get('new'): D=sp.eye(3)
 else:
  a,b=meta['deck_exponents'];D=GT**a*GU**b
 Et=[phi(D,p) for p in d['pts']]
 cand=[]
 for ci in ([0,1] if j==55 else [0]):
  for td in cycle_darts(j,ci):
   if td['kind']!='ORDINARY':continue
   cand.append((edist(Et,td['pts']),ci,td))
 cand.sort(key=lambda x:x[0])
 if not cand or cand[0][0]>1e-8:
  fails.append((d['key'],j,cand[:3]));continue
 dist,ci,td=cand[0]
 tsi=td['si'] if not(j==55 and ci==1) else mapBtoA[td['si']]
 tk=(j,tsi);union(d['key'],tk);links.append({'from':list(d['key']),'to':list(tk),'distance2':dist,'raw_target_cycle':ci,'raw_target_si':td['si'],'wall':w,'meta':meta})
assert not fails,fails
# Direct deck orbit identifications within each same face, including multiwall edges.
# Test D in V4 on endpoints; if exact-near match canonical cycle0 or 55 lobeB canonicalized, union.
decks=[('I',sp.eye(3)),('T',GT),('U',GU),('TU',GT*GU)]
deck_links=[]
for i in range(62):
 srcs=cycle_darts(i,0)
 targets0=cycle_darts(i,0)
 targets1=cycle_darts(i,1) if i==55 else []
 for d in srcs:
  for nm,D in decks[1:]:
   Et=[phi(D,p) for p in d['pts']]
   cand=[]
   for ci,arr in [(0,targets0),(1,targets1)]:
    for td in arr:
     if td['kind']!=d['kind']: continue
     cand.append((edist(Et,td['pts']),ci,td))
   if not cand:continue
   dist,ci,td=min(cand,key=lambda x:x[0])
   if dist<1e-10:
    tsi=td['si'] if not(i==55 and ci==1) else mapBtoA[td['si']]
    tk=(i,tsi);union(d['key'],tk);deck_links.append({'deck':nm,'from':list(d['key']),'to':list(tk),'distance2':dist,'raw_cycle':ci,'raw_si':td['si']})
# collect orbits
orb={}
for k in parent:orb.setdefault(find(k),[]).append(k)
orbits=sorted([sorted(v) for v in orb.values()],key=lambda x:x[0])
# stats
from collections import Counter
hist=Counter(len(x) for x in orbits)
ord_orb=sum(1 for O in orbits if any(bykey[k]['kind']=='ORDINARY' for k in O))
mix=[O for O in orbits if len({bykey[k]['kind'] for k in O})>1]
print('all face-edge darts',len(all_darts),'ordinary',len(ordinary),'multi',len(all_darts)-len(ordinary))
print('edge orbits',len(orbits),'size hist',dict(hist),'ordinary-containing',ord_orb,'mixed',len(mix))
for n in sorted(hist):
 if n>=3: print(' orbit size',n,[O for O in orbits if len(O)==n][:8])
out={'version':'v50','status':'PASS','face_edge_occurrences':len(all_darts),'ordinary_occurrences':len(ordinary),'multiwall_occurrences':len(all_darts)-len(ordinary),'edge_orbit_count':len(orbits),'orbit_size_histogram':dict(hist),'edge_orbits':orbits,'ordinary_transport_links':links,'direct_deck_links':deck_links,'lobeB_to_lobeA_under_TU':mapBtoA,'mixed_kind_orbits':mix}
(HERE/'selling_edge_orbits_v50.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
