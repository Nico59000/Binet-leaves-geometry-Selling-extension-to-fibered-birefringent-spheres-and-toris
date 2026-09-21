#!/usr/bin/env python3
import json,sys,math
from pathlib import Path
import sympy as sp
HERE=Path('/mnt/data/v50_work');sys.path.insert(0,str(HERE))
from selling_global_wall_linearization_v50 import GT,GU
BFS=json.load(open(HERE/'exact_TU_field_bfs_v50.json'));recs={r['id']:r for r in BFS['records']}
F=json.load(open(HERE/'selling_face_cycles_v50.json'));faces={x['field']:x for x in F['faces']}

def deck_for(meta):
 if meta.get('new'): return sp.eye(3)
 a,b=meta.get('deck_exponents',[0,0]);return (GT**a)*(GU**b)
def phi(D,p):
 q,r=p;v=sp.Matrix([1,-sp.Float(r,30),sp.Float(q,30)]);u=D.T*v
 return (float(sp.N(u[2]/u[0],20)),float(sp.N(-u[1]/u[0],20)))
def edist(E1,E2):
 # unordered endpoint set squared max best matching
 (a,b),(c,d)=E1,E2
 def ds(x,y):return (x[0]-y[0])**2+(x[1]-y[1])**2
 return min(max(ds(a,c),ds(b,d)),max(ds(a,d),ds(b,c)))
def cycle_darts(fid,cycle_index=0):
 f=faces[fid];cy=f['cycles'][cycle_index];out=[]
 for si,s in enumerate(cy['segments']):
  co=f['components'][s['component']]
  pts=[(f['vertices'][v]['q'],f['vertices'][v]['r']) for v in s['v']]
  out.append({'field':fid,'cycle':cycle_index,'si':si,'segment':s,'component':s['component'],'kind':co['kind'],'walls':co['walls'],'pts':pts})
 return out
# Map lobe B darts of 55 to lobe A darts under U deck coordinate action.
lobeA=cycle_darts(55,0);lobeB=cycle_darts(55,1)
mapBtoA={}
for db in lobeB:
 E=[phi(GU,p) for p in db['pts']]
 cand=[(edist(E,da['pts']),da) for da in lobeA]
 dist,da=min(cand,key=lambda x:x[0])
 if dist>1e-10:raise RuntimeError(('U lobe map failed',db,dist))
 mapBtoA[db['si']]=da['si']
print('U lobeB->A',mapBtoA)
# Canonical quotient darts: all cycles0, including lobe A only for 55.
darts=[];bykey={}
for i in range(62):
 for d in cycle_darts(i,0):
  if d['kind']=='ORDINARY':
   key=(i,d['si']);d['key']=key;darts.append(d);bykey[key]=d
print('ordinary darts',len(darts))
pairs={};fails=[]
for d in darts:
 i=d['field'];w=d['walls'][0];meta=recs[i]['neighbors'][w];j=meta['target'];D=deck_for(meta)
 Et=[phi(D,p) for p in d['pts']]
 # candidates in target: normal cycle0; for 55 both lobes then quotient B->A
 cand=[]
 for tj in ([0,1] if j==55 else [0]):
  for td in cycle_darts(j,tj):
   if td['kind']!='ORDINARY':continue
   dist=edist(Et,td['pts'])
   cand.append((dist,tj,td))
 cand.sort(key=lambda x:x[0])
 if not cand or cand[0][0]>1e-8:
  fails.append({'dart':d['key'],'wall':w,'target':j,'best':cand[:3]});continue
 dist,tj,td=cand[0]
 tsi=td['si'] if not (j==55 and tj==1) else mapBtoA[td['si']]
 tkey=(j,tsi)
 pairs[d['key']]={'target_dart':tkey,'distance2':dist,'target_raw_cycle':tj,'target_raw_si':td['si'],'deck':D.tolist(),'meta':meta}
print('fails',len(fails))
if fails:
 for f in fails[:20]:print(f)
# involution check after quotient
bad=[]
for k,p in pairs.items():
 kk=tuple(p['target_dart']);qv=pairs.get(kk)
 if not qv or tuple(qv['target_dart'])!=k:bad.append((k,kk,qv['target_dart'] if qv else None))
print('involution bad',len(bad),bad[:20])
# edge orbits
seen=set();edges=[]
for k in pairs:
 if k in seen:continue
 kk=tuple(pairs[k]['target_dart']);orb={k,kk};seen|=orb;edges.append(sorted(orb))
print('ordinary edge orbits',len(edges),'fixed darts',sum(1 for e in edges if len(e)==1))
out={'version':'v50','status':'PASS' if not fails and not bad else 'FAIL','ordinary_darts':len(darts),'ordinary_edge_orbits':len(edges),'fixed_ordinary_edge_orbits':sum(1 for e in edges if len(e)==1),'lobeB_to_lobeA_under_U':mapBtoA,'pairing':{str(k):v for k,v in pairs.items()},'edge_orbits':edges,'failures':fails,'involution_failures':bad}
(HERE/'selling_quotient_edge_pairing_v50.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
