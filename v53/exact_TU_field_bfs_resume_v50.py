#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp,sys,json,time
from pathlib import Path
HERE=Path('/mnt/data/v50_work');sys.path.insert(0,str(HERE))
from exact_boundary_arc_classifier_v50 import classify_cell
from selling_global_wall_linearization_v50 import R,GT,GU
from selling_reconstruction_engine_v49 import S,RELABEL
first={'g':9,'h':2,'k':1,'l':3,'m':5,'n':4}
def tup(A): return tuple(int(x) for x in list(A))
def mat(x): return sp.Matrix([[sp.Integer(v) for v in row] for row in x])
def powers(M,N):
 d={0:sp.eye(3)}
 for i in range(1,N+1):d[i]=d[i-1]*M
 Mi=M.inv()
 for i in range(1,N+1):d[-i]=d[-i+1]*Mi
 return d

def load_state():
 p=HERE/'exact_TU_field_bfs_checkpoint_v50.json';d=json.load(open(p));n=d['field_count'];reps=[None]*n;Ls=[None]*n
 reps[0]=sp.eye(3);Ls[0]=sp.eye(6)
 for rec in d['records']:
  i=rec['id']
  if reps[i] is None: reps[i]=mat(rec['path_matrix'])
  if Ls[i] is None: Ls[i]=mat(rec['wall_action_matrix'])
  for w,e in rec['neighbors'].items():
   if e.get('new'):
    j=e['target'];reps[j]=reps[i]*S[first[w]-1];Ls[j]=R[w]*Ls[i]
 assert all(x is not None for x in reps)
 return d,reps,Ls

def main(batch=12):
 d,reps,Ls=load_state();bound=d['exp_bound'];PT=powers(GT,bound);PU=powers(GU,bound);deck={(a,b):PT[a]*PU[b] for a in PT for b in PU};variants={}
 def addv(idx,A):
  for (a,b),D in deck.items():
   DA=D*A
   for bi,B in enumerate(RELABEL):
    X=DA*B;variants.setdefault(tup(X),(idx,a,b,bi,1));variants.setdefault(tup(-X),(idx,a,b,bi,-1))
 for i,A in enumerate(reps):addv(i,A)
 qi=d['processed'];records=d['records'];steps=0
 while qi<len(reps) and steps<batch:
  t=time.time();A,L=reps[qi],Ls[qi];cls=classify_cell(L);rec={'id':qi,'path_matrix':A.tolist(),'wall_action_matrix':L.tolist(),'ordinary_walls':cls['ordinary_walls'],'multiwall_components':[{'factor':x['factor'],'walls':x['walls']} for x in cls['multiwall_components']], 'neighbors':{}}
  for w in cls['ordinary_walls']:
   An=A*S[first[w]-1];Ln=R[w]*L;key=tup(An)
   if key in variants:
    j,a,b,bi,sgn=variants[key];rec['neighbors'][w]={'target':j,'new':False,'deck_exponents':[a,b],'relabel_index':bi,'sign':sgn}
   else:
    j=len(reps);reps.append(An);Ls.append(Ln);addv(j,An);rec['neighbors'][w]={'target':j,'new':True}
  records.append(rec);qi+=1;steps+=1
  state={'version':'v50','status':'RUNNING','exp_bound':bound,'field_count':len(reps),'processed':qi,'records':records}
  (HERE/'exact_TU_field_bfs_checkpoint_v50.json').write_text(json.dumps(state,indent=2,default=str),encoding='utf-8')
  print('field',qi-1,'processed',qi,'known',len(reps),'walls',cls['ordinary_walls'],'multi',[(x['walls'],x['factor']) for x in cls['multiwall_components']],'sec',round(time.time()-t,2),flush=True)
 closed=qi==len(reps)
 if closed:
  out={'version':'v50','status':'PASS-CLOSED','exp_bound':bound,'field_count':len(reps),'processed':qi,'closed':True,'records':records}
  (HERE/'exact_TU_field_bfs_v50.json').write_text(json.dumps(out,indent=2,default=str),encoding='utf-8')
 print(json.dumps({'processed':qi,'known':len(reps),'closed':closed}))
if __name__=='__main__':main()
