#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp, sys, json, hashlib, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from exact_boundary_arc_classifier_v50 import classify_cell
from selling_global_wall_linearization_v50 import R, GT, GU
from selling_reconstruction_engine_v49 import S, RELABEL
first={'g':9,'h':2,'k':1,'l':3,'m':5,'n':4}
def tup(A):return tuple(int(x) for x in list(A))
def powers(M,N):
 d={0:sp.eye(3)}
 for i in range(1,N+1):d[i]=d[i-1]*M
 Mi=M.inv()
 for i in range(1,N+1):d[-i]=d[-i+1]*Mi
 return d

def run(exp_bound=3,max_fields=120):
 PT=powers(GT,exp_bound);PU=powers(GU,exp_bound);deck={(a,b):PT[a]*PU[b] for a in PT for b in PU}
 reps=[sp.eye(3)]; Lreps=[sp.eye(6)];records=[];variants={}
 def add_variants(idx,A):
  for (a,b),D in deck.items():
   DA=D*A
   for bi,B in enumerate(RELABEL):
    X=DA*B
    variants.setdefault(tup(X),(idx,a,b,bi,1));variants.setdefault(tup(-X),(idx,a,b,bi,-1))
 add_variants(0,reps[0])
 qi=0
 ck=HERE/'exact_TU_field_bfs_checkpoint_v50.json'
 while qi<len(reps) and len(reps)<=max_fields:
  t0=time.time();A,L=reps[qi],Lreps[qi];cls=classify_cell(L)
  rec={'id':qi,'path_matrix':A.tolist(),'wall_action_matrix':L.tolist(),'ordinary_walls':cls['ordinary_walls'],
       'multiwall_components':[{'factor':x['factor'],'walls':x['walls']} for x in cls['multiwall_components']], 'neighbors':{}}
  for w in cls['ordinary_walls']:
   An=A*S[first[w]-1]; Ln=R[w]*L;key=tup(An)
   if key in variants:
    j,a,b,bi,sgn=variants[key];rec['neighbors'][w]={'target':j,'new':False,'deck_exponents':[a,b],'relabel_index':bi,'sign':sgn}
   else:
    j=len(reps);reps.append(An);Lreps.append(Ln);add_variants(j,An);rec['neighbors'][w]={'target':j,'new':True}
  records.append(rec);qi+=1
  state={'version':'v50','status':'RUNNING','exp_bound':exp_bound,'field_count':len(reps),'processed':qi,'records':records}
  ck.write_text(json.dumps(state,indent=2,default=str),encoding='utf-8')
  print('field',qi-1,'processed',qi,'known',len(reps),'walls',cls['ordinary_walls'],'multi',[(x['walls'],x['factor']) for x in cls['multiwall_components']],'sec',round(time.time()-t0,2),flush=True)
 closed=(qi==len(reps) and len(reps)<=max_fields)
 out={'version':'v50','status':'PASS-CLOSED' if closed else 'OPEN/TRUNCATED','exp_bound':exp_bound,'field_count':len(reps),'processed':qi,'closed':closed,'records':records}
 p=HERE/'exact_TU_field_bfs_v50.json';p.write_text(json.dumps(out,indent=2,default=str),encoding='utf-8')
 return out
if __name__=='__main__':
 out=run();print(json.dumps({'status':out['status'],'field_count':out['field_count'],'processed':out['processed'],'closed':out['closed']},indent=2))
