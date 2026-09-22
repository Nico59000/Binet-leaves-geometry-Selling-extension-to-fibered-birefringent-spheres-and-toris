#!/usr/bin/env python3
import sys,itertools,json,time,multiprocessing as mp,signal
from pathlib import Path
import sympy as sp
sys.path.insert(0,'/mnt/data/v50_work')
import selling_reconstruction_engine_v50 as e
import selling_fissure_rur_v50 as fr
EDGE={'g':(2,3),'h':(1,3),'k':(1,2),'l':(0,1),'m':(0,2),'n':(0,3)}
V=[sp.Matrix([-1,-1,-1]),sp.eye(3)[:,0],sp.eye(3)[:,1],sp.eye(3)[:,2]]

def pairmap(u,v):
 eu=set(EDGE[u]);ev=set(EDGE[v]);inter=eu&ev
 if len(inter)!=1:return None
 r=next(iter(inter));su=next(iter(eu-{r}));sv=next(iter(ev-{r}));rem=next(iter(set(range(4))-{r,su,sv}))
 pi=[rem,r,sv,su];return sp.Matrix.hstack(V[pi[1]],V[pi[2]],V[pi[3]])
MAP={(u,v):pairmap(u,v) for u,v in itertools.combinations(EDGE,2) if pairmap(u,v) is not None}

def crosspairs(M):
 return [(u,v) for (u,v),A in MAP.items() if e.standard_hk_crossing_feasible(A.T*M*A)]

def canon_rep(M):
 k=e.numerical_form_key(M)
 # materialize actual canonical representative, not merely key, for stable serialization
 best=None; bestk=None
 for A in e.RELABEL:
  N=A.T*M*A; nk=tuple(int(q) for q in list(N))
  sk=tuple(map(str,nk))
  if bestk is None or sk<bestk: bestk=sk;best=N
 return k,best

def fissure_worker(six):
 def _alarm(signum,frame): raise TimeoutError('per-form fissure timeout')
 signal.signal(signal.SIGALRM,_alarm); signal.alarm(8)
 try:
  M=e.gram_from_six(*map(sp.Integer,six))
  r=fr.classify_form(M,False)
  return {'six':list(six),'status':r['status'],'events':r['reduced_fissure_events'],
          'unresolved':r['unresolved'],
          'pair_summary':[(p['pair'],p.get('status'),p.get('reduced_event_count',0)) for p in r['pairs']]}
 except TimeoutError as exc:
  return {'six':list(six),'status':'PARTIAL_FAIL_CLOSED','events':0,
          'unresolved':['TIMEOUT_8S'],'pair_summary':[]}
 finally:
  signal.alarm(0)


def main():
 Qstart=sp.diag(-1,12,-5);M0=e.S[1].T*Qstart*e.S[1]
 k0,M0=canon_rep(M0)
 seen={k0:M0}; queue=[M0]; class_cache={k0:{'crosspairs':crosspairs(M0),'fissure_events':None,'status':'CROSSING'}}
 edge_rows=[]; rejected_final={}; unresolved={}; wave=0;t0=time.time()
 while queue:
  print('WAVE',wave,'front',len(queue),'seen',len(seen),flush=True)
  candidates={}
  for M in queue:
   ks=e.numerical_form_key(M)
   for wall,(i,j) in e.WALL_TO_MOVES.items():
    A=e.S[i-1];N=A.T*M*A;kn,Nc=canon_rep(N)
    edge_rows.append((ks,wall,kn))
    if kn in seen or kn in candidates or kn in rejected_final or kn in unresolved:continue
    cp=crosspairs(Nc)
    if cp:
     candidates[kn]={'M':Nc,'crosspairs':cp,'kind':'CROSSING'}
    else:
     candidates[kn]={'M':Nc,'crosspairs':[],'kind':'NEEDS_FISSURE'}
  need=[(k,r) for k,r in candidates.items() if r['kind']=='NEEDS_FISSURE']
  print(' candidates',len(candidates),'need_fissure',len(need),flush=True)
  fres={}
  if need:
   args=[tuple(map(int,e.six_from_gram(r['M']))) for _,r in need]
   with mp.Pool(processes=min(8,len(args))) as pool:
    vals=pool.map(fissure_worker,args)
   for (k,r),v in zip(need,vals):fres[k]=v
  nxt=[]; added_cross=added_fiss=rej=unk=0
  for k,r in candidates.items():
   if r['kind']=='CROSSING':
    seen[k]=r['M'];nxt.append(r['M']);class_cache[k]={'crosspairs':r['crosspairs'],'fissure_events':None,'status':'CROSSING'};added_cross+=1
   else:
    v=fres[k]
    if v['events']>0:
     seen[k]=r['M'];nxt.append(r['M']);class_cache[k]={'crosspairs':[],'fissure_events':v['events'],'status':'FISSURE_ONLY','fissure_detail':v};added_fiss+=1
    elif v['status']=='PASS':
     rejected_final[k]={'M':r['M'],'fissure_detail':v};rej+=1
    else:
     unresolved[k]={'M':r['M'],'fissure_detail':v};unk+=1
  print(' add_cross',added_cross,'add_fiss',added_fiss,'reject',rej,'unresolved',unk,'elapsed',round(time.time()-t0,1),flush=True)
  chk={'wave_completed':wave,'seen':len(seen),'next_frontier':len(nxt),'rejected':len(rejected_final),'unresolved':len(unresolved),'elapsed':time.time()-t0}
  Path('/mnt/data/v50_work/selling_reduced_census_wave_checkpoint_v50.json').write_text(json.dumps(chk,indent=2))
  queue=nxt;wave+=1
  if wave>100:raise RuntimeError('wave limit')
 # Dedup actual adjacency from all seen now, recompute all six walls and target membership.
 keys=sorted(seen,key=str); idx={k:i for i,k in enumerate(keys)}
 adj=[]; outside=[]
 for k in keys:
  M=seen[k]
  for wall,(i,j) in e.WALL_TO_MOVES.items():
   N=e.S[i-1].T*M*e.S[i-1];kn=e.numerical_form_key(N)
   rec=[idx[k],wall,idx[kn] if kn in idx else None]
   if kn in idx:adj.append(rec)
   else:outside.append([idx[k],wall, list(map(int,e.six_from_gram(N))), 'UNRESOLVED' if kn in unresolved else 'REJECTED'])
 forms=[]
 for k in keys:
  M=seen[k];cc=class_cache[k]
  forms.append({'id':idx[k],'six':list(map(int,e.six_from_gram(M))),'crosspairs':[list(p) for p in cc['crosspairs']],
                'admission':cc['status'],'reduced_fissure_events':cc.get('fissure_events')})
 out={'version':'v50','status':'PASS_CLOSED_REDUCED_FORM_CENSUS' if not unresolved else 'PARTIAL_FAIL_CLOSED',
      'waves':wave,'reduced_numerical_forms':len(seen),
      'crossing_admitted':sum(1 for c in class_cache.values() if c['status']=='CROSSING'),
      'fissure_only_admitted':sum(1 for c in class_cache.values() if c['status']=='FISSURE_ONLY'),
      'rejected_certified':len(rejected_final),'unresolved_candidates':len(unresolved),
      'adjacency_rows':len(adj),'outside_rows':len(outside),'forms':forms,'adjacency':adj,'outside':outside,
      'unresolved':[{'six':list(map(int,e.six_from_gram(r['M']))),'detail':r['fissure_detail']} for r in unresolved.values()],
      'elapsed_seconds':time.time()-t0,
      'guard':'Numerical forms are not yet historical C2 fields; connected boundary components and event order remain to be constructed.'}
 Path('/mnt/data/v50_work/selling_reduced_numerical_form_census_v50.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
 print(json.dumps({k:out[k] for k in ['status','waves','reduced_numerical_forms','crossing_admitted','fissure_only_admitted','rejected_certified','unresolved_candidates','adjacency_rows','outside_rows','elapsed_seconds']},indent=2))
if __name__=='__main__':main()
