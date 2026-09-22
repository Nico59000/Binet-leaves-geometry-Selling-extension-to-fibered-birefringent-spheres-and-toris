import sys,itertools,json,time
import sympy as sp
sys.path.insert(0,'/mnt/data/v50_work')
import selling_reconstruction_engine_v50 as e
EDGE={'g':(2,3),'h':(1,3),'k':(1,2),'l':(0,1),'m':(0,2),'n':(0,3)}
V=[sp.Matrix([-1,-1,-1]),sp.eye(3)[:,0],sp.eye(3)[:,1],sp.eye(3)[:,2]]

def pairmap(u,v):
 eu=set(EDGE[u]);ev=set(EDGE[v]);inter=eu&ev
 if len(inter)!=1:return None
 r=next(iter(inter));su=next(iter(eu-{r}));sv=next(iter(ev-{r}));rem=next(iter(set(range(4))-{r,su,sv}))
 pi=[rem,r,sv,su];A=sp.Matrix.hstack(V[pi[1]],V[pi[2]],V[pi[3]])
 return A
MAP={(u,v):pairmap(u,v) for u,v in itertools.combinations(EDGE,2) if pairmap(u,v) is not None}

def crosspairs(M):
 out=[]
 for (u,v),A in MAP.items():
  if e.standard_hk_crossing_feasible(A.T*M*A):out.append((u,v))
 return out

def hascross(M):return bool(crosspairs(M))
# start at generic known feasible M0
Qstart=sp.diag(-1,12,-5);M0=e.S[1].T*Qstart*e.S[1]
seen={e.numerical_form_key(M0):M0}; front=[M0]; edges=[]; rejected={}
for depth in range(40):
 print('depth',depth,'front',len(front),'seen',len(seen))
 nxt=[]
 for M in front:
  km=e.numerical_form_key(M)
  for wall,(i,j) in e.WALL_TO_MOVES.items():
   A=e.S[i-1];N=A.T*M*A;kn=e.numerical_form_key(N)
   if hascross(N):
    edges.append((km,wall,kn))
    if kn not in seen:seen[kn]=N;nxt.append(N)
   else:rejected.setdefault(kn,N)
 front=nxt
 if not front:break
print('FINAL',len(seen),'edges',len(edges),'rejected_unique',len(rejected))
# check candidates that may be fissure-only later
out={'count':len(seen),'forms':[],'edges':[],'rejected_count':len(rejected),'rejected_forms':[]}
idx={k:i for i,k in enumerate(sorted(seen,key=str))}
for k in sorted(seen,key=str):
 M=seen[k];out['forms'].append({'id':idx[k],'six':list(map(int,e.six_from_gram(M))),'crosspairs':[list(p) for p in crosspairs(M)]})
for a,w,b in edges:out['edges'].append([idx[a],w,idx[b]])
for k in sorted(rejected,key=str):
 M=rejected[k];out['rejected_forms'].append({'six':list(map(int,e.six_from_gram(M)))})
json.dump(out,open('/mnt/data/v50_work/forms_crossing_enum_v50_tmp.json','w'),indent=2)
