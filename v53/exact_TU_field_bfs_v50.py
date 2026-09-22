#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp, sys, json, hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
# importing classifier also imports global linearization; variables are exact
from exact_boundary_arc_classifier_v50 import classify_cell
from selling_global_wall_linearization_v50 import R, GT, GU, labs
from selling_reconstruction_engine_v49 import S, RELABEL

first={'g':9,'h':2,'k':1,'l':3,'m':5,'n':4}

def tup(A):return tuple(int(x) for x in list(A))
def wall_rep(A):
    # exact 6x6 action on homogeneous walls
    G,H,K,L,M,N=sp.symbols('g h k l m n'); sy=[G,H,K,L,M,N]
    a=-L-H-K;b=-M-G-K;c=-N-G-H
    X=sp.Matrix([[a,K,H],[K,b,G],[H,G,c]])
    Y=A.T*X*A
    aa,bb,cc=Y[0,0],Y[1,1],Y[2,2]; kk,hh,gg=Y[0,1],Y[0,2],Y[1,2]
    vals=[gg,hh,kk,-aa-hh-kk,-bb-gg-kk,-cc-gg-hh]
    Lm=sp.zeros(6,6)
    for i,e in enumerate(vals):
      e=sp.expand(e)
      for j,s in enumerate(sy):Lm[i,j]=e.coeff(s)
    return Lm

def powers(M,N):
    d={0:sp.eye(3)}
    for i in range(1,N+1):d[i]=d[i-1]*M
    Mi=M.inv()
    for i in range(1,N+1):d[-i]=d[-i+1]*Mi
    return d

def run(exp_bound=3,max_fields=200):
    PT=powers(GT,exp_bound);PU=powers(GU,exp_bound)
    deck={(a,b):PT[a]*PU[b] for a in PT for b in PU}
    reps=[sp.eye(3)]; records=[]; equivalences=[]
    # variant hash: equivalent path matrix -> representative and witness
    variants={}
    def add_variants(idx,A):
      for (a,b),D in deck.items():
       DA=D*A
       for bi,B in enumerate(RELABEL):
        X=DA*B
        variants.setdefault(tup(X),(idx,a,b,bi,1))
        variants.setdefault(tup(-X),(idx,a,b,bi,-1))
    add_variants(0,reps[0])
    qi=0
    while qi<len(reps):
      if len(reps)>max_fields:break
      A=reps[qi];L=wall_rep(A);cls=classify_cell(L)
      rec={'id':qi,'path_matrix':A.tolist(),'ordinary_walls':cls['ordinary_walls'],
           'multiwall_components':[{'factor':x['factor'],'walls':x['walls']} for x in cls['multiwall_components']],
           'neighbors':{}}
      for w in cls['ordinary_walls']:
        An=A*S[first[w]-1]; key=tup(An)
        if key in variants:
          j,a,b,bi,sgn=variants[key]
          rec['neighbors'][w]={'target':j,'new':False,'deck_exponents':[a,b],'relabel_index':bi,'sign':sgn}
          equivalences.append({'source':qi,'wall':w,'target':j,'a':a,'b':b,'relabel_index':bi,'sign':sgn})
        else:
          j=len(reps);reps.append(An);add_variants(j,An)
          rec['neighbors'][w]={'target':j,'new':True,'deck_exponents':[0,0],'relabel_index':0,'sign':1}
      records.append(rec);print('field',qi,'of',len(reps),'walls',cls['ordinary_walls'],'multi',[(x['walls'],x['factor']) for x in cls['multiwall_components']],flush=True);qi+=1
    # closure means all reps processed and no max hit
    closed=(qi==len(reps) and len(reps)<=max_fields)
    return {'version':'v50','status':'PASS-CLOSED' if closed else 'OPEN/TRUNCATED','exp_bound':exp_bound,
            'field_count':len(reps),'processed':qi,'closed':closed,'records':records,'equivalences':equivalences}

if __name__=='__main__':
 out=run(3,120)
 p=HERE/'exact_TU_field_bfs_v50.json';p.write_text(json.dumps(out,indent=2),encoding='utf-8')
 print(json.dumps({'status':out['status'],'field_count':out['field_count'],'processed':out['processed'],'closed':out['closed'],'sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
