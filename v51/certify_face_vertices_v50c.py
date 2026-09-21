#!/usr/bin/env python3
from __future__ import annotations
import json,sys,hashlib,math
from pathlib import Path
import sympy as sp
HERE=Path('/mnt/data/v50_work');sys.path.insert(0,str(HERE))
from selling_global_wall_linearization_v50 import w0,labs,D0,q,r
BFS=json.load(open(HERE/'exact_TU_field_bfs_v50.json'));recs={x['id']:x for x in BFS['records']}
F=json.load(open(HERE/'selling_face_cycles_v50.json'));faces={x['field']:x for x in F['faces']}

def poly_prim(P,var):
 P=sp.Poly(P,var,domain=sp.QQ); cont,expr=sp.polys.polytools.primitive(P.as_expr(),var); P=sp.Poly(expr,var,domain=sp.QQ)
 if P.LC()<0:P=-P
 return P.sqf_part()
def intervals(P,var):
 P=poly_prim(P,var)
 return [(iv,m) for iv,m in sp.intervals(P,eps=sp.Rational(1,10)**24)]
def contains(iv,x):
 a,b=map(float,iv);return a-1e-12<=x<=b+1e-12
def root_sign(G,P,iv,var):
 G=sp.Poly(sp.cancel(G),var,domain=sp.QQ);P=poly_prim(P,var)
 if G.is_zero:return 0
 if iv[0]==iv[1]:
  v=G.eval(iv[0]);return 0 if v==0 else (1 if v>0 else -1)
 gd=sp.gcd(P,G)
 if gd.degree()>0:
  for giv,m in sp.intervals(gd,eps=sp.Rational(1,10)**28):
   if not (giv[1]<iv[0] or iv[1]<giv[0]):return 0
 groots=[iv2 for iv2,m in sp.intervals(G,eps=sp.Rational(1,10)**28)] if G.degree()>0 else []
 a,b=iv
 for _ in range(20):
  if not any(not (gb<a or b<ga) for ga,gb in groots):break
  if a==b:break
  a,b=sp.refine_root(P,a,b,eps=(b-a)/10)
 mid=(a+b)/2;val=G.eval(mid)
 if val>0:return 1
 if val<0:return -1
 raise RuntimeError(('sign',G,P,iv))

def triangular(f,g,qa,ra):
 for mode in ['q','r']:
  if mode=='q':
   G=sp.groebner([f,g],r,q,order='lex',domain=sp.QQ);var=q;other=r;approx=qa
  else:
   G=sp.groebner([f,g],q,r,order='lex',domain=sp.QQ);var=r;other=q;approx=ra
  basis=[sp.expand(p.as_expr()) for p in G.polys]
  unis=[];lins=[]
  for p in basis:
   P=sp.Poly(p,other,var,domain=sp.QQ)
   if P.degree(other)==0 and sp.Poly(p,var,domain=sp.QQ).degree()>0:unis.append(sp.Poly(p,var,domain=sp.QQ))
   if P.degree(other)==1:lins.append(p)
  if not unis or not lins:continue
  P=poly_prim(unis[-1],var)
  # choose isolating interval nearest approx
  ints=intervals(P,var)
  if not ints:continue
  cand=[(abs((float(iv[0])+float(iv[1]))/2-approx),iv,m) for iv,m in ints]
  _,iv,m=min(cand,key=lambda x:x[0])
  if not contains(iv,approx) and cand[0][0]>1e-5:continue
  # choose a linear relation with denominator nonzero at selected root and matching other approx
  for lin in lins:
   po=sp.Poly(lin,other,domain=sp.QQ[var])
   a=sp.Poly(po.coeff_monomial(other),var,domain=sp.QQ)
   b=sp.Poly(po.coeff_monomial(1),var,domain=sp.QQ)
   sa=root_sign(a,P,iv,var)
   if sa==0:continue
   expr=sp.cancel(-b.as_expr()/a.as_expr())
   mid=(iv[0]+iv[1])/2
   try:val=float(sp.N(expr.subs(var,mid),30))
   except:continue
   otherapprox=ra if mode=='q' else qa
   if abs(val-otherapprox)>1e-3:continue
   return {'mode':mode,'var':var,'other':other,'P':P,'iv':iv,'other_expr':expr,'groebner':basis}
 raise RuntimeError(('no triangular',f,g,qa,ra))

def eval_sign(expr,T):
 var=T['var'];other=T['other'];P=T['P'];iv=T['iv'];oe=T['other_expr']
 sub=sp.cancel(expr.subs(other,oe));num,den=sp.fraction(sp.together(sub));num=sp.Poly(num,var,domain=sp.QQ);den=sp.Poly(den,var,domain=sp.QQ)
 sn=root_sign(num,P,iv,var);sd=root_sign(den,P,iv,var)
 if sd==0:raise RuntimeError(('zero denominator',expr,T))
 return sn*sd

def rem_zero(expr,T):
 var=T['var'];other=T['other'];P=T['P'];oe=T['other_expr']
 sub=sp.cancel(expr.subs(other,oe));num,den=sp.fraction(sp.together(sub));
 return sp.Poly(num,var,domain=sp.QQ).rem(P).is_zero

cert=[];fails=[]
for i in range(62):
 face=faces[i];cy=face['cycles'][0];L=sp.Matrix(recs[i]['wall_action_matrix']);vec=L*w0;polys={labs[k]:sp.factor(vec[k]) for k in range(6)}
 n=len(cy['vertices'])
 for k,vlocal in enumerate(cy['vertices']):
  # vertex between previous and next boundary segments
  sprev=cy['segments'][(k-1)%n];snext=cy['segments'][k]
  f=sp.sympify(face['components'][sprev['component']]['factor'],locals={'q':q,'r':r})
  g=sp.sympify(face['components'][snext['component']]['factor'],locals={'q':q,'r':r})
  vv=face['vertices'][vlocal];qa=float(vv['q']);ra=float(vv['r'])
  try:
   T=triangular(f,g,qa,ra)
   assert eval_sign(f,T)==0 and eval_sign(g,T)==0
   signs={lab:eval_sign(p,T) for lab,p in polys.items()};sd=eval_sign(D0,T)
   ok=sd>0 and all(s<=0 for s in signs.values())
   if not ok:raise RuntimeError(('signs',i,k,signs,sd))
   rec={'field':i,'cycle_vertex_index':k,'local_vertex_id':vlocal,'approx':[qa,ra],
        'incident_component_indices':[sprev['component'],snext['component']],
        'incident_factors':[str(f),str(g)],'representation_variable':T['mode'],
        'root_polynomial':str(T['P'].as_expr()),'root_interval':[str(T['iv'][0]),str(T['iv'][1])],
        'other_coordinate_expr':str(T['other_expr']),'groebner_basis':[str(x) for x in T['groebner']],
        'wall_signs':signs,'D0_sign':sd,'active_walls':[lab for lab,s in signs.items() if s==0],
        'status':'PASS-EXACT'}
   cert.append(rec)
  except Exception as e:
   # exact rational fallback for reducible elimination/multiple projection roots
   done=False
   try:
    sols=sp.solve_poly_system([f,g],q,r)
    cand=[]
    for Q,Rr in sols:
     if Q.is_Rational and Rr.is_Rational:
      cand.append(((float(Q)-qa)**2+(float(Rr)-ra)**2,Q,Rr))
    if cand:
     _,Q,Rr=min(cand,key=lambda x:x[0])
     if (float(Q)-qa)**2+(float(Rr)-ra)**2<1e-8:
      signs={lab:(0 if sp.simplify(pp.subs({q:Q,r:Rr}))==0 else (1 if sp.simplify(pp.subs({q:Q,r:Rr}))>0 else -1)) for lab,pp in polys.items()}
      dv=sp.simplify(D0.subs({q:Q,r:Rr}));sd=0 if dv==0 else (1 if dv>0 else -1)
      assert sd>0 and all(v<=0 for v in signs.values())
      rec={'field':i,'cycle_vertex_index':k,'local_vertex_id':vlocal,'approx':[qa,ra],
           'incident_component_indices':[sprev['component'],snext['component']],
           'incident_factors':[str(f),str(g)],'representation_variable':'rational',
           'coordinates_exact':[str(Q),str(Rr)],'wall_signs':signs,'D0_sign':sd,
           'active_walls':[lab for lab,s0 in signs.items() if s0==0],'status':'PASS-EXACT-RATIONAL-FALLBACK'}
      cert.append(rec);done=True
   except Exception:
    pass
   if not done:
    fails.append({'field':i,'k':k,'local_vertex_id':vlocal,'error':repr(e),'f':str(f),'g':str(g),'approx':[qa,ra]});print('FAIL',fails[-1],flush=True)
 print('field',i,'certified',n,flush=True)
out={'version':'v50','status':'PASS' if not fails else 'FAIL','vertex_occurrences':len(cert),'failures':fails,'certificates':cert}
raw=json.dumps(out,indent=2).encode();out['sha256_pre_field']=hashlib.sha256(raw).hexdigest()
(HERE/'selling_face_vertex_exact_certificates_v50.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print('TOTAL',len(cert),'FAIL',len(fails))
