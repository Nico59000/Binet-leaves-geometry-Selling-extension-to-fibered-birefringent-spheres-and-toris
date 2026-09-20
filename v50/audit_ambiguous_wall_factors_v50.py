#!/usr/bin/env python3
import json,sys,sympy as sp,time
from pathlib import Path
HERE=Path('/mnt/data/v50_work');sys.path.insert(0,str(HERE))
from selling_global_wall_linearization_v50 import w0,labs,D0
import exact_boundary_arc_classifier_v50 as c
q,r=c.q,c.r
D=json.load(open(HERE/'exact_TU_field_bfs_v50.json'))

def norm(expr):
 F=sp.Poly(expr,q,r,domain=sp.QQ)
 if F.LC()<0:F=-F
 return sp.factor(F.monic().as_expr())

def factor_valid(f,polys):
 vv=c.classify_vertical_factor(f,polys,[])
 if vv:return {'factor':str(f),'valid':True,'method':'vertical','samples':vv}
 C=c.critical_q_poly(f,list(polys.values())+[D0]);qs=c.q_samples_from_critical(C);valid=[]
 for q0 in qs:
  fq=sp.Poly(sp.expand(f.subs(q,q0)),r,domain=sp.QQ)
  if fq.is_zero or fq.degree()<=0:continue
  for iv,mult in sp.intervals(fq,eps=sp.Rational(1,10)**12):
   signs={lab:c.root_sign(sp.expand(p.subs(q,q0)),fq,iv) for lab,p in polys.items()}
   sd=c.root_sign(sp.expand(D0.subs(q,q0)),fq,iv)
   if sd>0 and all(v<=0 for v in signs.values()):valid.append({'q':str(q0),'r_interval':[str(iv[0]),str(iv[1])],'signs':signs})
 return {'factor':str(f),'valid':bool(valid),'method':'sample','samples':valid}

out=[]
for rec in D['records']:
 L=sp.Matrix(rec['wall_action_matrix']);vec=L*w0;polys={labs[i]:sp.factor(vec[i]) for i in range(6)}
 for w in rec['ordinary_walls']:
  fac=c.primitive_factor_list(polys[w])
  if len(fac)<=1:continue
  t=time.time();R=[]
  for f,m in fac:R.append(factor_valid(f,polys))
  rr={'field':rec['id'],'wall':w,'poly':str(polys[w]),'factors':R,'seconds':time.time()-t}
  out.append(rr);print(rec['id'],w,[(x['factor'],x['valid']) for x in R],'sec',round(rr['seconds'],2),flush=True)
(HERE/'selling_ambiguous_wall_factor_audit_v50.json').write_text(json.dumps({'version':'v50','records':out},indent=2),encoding='utf-8')
