#!/usr/bin/env python3
import json, sympy as sp
from pathlib import Path
ROOT=Path(__file__).resolve().parent
cp=json.loads((ROOT/'v54_GATE6_RESEARCH_CHECKPOINT.json').read_text(encoding='utf-8'))
checks=[]
def ok(name, cond):
    checks.append((name,bool(cond)))
    if not cond: raise AssertionError(name)
# Fig2 exact source algebra
Q=sp.Matrix(cp['Fig2']['source_form_Q']); F=sp.Matrix(cp['Fig2']['adjoint_F'])
R=sp.Matrix(cp['Fig2']['R']); S=sp.Matrix(cp['Fig2']['S'])
D=sp.Matrix(cp['Fig2']['Delta']); P=sp.Matrix(cp['Fig2']['Phi'])
I=sp.eye(3)
ok('Fig2_det5', Q.det()==5)
ok('Fig2_adjoint', F==Q.adjugate())
for name,A in [('R',R),('S',S)]:
    ok(name+'_stabilizes_Q', A.T*Q*A==Q)
    ok(name+'_det1', A.det()==1)
ok('R2',R**2==I); ok('S2',S**2==I); ok('D2',D**2==I); ok('P2',P**2==I)
ok('RS_order3',(R*S)**3==I)
ok('S_Phi_commute',S*P==P*S)
ok('R_Delta_commute',R*D==D*R)
# point transport f(u)/F(v)
u0=sp.Matrix([1,0,0]);v0=sp.Matrix([1,0,0])
u1=R*u0;v1=R.inv().T*v0
# signs irrelevant for quadratic evaluation
ok('P1_num_projective',u1 in (sp.Matrix([1,1,0]),sp.Matrix([-1,-1,0])))
ok('P1_den_projective',v1 in (sp.Matrix([1,0,0]),sp.Matrix([-1,0,0])))
u2=S*u1;v2=S.inv().T*v1
U2=sp.Matrix(cp['Fig2']['derived_P2']['numerator_vector']);V2=sp.Matrix(cp['Fig2']['derived_P2']['denominator_vector'])
ok('P2_num_projective',u2 in (U2,-U2));ok('P2_den_projective',v2 in (V2,-V2))
def qv(M,v):return sp.expand((v.T*M*v)[0])
for name,u,v in [('P0',u0,v0),('P1',u1,v1),('P2',u2,v2)]:
    ok(name+'_ratio_2_2',qv(Q,u)==2 and qv(F,v)==2)
# Pascal lane B exact reduction
QP=sp.Matrix(cp['Pascal_Noether_lane_B']['Q_primitive'])
ok('Pascal_det4',QP.det()==4)
ok('Pascal_signature',sorted([x for x in QP.eigenvals().keys()])==[-2,-1,2])
x,y,z=sp.symbols('x y z', real=True, nonzero=True)
v=sp.Matrix([x,y,z]); hyp=sp.expand((v.T*QP.inv()*v)[0])
ok('Pascal_hyperboloid',sp.expand(hyp-(x*y-z**2))==0)
Qplus=2*v*v.T-QP
a,b,c=Qplus[0,0],Qplus[1,1],Qplus[2,2]
g,h,k=Qplus[1,2],Qplus[0,2],Qplus[0,1]
l=sp.expand(-a-h-k);m=sp.expand(-b-g-k);n=sp.expand(-c-g-h)
ok('Pascal_k_on_hyperboloid',sp.rem(sp.Poly(k-2*z**2,x*y if False else z),sp.Poly(1,z))==sp.Poly(k-2*z**2,z) if False else sp.expand(k.subs(x*y,1+z**2)-2*z**2)==0)
# direct substitution z=0,y=1/x after reduction
subs={z:0,y:1/x}
vals=[sp.simplify(t.subs(subs)) for t in (g,h,k,l,m,n)]
ok('Pascal_reduced_walls',vals==[0,0,0,-2*x**2,-2/x**2,-1])
Qred=sp.simplify(Qplus.subs(subs))
ok('Pascal_Qplus_diag',Qred==sp.diag(2*x**2,2/x**2,1))
# guards/statuses
ok('Fig3_2of91',cp['Fig3']['direct_raster_bound']=='2/91' and cp['Fig3']['not_propagated']=='89/91')
ok('Fig2_raster_stays_NT',cp['Fig2']['S_specific_raster_edge'].startswith('NT') and cp['Fig2']['Phi_raster_species'].startswith('NT'))
ok('laneA_coefficient_NT',cp['Pascal_Noether_lane_A']['exact_wrong_coefficient'].startswith('NT'))
ok('R2_r4_closed',cp['R2PURE_13']['radius4']=='NOT-OPENED' and cp['R2PURE_13']['B6']=='NOT-OPENED')
ok('Gate6_nonpublishable',cp['publication_state']=='RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE')
print(json.dumps({'status':'PASS','checks':len(checks),'passed':sum(v for _,v in checks),'failed':[n for n,v in checks if not v]},indent=2))
