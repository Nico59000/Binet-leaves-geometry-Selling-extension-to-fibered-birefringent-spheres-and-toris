#!/usr/bin/env python3
from pathlib import Path
import json
import sympy as sp
ROOT=Path(__file__).resolve().parent
CP=json.loads((ROOT/"v54_GATE6_HYPERTENSOR_STRICT_BRIDGE_CHECKPOINT.json").read_text(encoding="utf-8"))
Q=sp.Matrix([[0,2,0],[2,0,0],[0,0,-1]])
S=sp.Matrix([[0,1,0],[1,0,0],[0,0,1]])
assert Q.det()==4 and S.T*Q*S==Q
x=sp.symbols("x",nonzero=True,real=True)
v=sp.Matrix([x,1/x,0])
Qp=sp.simplify(2*v*v.T-Q)
assert Qp==sp.diag(2*x**2,2/x**2,1)
Qpol=sp.diag(-1,-1,2)
Qp1=Qp.subs(x,1)
Q0=sp.simplify(Qp1-sp.trace(Qp1)/3*sp.eye(3))
assert Q0==-sp.Rational(1,3)*Qpol
n0=sp.sqrt(sp.trace(Q0.T*Q0)); np=sp.sqrt(sp.trace(Qpol.T*Qpol))
assert sp.simplify(Q0/n0+Qpol/np)==sp.zeros(3)
u=sp.symbols("u")
Qstruct=sp.expand(-4*u**4+4*(u+2)**2*u**2-12*u**2*(u+2))
assert sp.factor(Qstruct)==4*u**2*(u-2)
Frat=sp.simplify(-4*u**2/(u+2)+4*(u+2)-12)
assert sp.simplify(Frat-4*(u-2)/(u+2))==0
assert Frat.subs(u,2)==0 and sp.diff(Frat,u).subs(u,2)!=0
phi=(1+sp.sqrt(5))/2; q=phi/2
assert sp.simplify(q**2-q**3-sp.Rational(1,8))==0
assert CP["nucleation_test"]["direct_Pascal_point_to_A2_nucleation"]=="REFUTED-TYPED"
assert CP["Pascal_p46"]["exact_wrong_coefficient"].startswith("NT")
print("PASS")
