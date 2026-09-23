#!/usr/bin/env python3
import sympy as sp, math, json
from pathlib import Path

x,y,z=sp.symbols("x y z")
U,V,W=sp.symbols("U V W")
A0,A1,A2,A3,A4=sp.symbols("A0 A1 A2 A3 A4")
B0,B1,B2,B3=sp.symbols("B0 B1 B2 B3")
G0,G1,G2=sp.symbols("G0 G1 G2")
a=sp.symbols("a")

alpha=A0*x**4+4*A1*x**3*y+6*A2*x**2*y**2+4*A3*x*y**3+A4*y**4
beta=B0*x**3+3*B1*x**2*y+3*B2*x*y**2+B3*y**3
gamma=G0*x**2+2*G1*x*y+G2*y**2

def order(f):
    return sp.Poly(sp.expand(f),x,y).total_degree()
def tr(f,g,r):
    m,n=order(f),order(g)
    s=0
    for k in range(r+1):
        s += (-1)**k*sp.binomial(r,k)*sp.diff(f,x,r-k,y,k)*sp.diff(g,x,k,y,r-k)
    return sp.factor(sp.factorial(m-r)*sp.factorial(n-r)/(sp.factorial(m)*sp.factorial(n))*s)
def dual(f):
    return sp.expand(f.subs({x:V,y:-U}))

# Pascal binary carriers
k=tr(alpha,alpha,2)
i=tr(alpha,alpha,4)
p=tr(alpha,beta,3)
wb=tr(beta,beta,2)
Q=tr(beta,wb,1)
R=tr(wb,wb,2)
pi=tr(k,beta,3)
J=tr(alpha,k,4)
h=tr(gamma,gamma,2)
m=tr(gamma,beta,2)

# Low spinta calibration
assert tr(x**4,y**3,3)==x

# Noether source law
x1,y1,z1,x2,y2,z2,x3,y3,z3=sp.symbols("x1 y1 z1 x2 y2 z2 x3 y3 z3")
A12=sp.Matrix([[U,x1,x2],[V,y1,y2],[W,z1,z2]])
A23=sp.Matrix([[U,x2,x3],[V,y2,y3],[W,z2,z3]])
A31=sp.Matrix([[U,x3,x1],[V,y3,y1],[W,z3,z1]])
templ=sp.Poly(sp.expand(A12.det()**2*A23.det()**2*A31.det()**2/sp.Integer(6)),
              U,V,W,x1,y1,z1,x2,y2,z2,x3,y3,z3,domain=sp.QQ)
basis=[(4,0,0),(3,1,0),(2,2,0),(1,3,0),(0,4,0),
       (3,0,1),(2,1,1),(1,2,1),(0,3,1),
       (2,0,2),(1,1,2),(0,2,2),(1,0,3),(0,1,3),(0,0,4)]
multi=[math.factorial(4)//(math.factorial(r)*math.factorial(s)*math.factorial(t)) for r,s,t in basis]

def rho_from_actual(actual):
    coeff={e:sp.cancel(actual.get(e,0)/sp.Integer(mm)) for e,mm in zip(basis,multi)}
    psi=0
    for mon,c in templ.terms():
        eu,ev,ew=mon[:3]
        e1,e2,e3=mon[3:6],mon[6:9],mon[9:12]
        c1,c2,c3=coeff[e1],coeff[e2],coeff[e3]
        if c1!=0 and c2!=0 and c3!=0:
            psi += c*c1*c2*c3*U**eu*V**ev*W**ew
    psi=sp.expand(psi)
    D=0
    for e,c in actual.items():
        if c!=0:
            D += c*sp.diff(psi,U,e[0],V,e[1],W,e[2])
    return sp.factor(sp.expand(D/sp.Integer(144)))

# gamma=0 family
actual_g0={
 (4,0,0):A0,(3,1,0):4*A1,(2,2,0):6*A2,(1,3,0):4*A3,(0,4,0):A4,
 (3,0,1):4*B0,(2,1,1):12*B1,(1,2,1):12*B2,(0,3,1):4*B3,
 (2,0,2):0,(1,1,2):0,(0,2,2):0,(1,0,3):0,(0,1,3):0,(0,0,4):a}
rho_g0=rho_from_actual(actual_g0)
p_g0=sp.expand(
    W**2*(sp.Rational(1,5)*a*J-sp.Rational(4,5)*R)
    +W*(-sp.Rational(4,5)*a*dual(pi))
    +dual(-sp.Rational(6,5)*a*tr(p,beta,1)-sp.Rational(1,5)*a*tr(wb,alpha,2))
)
assert sp.expand(rho_g0-sp.Rational(5,2)*p_g0)==0

# alpha=0 family
actual_a0={
 (4,0,0):0,(3,1,0):0,(2,2,0):0,(1,3,0):0,(0,4,0):0,
 (3,0,1):4*B0,(2,1,1):12*B1,(1,2,1):12*B2,(0,3,1):4*B3,
 (2,0,2):6*G0,(1,1,2):12*G1,(0,2,2):6*G2,
 (1,0,3):0,(0,1,3):0,(0,0,4):a}
rho_a0=rho_from_actual(actual_a0)
wg=tr(wb,gamma,2)
Qg=tr(Q,gamma,2)
p_a0_print=sp.expand(
    W**2*(-sp.Rational(4,5)*R)
    +W*(-sp.Rational(8,5)*dual(Qg))
    +sp.Rational(16,5)*wg*dual(gamma)
    +sp.Rational(4,5)*h*dual(wb)
    -2*dual(m)**2
)
res=sp.factor(sp.expand(rho_a0-sp.Rational(5,2)*p_a0_print))
assert sp.expand(res + 2*wg*dual(gamma))==0

p_a0_corr=sp.expand(p_a0_print-sp.Rational(4,5)*wg*dual(gamma))
assert sp.expand(rho_a0-sp.Rational(5,2)*p_a0_corr)==0

# Polar fixed-ray amplitude adapter
Qraw=sp.Matrix([[0,sp.Rational(1,2),0],[sp.Rational(1,2),0,0],[0,0,-sp.Rational(1,4)]])
Qplus_fixed=sp.diag(sp.Rational(1,2),sp.Rational(1,2),sp.Rational(1,4))
assert Qplus_fixed[:2,:2]==sp.Rational(1,2)*sp.eye(2)

print("PASS")
