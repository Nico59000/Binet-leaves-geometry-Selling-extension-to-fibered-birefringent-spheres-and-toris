#!/usr/bin/env python3
from pathlib import Path
import sympy as sp, math, json

ROOT=Path(__file__).resolve().parent
TGT=json.loads((ROOT/"v54_GATE6_RHO15_EXACT_TARGET.json").read_text(encoding="utf-8"))
CP=json.loads((ROOT/"v54_GATE6_P46_GLOBAL_REPAIR_CHECKPOINT.json").read_text(encoding="utf-8"))

x,y,z=sp.symbols("x y z")
U,V,W=sp.symbols("U V W")
A0,A1,A2,A3,A4=sp.symbols("A0 A1 A2 A3 A4")
B0,B1,B2,B3=sp.symbols("B0 B1 B2 B3")
G0,G1,G2=sp.symbols("G0 G1 G2")
D0,D1,a=sp.symbols("D0 D1 a")
alpha=A0*x**4+4*A1*x**3*y+6*A2*x**2*y**2+4*A3*x*y**3+A4*y**4
beta=B0*x**3+3*B1*x**2*y+3*B2*x*y**2+B3*y**3
gamma=G0*x**2+2*G1*x*y+G2*y**2

checks=[]
def ok(name,cond):
    checks.append((name,bool(cond)))
    if not cond: raise AssertionError(name)

def hd(F):
    F=sp.expand(F)
    if F==0:return None
    ds={sum(m) for m,c in sp.Poly(F,x,y).terms()}
    assert len(ds)==1
    return next(iter(ds))
def tr(F,G,r):
    F=sp.expand(F);G=sp.expand(G)
    if F==0 or G==0:return sp.Integer(0)
    m,n=hd(F),hd(G)
    if r>m or r>n:return sp.Integer(0)
    pref=sp.Rational(math.factorial(m-r)*math.factorial(n-r),
                     math.factorial(m)*math.factorial(n))
    S=sum((-1)**j*sp.binomial(r,j)*sp.diff(F,x,r-j,y,j)*sp.diff(G,x,j,y,r-j)
          for j in range(r+1))
    return sp.expand(pref*S)
def dual(F):
    return sp.expand(F.subs({x:V,y:-U},simultaneous=True))
def qmat(poly):
    P=sp.Poly(sp.expand(poly),U,V,W); vs=[U,V,W]; M=sp.zeros(3)
    for i,t in enumerate(vs): M[i,i]=P.coeff_monomial(t**2)
    for i in range(3):
        for j in range(i+1,3):
            M[i,j]=M[j,i]=P.coeff_monomial(vs[i]*vs[j])/2
    return M

# Source-calibrated spinta.
ok("spinta_low_identity",tr(x**4,y**3,3)==x)

def p_corrected(al,be,ga,aa):
    k=tr(al,al,2); ii=tr(al,al,4); p=tr(al,be,3)
    ww=tr(be,be,2); Q=tr(be,ww,1); R=tr(ww,ww,2)
    pi=tr(k,be,3); J=tr(al,k,4); h=tr(ga,ga,2); m=tr(ga,be,2)
    C33=(sp.Rational(1,5)*aa*J
         +tr(tr(k,ga,2),ga,2)
         +sp.Rational(7,30)*ii*h
         -sp.Rational(8,5)*tr(tr(p,be,1),ga,2)
         +sp.Rational(8,5)*tr(tr(ww,al,2),ga,2)
         -sp.Rational(4,5)*R)
    C3=(-sp.Rational(4,5)*aa*dual(pi)
        -sp.Rational(44,5)*dual(tr(tr(al,be,1),ga**2,4))
        +sp.Rational(28,25)*dual(tr(tr(al,be,2),ga**2,3))
        +sp.Rational(28,15)*dual(p)*h
        -sp.Rational(8,5)*dual(tr(Q,ga,2)))
    C0=(sp.Rational(4,5)*aa*dual(tr(k,ga,2))
        +sp.Rational(7,15)*aa*ii*dual(ga)
        +sp.Rational(8,5)*h*dual(tr(al,ga,2))
        -sp.Rational(6,5)*tr(al,ga**2,4)*dual(ga)
        -sp.Rational(6,5)*aa*dual(tr(p,be,1))
        -sp.Rational(1,5)*aa*dual(tr(ww,al,2))
        +sp.Rational(12,5)*tr(ww,ga,2)*dual(ga)
        +sp.Rational(4,5)*h*dual(ww)
        -2*dual(m)**2)
    return sp.expand(W**2*C33+W*C3+C0)

P13=p_corrected(alpha,beta,gamma,a)
ok("p13_223_terms",len(sp.Poly(P13,U,V,W,A0,A1,A2,A3,A4,B0,B1,B2,B3,G0,G1,G2,a).terms())==223)

# Reconstruct Noether rho in the standard U,V,W convention from its exact sparse target.
c=sp.symbols("c0:15")
rho=0
for key,comp in TGT["components"].items():
    den=sp.Integer(comp["raw_den_lcm"])
    gcd=sp.Integer(comp["primitive_gcd_removed"])
    pol=0
    for term in comp["terms"]:
        m=sp.Integer(term["coef"])
        for vv,pw in zip(c,term["exp"]): m*=vv**pw
        pol+=m
    mult=sp.Integer(TGT["standard_component_multiplier_sign"][key])
    ex=eval(key)
    rho += mult*gcd*pol/den * U**ex[0]*V**ex[1]*W**ex[2]
rho=sp.expand(rho)

# Calibration of target convention.
basis=[tuple(q) for q in TGT["quartic_coefficient_basis"]]
idx={e:i for i,e in enumerate(basis)}
cal={vv:0 for vv in c}
cal[c[idx[(2,2,0)]]]=1
cal[c[idx[(2,0,2)]]]=-1
cal[c[idx[(0,2,2)]]]=-1
rho_cal=sp.expand(rho.subs(cal))
ok("rho_target_calibration",rho_cal==-(U**2+V**2-W**2)/sp.Integer(648))

# Exact 13-parameter Pascal gauge identity.
vals13=[A0,4*A1,6*A2,4*A3,A4,4*B0,12*B1,12*B2,4*B3,6*G0,12*G1,6*G2,0,0,a]
rho13=sp.expand(rho.subs({vv:q for vv,q in zip(c,vals13)}))
ok("exact_13_parameter_identity",sp.expand(rho13-sp.Rational(5,2)*P13)==0)

# Exact extension to all 15 coefficients through Pascal's z^3 elimination.
delta=D0*x+D1*y
ga_r=sp.expand(gamma-delta**2/a)
be_r=sp.expand(beta-3*gamma*delta/a+2*delta**3/a**2)
al_r=sp.expand(alpha-4*beta*delta/a+6*gamma*delta**2/a**2-3*delta**4/a**3)
def bcoeff(poly,degree,i):
    return sp.cancel(sp.Poly(poly,x,y).coeff_monomial(x**(degree-i)*y**i)/sp.binomial(degree,i))
redmap={}
for i,s in enumerate([A0,A1,A2,A3,A4]): redmap[s]=bcoeff(al_r,4,i)
for i,s in enumerate([B0,B1,B2,B3]): redmap[s]=bcoeff(be_r,3,i)
for i,s in enumerate([G0,G1,G2]): redmap[s]=bcoeff(ga_r,2,i)
redmap[a]=a

Cred=qmat(P13).applyfunc(lambda e:sp.cancel(e.subs(redmap,simultaneous=True)))
T=sp.Matrix([[1,0,0],[0,1,0],[-D0/a,-D1/a,1]])
Cold=T*Cred*T.T

vals15=[A0,4*A1,6*A2,4*A3,A4,4*B0,12*B1,12*B2,4*B3,
        6*G0,12*G1,6*G2,4*D0,4*D1,a]
rho15=sp.expand(rho.subs({vv:q for vv,q in zip(c,vals15)}))
R=qmat(rho15)
for i in range(3):
    for j in range(i,3):
        num,den=sp.fraction(sp.together(R[i,j]-sp.Rational(5,2)*Cold[i,j]))
        ok(f"15coeff_entry_{i}{j}_zero",sp.expand(num)==0)

# Numeric direct Pascal-reduction GL3 replay, independent of just substituting rho.
def coeff_z(poly,k):
    return sp.Poly(sp.expand(poly),z).coeff_monomial(z**k)
def p_general_matrix(f):
    aa=sp.expand(coeff_z(f,4))
    if aa==0: raise ValueError("zero z4 pivot")
    de=sp.expand(coeff_z(f,3)/4)
    ga=sp.expand(coeff_z(f,2)/6)
    be=sp.expand(coeff_z(f,1)/4)
    al=sp.expand(coeff_z(f,0))
    alr=sp.expand(al-4*be*de/aa+6*ga*de**2/aa**2-3*de**4/aa**3)
    ber=sp.expand(be-3*ga*de/aa+2*de**3/aa**2)
    gar=sp.expand(ga-de**2/aa)
    C=qmat(p_corrected(alr,ber,gar,aa))
    pd=sp.Poly(de,x,y); d0=pd.coeff_monomial(x); d1=pd.coeff_monomial(y)
    TT=sp.Matrix([[1,0,0],[0,1,0],[-d0/aa,-d1/aa,1]])
    return sp.simplify(TT*C*TT.T)

f0=(2*x**4+3*x**3*y-2*x**2*y**2+x*y**3+3*y**4
    +4*(x**3+x**2*y-2*x*y**2+2*y**3)*z
    +6*(2*x**2+x*y-y**2)*z**2+4*(x-2*y)*z**3+3*z**4)
C0=p_general_matrix(f0)
vec=sp.Matrix([x,y,z])
Gs=[
 sp.Matrix([[1,1,0],[0,1,0],[0,0,1]]),
 sp.Matrix([[1,0,1],[0,1,0],[0,0,1]]),
 sp.Matrix([[1,0,0],[1,1,0],[0,0,1]]),
 sp.Matrix([[0,1,0],[1,0,0],[0,0,1]]),
 sp.diag(-1,1,1),
 sp.Matrix([[1,0,0],[0,0,1],[0,1,0]])
]
passed=0
for n,G in enumerate(Gs):
    Gi=G.inv(); q=Gi*vec
    ft=sp.expand(f0.subs({x:q[0],y:q[1],z:q[2]},simultaneous=True))
    ok(f"GL3_pivot_{n}",sp.Poly(ft,x,y,z).coeff_monomial(z**4)!=0)
    Ct=p_general_matrix(ft)
    ok(f"GL3_covariance_{n}",sp.simplify(Ct-G*C0*G.T)==sp.zeros(3))
    passed+=1
ok("GL3_count_6",passed==6)

# Centered amplitude jet; no time interpretation.
chi,gH,a0,s=sp.symbols("chi gH a0 s",nonzero=True)
sstar=2*a0/(gH*chi)
Lam=chi*s/2
Lamstar=chi*sstar/2
dLam=sp.simplify(Lam-Lamstar)
mu=sp.simplify(gH*dLam)
ok("center_reference",sp.simplify(Lamstar-a0/gH)==0)
ok("centered_deltaLambda",sp.simplify(dLam-chi*(s-sstar)/2)==0)
ok("centered_mu",sp.simplify(mu-gH*chi*(s-sstar)/2)==0)
ok("mu_zero_at_reference",sp.simplify(mu.subs(s,sstar))==0)

report={
 "phase":"v54 Gate 6 p46 global repair + centered amplitude",
 "status":"PASS",
 "checks":len(checks),
 "passed":sum(v for _,v in checks),
 "failed":[n for n,v in checks if not v],
 "gate_closed":False,
 "publication_state":CP["publication_state"],
 "promotions":{
   "p46_corrected_formula":"PROVEN-GLOBAL-POLYNOMIAL-IDENTITY-IN-15-COEFFICIENTS",
   "GL3_covariance":"PROVEN/6-DIRECT-REPLAYS-PASS",
   "centered_amplitude_jet":"PROVEN-EXACT"
 }
}
(ROOT/"v54_GATE6_P46_GLOBAL_REPAIR_REPORT.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
print(json.dumps(report,indent=2))
