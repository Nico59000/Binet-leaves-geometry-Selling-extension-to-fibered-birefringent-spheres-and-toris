import sympy as sp

s5=sp.sqrt(5)
phi=(1+s5)/2
q=sp.simplify(phi/2)
assert sp.simplify(q**2-q**3-sp.Rational(1,8))==0

rho=sp.symbols('rho', positive=True)
alpha2=sp.simplify(rho*q**2)
alpha3=sp.simplify(rho*q**3)
assert sp.simplify(alpha2-alpha3-rho/sp.Integer(8))==0

# Exact increment law from definitions.
chi,zeta,w,R,b,z=sp.symbols('chi zeta w R b z', positive=True)
lam=sp.simplify(chi*w/2*(R**2+b**2/2))
xi=sp.simplify(zeta*w*(z**3+sp.Rational(3,2)*b**2*z))
lam2=sp.simplify(lam.subs({R:q*R,b:q*b,w:rho*w}))
xi2=sp.simplify(xi.subs({z:q*z,b:q*b,w:rho*w}))
assert sp.simplify(lam2-alpha2*lam)==0
assert sp.simplify(xi2-alpha3*xi)==0

# A2 discriminant covariance for rho=1.
mu,nu,beta=sp.symbols('mu nu beta')
D=4*mu**3-27*beta*nu**2
Dq=sp.expand(D.subs({mu:q**2*mu,nu:q**3*nu}, simultaneous=True))
assert sp.simplify(Dq-q**6*D)==0

# Phase-cycle lengths q_o=2,3 and block covariance.
for qo,L in [(2,4),(3,8)]:
    assert L==8//sp.gcd(8,qo)
    Dblock=sp.expand(D.subs({mu:q**(2*L)*mu,nu:q**(3*L)*nu}, simultaneous=True))
    assert sp.simplify(Dblock-q**(6*L)*D)==0

# Exact discriminant increment identity.
dmu,dnu=sp.symbols('dmu dnu')
Dnext=sp.expand(4*(mu+dmu)**3-27*beta*(nu+dnu)**2)
DeltaD=sp.expand(Dnext-D)
expected=12*mu**2*dmu+12*mu*dmu**2+4*dmu**3-54*beta*nu*dnu-27*beta*dnu**2
assert sp.expand(DeltaD-expected)==0

# Closed form geometric accumulation.
ell0,a,N=sp.symbols('ell0 a N', positive=True)
# symbolic finite sum checked for integer samples
for n in range(0,8):
    lhs=sum(q**(2*j) for j in range(n+1))
    rhs=(1-q**(2*(n+1)))/(1-q**2)
    assert sp.simplify(lhs-rhs)==0

print('PASS')
print('q=', sp.simplify(q))
print('q2-q3=', sp.simplify(q**2-q**3))
print('alpha2-alpha3=', sp.simplify(alpha2-alpha3))
print('D one-step factor (rho=1)=', sp.simplify(q**6))
print('D 4-step factor=', sp.simplify(q**24))
print('D 8-step factor=', sp.simplify(q**48))
