import sympy as sp
import json, math

u=sp.symbols('u')
g,h,k=sp.symbols('g h k')
A,B,C,G,H,K,I0=sp.symbols('A B C G H K I0')
P=(u+g)*(u+h)*(u+k)
Q=(A*(u+h)**2*(u+k)**2
   +B*(u+k)**2*(u+g)**2
   +C*(u+g)**2*(u+h)**2
   +2*G*(u+g)**2*(u+h)*(u+k)
   +2*H*(u+h)**2*(u+g)*(u+k)
   +2*K*(u+k)**2*(u+g)*(u+h)
   +I0*P)
poly=sp.Poly(sp.expand(Q),u)
a4,a3,a2,a1,a0=poly.all_coeffs()

# Generic binary quartic invariant identity (kept independent of Selling expansion).
q4,q3,q2,q1,q0=sp.symbols('q4 q3 q2 q1 q0')
F=q4*u**4+q3*u**3+q2*u**2+q1*u+q0
I2=12*q4*q0-3*q3*q1+q2**2
J3=72*q4*q2*q0+9*q3*q2*q1-27*q4*q1**2-27*q3**2*q0-2*q2**3
Disc=sp.discriminant(F,u)
assert sp.expand(27*Disc-(4*I2**3-J3**2))==0

# Triple root lies at cusp vertex in invariant plane.
r,s=sp.symbols('r s')
tr=sp.Poly((u-r)**3*(u-s),u)
t4,t3,t2,t1,t0=tr.all_coeffs()
Itr=sp.expand(12*t4*t0-3*t3*t1+t2**2)
Jtr=sp.expand(72*t4*t2*t0+9*t3*t2*t1-27*t4*t1**2-27*t3**2*t0-2*t2**3)
assert sp.simplify(Itr)==0 and sp.simplify(Jtr)==0

phi=(1+sp.sqrt(5))/2
q=phi/2
assert sp.simplify(phi**2-(phi+1))==0
assert sp.simplify(q**2-q**3-sp.Rational(1,8))==0
mu,nu,b,ss=sp.symbols('mu nu b ss')
D=4*mu**3-27*b*nu**2
assert sp.expand(D.subs({mu:ss**2*mu,nu:ss**3*nu})-ss**6*D)==0
Delta=lambda n: sp.simplify(q**n-q**(n+1))
orth={1:sp.simplify(phi**2*Delta(0)),2:sp.simplify(phi*Delta(1)),3:sp.simplify(Delta(2))}
assert orth=={1:sp.Rational(1,2),2:sp.Rational(1,4),3:sp.Rational(1,8)}

phif=float(phi.evalf(30)); qf=float(q.evalf(30)); ef=math.e
samples={'q^2':qf**2,'q^(phi^2)':qf**(phif**2),'q^e':qf**ef,'q^3':qf**3}
assert samples['q^2'] > samples['q^(phi^2)'] > samples['q^e'] > samples['q^3']
ann=samples['q^2']-samples['q^3']
positions={'phi^2':(samples['q^2']-samples['q^(phi^2)'])/ann,'e':(samples['q^2']-samples['q^e'])/ann}

coeff_names=['a4','a3','a2','a1','a0']
coeff_exprs=[sp.factor(x) for x in (a4,a3,a2,a1,a0)]
out={
 'status':'PASS',
 'quartic_coefficients':dict(zip(coeff_names,map(str,coeff_exprs))),
 'quartic_invariants':{
   'I2':'12*a4*a0 - 3*a3*a1 + a2^2',
   'J3':'72*a4*a2*a0 + 9*a3*a2*a1 - 27*a4*a1^2 - 27*a3^2*a0 - 2*a2^3',
   'discriminant':'(4*I2^3 - J3^2)/27',
   'triple_root':'I2=J3=0 for (u-r)^3(u-s)'
 },
 'golden':{
   'phi^2_equals_phi_plus_1':True,'q':'phi/2','q2_minus_q3':'1/8',
   'Delta_n':'phi^(n-2)/2^(n+1)',
   'orthant_ladder':{str(d):str(v) for d,v in orth.items()},
   'samples':samples,'annulus_positions':positions
 },
 'A2_weighted_scaling':{'map':'(mu,nu)->(s^2 mu,s^3 nu)','D_scaling':'D->s^6 D','golden_step':'s=q gives q^2,q^3'}
}
with open('/mnt/data/selling_quartic_A2_golden_v4_certificate.json','w',encoding='utf-8') as f: json.dump(out,f,indent=2,ensure_ascii=False)
print(json.dumps({'status':'PASS','q2-q3':'1/8','orthant':out['golden']['orthant_ladder'],'samples':samples},indent=2))
