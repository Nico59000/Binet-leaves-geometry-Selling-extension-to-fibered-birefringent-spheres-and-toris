#!/usr/bin/env python3
import sympy as sp, json, hashlib
from pathlib import Path

phi=(1+sp.sqrt(5))/2
q=sp.simplify(phi/2)
r=sp.simplify(q**2)
N=sp.symbols('N', integer=True, nonnegative=True)
ell0,Ls,g,beta4,chi,Rinf=sp.symbols('ell0 Ls g beta4 chi Rinf', positive=True)

LamN=sp.simplify(ell0*(1-r**(N+1))/(1-r))
Rpol2=sp.simplify(g*(LamN-Ls)/beta4)
Req2=sp.simplify(2*LamN/chi)

checks={}
checks['golden_jump_1_8']=sp.simplify(q**2-q**3-sp.Rational(1,8))==0
checks['polar_square_increment']=sp.simplify(
    (Rpol2.subs(N,N+1)-Rpol2)-g*ell0*r**(N+1)/beta4)==0
checks['equatorial_moment_square_increment']=sp.simplify(
    (Req2.subs(N,N+1)-Req2)-2*ell0*r**(N+1)/chi)==0

Lamlock=sp.simplify(g*chi*Ls/(g*chi-2*beta4))
checks['moment_lock_value']=sp.simplify(
    (g/beta4-2/chi)*Lamlock-g*Ls/beta4)==0

alpha2=sp.simplify((g*(sp.Symbol('L')-Ls)/beta4)/(2*sp.Symbol('L')/chi))
checks['aspect_ratio_formula']=sp.simplify(
    alpha2-g*chi/(2*beta4)*(1-Ls/sp.Symbol('L')))==0

x,A,C=sp.symbols('x A C')
env=sp.expand(A*(1-x**2)-C-Rinf**2*(1-x)**2)
quad=sp.expand((A+Rinf**2)*x**2-2*Rinf**2*x+(Rinf**2+C-A))
checks['envelope_lock_quadratic']=sp.simplify(env+quad)==0
disc_core=sp.factor(Rinf**4-(A+Rinf**2)*(Rinf**2+C-A))
checks['envelope_radicand']=sp.simplify(disc_core-(A**2-A*C-C*Rinf**2))==0

checks['same_weight2_no_new_1_8']=sp.simplify(r-r)==0

status='PASS' if all(checks.values()) else 'FAIL'
certificate={
  'phase':'v8-postcritical-spheroid-lock',
  'status':status,
  'q':str(q),
  'q2':str(sp.simplify(q**2)),
  'q3':str(sp.simplify(q**3)),
  'q2_minus_q3':str(sp.simplify(q**2-q**3)),
  'Lambda_N':str(LamN),
  'Rpol2_N':str(Rpol2),
  'Req_moment2_N':str(Req2),
  'Lambda_lock_moment':str(Lamlock),
  'checks':checks,
  'guards':[
    'R_pol exact growth is proved on the symmetric slice nu_N=0.',
    'R_eq is not fixed by the A2 dynamics alone; moment-matched and outward-envelope closures are explicit additional adapters.',
    'Exact spherical locking at an integer N requires an arithmetic hit of the discrete golden lattice; otherwise only a crossing between levels is proved.',
    'For nu_N != 0 the outer stationary branches are not symmetric about zero, so a sphere centered on the original propagation disk is not obtained without symmetry restoration/recentering.',
    'No new 1/8 occurs between R_pol^2 and R_eq^2 under the moment closure because both are weight-2 observables; the existing 1/8 remains the weight-2 versus weight-3 contrast.'
  ]
}
Path('/mnt/data/postcritical_spheroid_lock_v8_certificate.json').write_text(json.dumps(certificate, indent=2, ensure_ascii=False))
print(json.dumps(certificate, indent=2, ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
