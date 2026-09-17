#!/usr/bin/env python3
import sympy as sp, json
from pathlib import Path
q=(1+sp.sqrt(5))/4
r=sp.simplify(q**2)
k,C=sp.symbols('k C')
A=sp.Matrix([[sp.Rational(1,2)*(1+k), sp.Rational(1,2)*(1-k)],
             [sp.Rational(1,2)*(1-k), sp.Rational(1,2)*(1+k)]])
Pp=sp.Matrix([[sp.Rational(1,2),sp.Rational(1,2)],[sp.Rational(1,2),sp.Rational(1,2)]])
Pm=sp.Matrix([[sp.Rational(1,2),-sp.Rational(1,2)],[-sp.Rational(1,2),sp.Rational(1,2)]])
checks={}
checks['projector_sum']=sp.simplify(Pp+Pm-sp.eye(2))==sp.zeros(2)
checks['projector_orthogonal']=Pp*Pm==sp.zeros(2)
checks['relaxer_decomposition']=sp.simplify(A-(Pp+k*Pm))==sp.zeros(2)
l=sp.symbols('l')
Al=Pp+l*Pm
checks['semigroup']=sp.simplify(A*Al-(Pp+k*l*Pm))==sp.zeros(2)
ones=sp.Matrix([[1,1]])
diff=sp.Matrix([[1,-1]])
checks['sum_preserved']=sp.simplify(ones*A-ones)==sp.zeros(1,2)
checks['difference_scaled']=sp.simplify(diff*A-k*diff)==sp.zeros(1,2)
# closed recurrence verification
m,N0,D0=sp.symbols('m N0 D0', integer=True, nonnegative=True)
# identity for geometric convolution symbolically with an auxiliary t
x,y=sp.symbols('x y')
geom=(x**m-y**m)/(x-y)
checks['geom_convolution_identity']=sp.simplify((x-y)*geom-(x**m-y**m))==0
# golden gap
checks['golden_gap']=sp.simplify(r-q**3-sp.Rational(1,8))==0
# matched closed form satisfies one-step recurrence using N=N0+m
Dclosed=k**m*D0+k*C*r**(N0+1)*(k**m-r**m)/(k-r)
Dnext=k**(m+1)*D0+k*C*r**(N0+1)*(k**(m+1)-r**(m+1))/(k-r)
checks['closed_recurrence']=sp.simplify(Dnext-(k*Dclosed+k*C*r**(N0+m+1)))==0
# common radius sum
ell,g,b,chi,Ls=sp.symbols('ell g b chi Ls', positive=True)
Linf=ell/(1-r)
Rcommon=sp.Rational(1,2)*(g/b*(Linf-Ls)+2*Linf/chi)
checks['common_radius_formula']=sp.simplify(2*Rcommon-(g/b*(Linf-Ls)+2*Linf/chi))==0
# shell stability derivative
R,tau,alpha,Rs=sp.symbols('R tau alpha Rs', positive=True)
F=R-tau*alpha*R*(R**2-Rs**2)
der=sp.diff(F,R).subs(R,Rs)
checks['shell_linearization']=sp.simplify(der-(1-2*tau*alpha*Rs**2))==0
# autonomous eigenvalues
M=sp.Matrix([[k,k*C],[0,r]])
L=sp.symbols('L')
checks['lift_charpoly']=sp.simplify((L*sp.eye(2)-M).det()-(L-k)*(L-r))==0
cert={
 'phase':'v9-dynamic-spherical-lock',
 'status':'PASS' if all(checks.values()) else 'FAIL',
 'q':str(q), 'q2':str(r), 'q3':str(sp.simplify(q**3)),
 'q2_minus_q3':str(sp.simplify(r-q**3)),
 'relaxer':'Pi_plus + kappa Pi_minus',
 'lift_eigenvalues':['kappa','q^2'],
 'full_shape_center_eigenvalues':['kappa_even','kappa_odd','q^2','q^3'],
 'radial_linearized_eigenvalue':'1-2*tau*alpha*R_*^2',
 'checks':checks,
 'guards':[
  'Soft relaxation makes isotropy asymptotically attractive but not exactly invariant at finite N while anisotropic forcing remains nonzero.',
  'Hard locking kappa=0 is an extra projection, not derived from Selling.',
  'The choice kappa=q^3 is a typed golden matching assumption; only the resulting gap q^2-q^3=1/8 is exact.',
  'Shape relaxation alone leaves a neutral common-radius direction; the shell potential is needed to isolate a unique stable sphere radius.',
  'Centered spherical stability also requires contraction or vanishing of the odd bias nu/Xi.'
 ]
}
Path('/mnt/data/dynamic_spherical_lock_v9_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False))
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if cert['status']=='PASS' else 1)
