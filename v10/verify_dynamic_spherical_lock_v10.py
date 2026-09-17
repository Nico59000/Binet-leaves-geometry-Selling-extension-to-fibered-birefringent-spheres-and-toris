#!/usr/bin/env python3
import sympy as sp, json
from pathlib import Path

# Symbols
K,Y,a,b,k,k3,C,h,H,rho,alpha,tau,q = sp.symbols('K Y a b k k3 C h H rho alpha tau q', positive=True)
D,nu,G,g = sp.symbols('D nu G g')

checks={}

# Exact shell map parameterized by a=e^{-2 alpha K t}
Phi=lambda y,aa: sp.simplify(K*y/(y+(K-y)*aa))
Z=lambda y: sp.simplify(K/y-1)
checks['shell_Z_linear']=sp.simplify(Z(Phi(Y,a))-a*Z(Y))==0
checks['shell_semigroup']=sp.simplify(Phi(Phi(Y,b),a)-Phi(Y,a*b))==0
checks['shell_fixed_difference']=sp.factor(Phi(Y,a)-K)==sp.factor(a*K*(Y-K)/(a*K+(1-a)*Y))
# fixed points: factor Phi-Y
fixed_factor=sp.factor(Phi(Y,a)-Y)
checks['shell_fixed_factor']=sp.simplify(fixed_factor - sp.factor(Y*(1-a)*(K-Y)/(a*K+(1-a)*Y)))==0

# shape normalized coordinate under relaxer
P,E=sp.symbols('P E', positive=True)
Pp=sp.Rational(1,2)*((1+k)*P+(1-k)*E)
Ep=sp.Rational(1,2)*((1-k)*P+(1+k)*E)
delta=sp.simplify((P-E)/(P+E))
deltap=sp.simplify((Pp-Ep)/(Pp+Ep))
checks['delta_scaled']=sp.simplify(deltap-k*delta)==0
checks['sum_preserved']=sp.simplify(Pp+Ep-P-E)==0

# even compensator: pre-relax D gets + C h, relax scales k, compensation removes k C h
D_before=k*(D+C*h)
D_after=sp.simplify(D_before-k*C*h)
checks['even_compensation']=sp.simplify(D_after-k*D)==0
# odd compensator
nu_after=sp.simplify(k3*(nu+G*g)-k3*G*g)
checks['odd_compensation']=sp.simplify(nu_after-k3*nu)==0

# radial compensation no-go: input K+Hh-rho; Phi=K iff rho=Hh via numerator
x=K+H*h-rho
expr=sp.factor(Phi(x,a)-K)
expected=sp.factor(a*K*(H*h-rho)/(K+(1-a)*(H*h-rho)))
checks['radial_forcing_identity']=sp.simplify(expr-expected)==0

# golden identity q = phi/2
phi=(1+sp.sqrt(5))/2
qg=sp.simplify(phi/2)
checks['golden_gap']=sp.simplify(qg**2-qg**3-sp.Rational(1,8))==0
# radial matched time tau2=-log q/(alpha K) -> exp(-2 alpha K tau2)=q^2
# verify algebraically by substituting exp(log) using positive q symbol structurally: exponent is 2 log q
checks['radial_match_exponent']=sp.simplify(-2*alpha*K*(-sp.log(q)/(alpha*K))-2*sp.log(q))==0

status='PASS' if all(checks.values()) else 'FAIL'
cert={
  'phase':'v10-global-dynamic-spherical-lock',
  'status':status,
  'checks':checks,
  'exact_shell_map':'Phi_t(Y)=K Y/[Y+(K-Y) exp(-2 alpha K t)]',
  'global_linear_coordinate':'Z=K/Y-1; Z(Phi_t(Y))=exp(-2 alpha K t) Z(Y)',
  'shape_coordinate':'delta=(P-E)/(P+E); delta -> kappa delta',
  'even_compensator':'subtract (kappa C h/2)*(1,-1) after relaxed injection',
  'odd_compensator':'subtract kappa_3 G g from the odd channel',
  'radial_fixed_point_guard':'under additive size forcing Hh, exact R_* locking requires rho=Hh',
  'golden_matching':'a_R=q^2, kappa=q^3 gives a_R-kappa=1/8',
  'guards':[
    'The compensators are declared feedback/feed-forward operators; they are not derived from Selling.',
    'Global radial attraction is proved on the positive radius sector for the exact shell flow, not for arbitrary discrete integrators.',
    'Shape and odd compensation can preserve centered isotropy while allowing the common radius to continue growing.',
    'An exact fixed radius during nonzero isotropic injection requires cancellation or rerouting of that radial injection.',
    'The golden choices a_R=q^2 and kappa=q^3 are typed closures, not canonical consequences of Binet/Selling/octagram resonance.'
  ]
}
Path('/mnt/data/dynamic_spherical_lock_v10_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False))
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
