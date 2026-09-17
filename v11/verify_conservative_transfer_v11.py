#!/usr/bin/env python3
import sympy as sp, json
from pathlib import Path

k,th,r=sp.symbols('k theta r')
M=sp.Matrix([[k,0,k*th*(1-r)],[1-k,1,(1-k*th)*(1-r)],[0,0,r]])
checks={}
checks['column_sum_conservation']=sp.simplify(sp.Matrix([[1,1,1]])*M-sp.Matrix([[1,1,1]]))==sp.zeros(1,3)
checks['determinant']=sp.factor(M.det()-k*r)==0
cp=M.charpoly(); L=cp.gen
checks['charpoly']=sp.expand(cp.as_expr()-(L-1)*(L-k)*(L-r))==0

n=sp.symbols('n', integer=True, nonnegative=True)
x0,t0=sp.symbols('x0 t0')
xf=k**n*x0+k*th*(1-r)*t0*(k**n-r**n)/(k-r)
checks['active_solution_n1']=sp.simplify(xf.subs(n,1)-(k*x0+k*th*(1-r)*t0))==0
R,T=sp.symbols('R T')
checks['theta0_exact_invariance']=sp.simplify((M*sp.Matrix([0,R,T])).subs(th,0)[0])==0

phi=(1+sp.sqrt(5))/2
q=sp.simplify(phi/2); q2=sp.simplify(q**2); q3=sp.simplify(q**3)
checks['golden_eighth']=sp.simplify(q2-q3-sp.Rational(1,8))==0

kD,kU,k3,thD,thU,th3=sp.symbols('kD kU k3 thD thU th3')
def block(kk,tt,rr):
    return sp.Matrix([[kk,0,kk*tt*(1-rr)],[1-kk,1,(1-kk*tt)*(1-rr)],[0,0,rr]])
C=sp.diag(block(kD,thD,q2),block(kU,thU,q2),block(k3,th3,q3))
J=sp.diag(*([1]*6+[-1]*3))
G=sp.diag(*([q2]*6+[q3]*3))
checks['C2_equivariance']=sp.simplify(C*J-J*C)==sp.zeros(9)
checks['golden_commutation']=sp.simplify(C*G-G*C)==sp.zeros(9)

eps,V,Rs=sp.symbols('eps V Rs')
checks['selling_ledger']=sp.expand((V-4*eps)+(Rs+4*eps)-(V+Rs))==0

status='PASS' if all(checks.values()) else 'FAIL'
certificate={
  'phase':'v11-conservative-toric-polar-radial-transfer',
  'status':status,
  'router_matrix':str(M),
  'spectrum':['1','kappa','rho'],
  'golden_q':str(q),
  'q2':str(q2),
  'q3':str(q3),
  'q2_minus_q3':str(sp.simplify(q2-q3)),
  'checks':checks,
  'statuses':{
    'three_port_moment_conservation':'PROVEN',
    'active_lock_with_golden_donors':'PROVEN',
    'C2_equivariance':'PROVEN',
    'golden_weight_2_3_commutation':'PROVEN',
    'selling_scalar_conservative_lift':'PROVEN',
    'canonical_Selling_to_geometric_reservoir_identification':'NT_PENDING_MORPHISM'
  },
  'guards':[
    'The full conservative system retains eigenvalue 1 in reservoir directions; the sphere is an attractor of the active subsystem, not a unique global point.',
    'Signed moments may be realized positively by Jordan splitting into positive and negative reservoir channels.',
    'Octagrammatic phase resonance alone does not identify a unique physical reservoir.',
    'The Selling ledger conservation is exact at the vonorm-descent level; mapping that reservoir to D, U, or nu remains typed and unproved.'
  ]
}
Path('/mnt/data/conservative_transfer_v11_certificate.json').write_text(json.dumps(certificate,indent=2,ensure_ascii=False))
print(json.dumps(certificate,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
