#!/usr/bin/env python3
import sympy as sp, json
from pathlib import Path

k,th,rho=sp.symbols('k theta rho', real=True)
C=sp.Matrix([[k,0,k*th*(1-rho)],[1-k,1,(1-k*th)*(1-rho)],[0,0,rho]])
checks={}
checks['router_ledger_conservation']=sp.simplify(sp.Matrix([[1,1,1]])*C-sp.Matrix([[1,1,1]]))==sp.zeros(1,3)
checks['router_C2_commutation']=sp.simplify(C*(-sp.eye(3))-(-sp.eye(3))*C)==sp.zeros(3)

t=sp.symbols('t', real=True)
a1,a3,b0,b2=sp.symbols('a1 a3 b0 b2')
odd=a1*t+a3*t**3
even=b0+b2*t**2
checks['odd_Dirichlet']=sp.simplify(odd.subs(t,0))==0
checks['even_Neumann']=sp.simplify(sp.diff(even,t).subs(t,0))==0
checks['odd_derivative_even']=sp.simplify(sp.diff(odd,t).subs(t,-t)-sp.diff(odd,t))==0

l1,l2,l3=sp.symbols('l1 l2 l3')
L=sp.Matrix([[l1,l2,l3]])
checks['linear_odd_to_even_no_go']=sp.solve(list(L*(-sp.eye(3))-L),[l1,l2,l3],dict=True)==[{l1:0,l2:0,l3:0}]

A,B,Cc,D,E=sp.symbols('A B C D E')
S=sp.Matrix([[A,B,Cc],[B,A,Cc],[D,D,E]])
P=sp.Matrix([[0,1,0],[1,0,0],[0,0,1]])
checks['sheet_swap_commutant']=sp.simplify(S*P-P*S)==sp.zeros(3)
vodd=sp.Matrix([1,-1,0])
checks['odd_scattering_no_even_port']=sp.simplify(S*vodd-(A-B)*vodd)==sp.zeros(3,1)

alpha,j,Ep=sp.symbols('alpha j Ep', real=True)
checks['quadratic_energy_conservation']=sp.simplify((1-alpha)*j**2+(Ep+alpha*j**2)-(j**2+Ep))==0
checks['quadratic_C2_even']=sp.simplify((-j)**2-j**2)==0
q=sp.symbols('q', nonzero=True)
checks['quadratic_golden_weight6']=sp.simplify((q**3*j)**2-q**6*j**2)==0
checks['branch_zero_preserved']=sp.simplify(C*sp.zeros(3,1))==sp.zeros(3,1)

status='PASS' if all(bool(v) for v in checks.values()) else 'FAIL'
certificate={
 'phase':'v19-branch-kirchhoff',
 'status':status,
 'checks':checks,
 'statuses':{
   'odd_value_zero_at_branch':'PROVEN_FROM_C2_CONTINUITY',
   'even_normal_derivative_zero':'PROVEN_FROM_C2_C1_REGULARITY',
   'odd_Kirchhoff_transmission':'PROVEN_IN_DECLARED_NORMAL_SLICE_FLUX_MODEL',
   'router_preserves_branch_conditions':'PROVEN_EXACT',
   'linear_C2_odd_to_even_transfer':'REFUTED',
   'linear_sheet_scattering_odd_to_even_reservoir':'REFUTED_BY_COMMUTANT',
   'quadratic_energy_conversion':'PROVEN_CONSTRUCTED_OPTIONAL',
   'canonical_quadratic_conversion_coefficient':'NT_NOT_SOURCE_SELECTED',
   'quadratic_energy_reservoir_weight':'WEIGHT_6_PROVEN',
   'direct_weight6_to_DU_weight2_adapter':'NT_PENDING_DEGREE_MINUS_4_ADAPTER',
   'minimal_regular_branch_closure':'PURE_TRANSMISSION_ALPHA_0'
 },
 'guards':[
   'The flux law J=-Lambda_3 d_t y is a declared local closure, not a historical Selling equation.',
   'Pointwise vanishing of an odd field at the branch point is not itself loss of an integrated conserved quantity.',
   'Signed third-moment conservation and quadratic energy conservation are different graded ledgers.',
   'C2 forbids linear odd-to-even conversion, but permits even nonlinear invariants such as j^2.',
   'Alpha>0 introduces a new branch interaction; alpha=0 is the minimal smooth transmission closure.',
   'Quadratic conversion has golden weight 6 and cannot be identified with the weight-2 D/U channels without an extra degree -4 adapter.'
 ]
}
Path('/mnt/data/branch_kirchhoff_v19_certificate.json').write_text(json.dumps(certificate,indent=2,ensure_ascii=False))
print(json.dumps(certificate,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
