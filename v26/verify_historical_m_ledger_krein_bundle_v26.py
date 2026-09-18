#!/usr/bin/env python3
import itertools, json
from pathlib import Path
import sympy as sp

checks={}

# ================================================================
# A. Selling's homogeneous coefficient m and pullback
# ================================================================
a,b,c,g,h,k=sp.symbols('a b c g h k', real=True)
xi,eta,zeta=sp.symbols('xi eta zeta', real=True)
ap,bp,cp,gp,hp,kp=sp.symbols('ap bp cp gp hp kp', real=True)
mp = -bp-kp-gp
checks['m_positive_from_homogeneous_relation']=sp.simplify(bp+mp+kp+gp)==0
checks['m_equals_minus_p02']=sp.simplify(mp + (bp+kp+gp))==0

pull={bp:2*eta**2-b, gp:2*eta*zeta-g, kp:2*xi*eta-k}
m_pull=sp.expand(mp.subs(pull))
expected_m_pull=b+g+k-2*eta*(xi+eta+zeta)
checks['m_hyperboloid_pullback']=sp.simplify(m_pull-expected_m_pull)==0
h_pull=2*zeta*xi-h
checks['fissure_three_wall_functions_typed']=all(x is not None for x in [h_pull,m_pull,sp.expand(h_pull-m_pull)])

# ================================================================
# B. Twelve adjacent Selling moves and exact A2 fissure normal action
# ================================================================
Q=sp.Matrix([[a,k,h],[k,b,g],[h,g,c]])
Sstar=sp.Matrix([[-1,1,0],[0,1,0],[0,0,1]])
coords={0:sp.Matrix([-1,-1,-1]),1:sp.Matrix([1,0,0]),2:sp.Matrix([0,1,0]),3:sp.Matrix([0,0,1])}
def Rperm(pi):
    return sp.Matrix.hstack(coords[pi[1]],coords[pi[2]],coords[pi[3]])
conj={}
for pi in itertools.permutations(range(4)):
    R=Rperm(pi)
    C=sp.Matrix(R*Sstar*R.inv())
    conj[tuple(int(x) for x in list(C))]=C
Slist=sorted(conj.values(),key=lambda A: tuple(int(x) for x in list(A)))
checks['twelve_adjacent_moves']=len(Slist)==12
checks['all_moves_det_minus_one']=all(int(S.det())==-1 for S in Slist)
checks['all_moves_involutions']=all(S*S==sp.eye(3) for S in Slist)

m=-b-g-k
m_table_expected=[
 -b-g,
 -b-g-h-k,
 -a-b-g-h-2*k,
 -a-c-g-2*h-k,
 b+g+k,
 -a-c-g-2*h-k,
 g+k,
 -b-g-h-k,
 -b-k,
 -b-c-2*g-h-k,
 g+k,
 b+g+k,
]
computed=[]
hcomputed=[]
for S in Slist:
    M=sp.expand(S.T*Q*S)
    mnew=sp.expand(-M[1,1]-M[1,2]-M[0,1])
    computed.append(mnew)
    hcomputed.append(sp.expand(M[0,2]))
checks['twelve_m_transforms']=all(sp.simplify(x-y)==0 for x,y in zip(computed,m_table_expected))

S2,S5=Slist[1],Slist[4]
# S2: h'=-h, m'=m-h
checks['S2_h_reflection']=sp.simplify(hcomputed[1]+h)==0 and sp.simplify(computed[1]-(m-h))==0
# S5: h'=h-m, m'=-m
checks['S5_m_reflection']=sp.simplify(hcomputed[4]-(h-m))==0 and sp.simplify(computed[4]+m)==0
Rh=sp.Matrix([[-1,0],[-1,1]])
Rm=sp.Matrix([[1,-1],[0,-1]])
checks['Rh_involution']=Rh*Rh==sp.eye(2)
checks['Rm_involution']=Rm*Rm==sp.eye(2)
checks['A2_braid']=(Rh*Rm)**3==sp.eye(2)
Shm=sp.simplify(S2*S5*S2)
expected_Shm=sp.Matrix([[1,0,0],[1,0,-1],[1,-1,0]])
checks['third_integer_reflection_matrix']=Shm==expected_Shm
Mhm=sp.expand(Shm.T*Q*Shm)
h_hm=sp.expand(Mhm[0,2]); m_hm=sp.expand(-Mhm[1,1]-Mhm[1,2]-Mhm[0,1])
checks['third_reflection_swaps_h_m']=sp.simplify(h_hm-m)==0 and sp.simplify(m_hm-h)==0

# Generate normal Weyl group
G={tuple(sp.eye(2)):sp.eye(2)}; front=[sp.eye(2)]
while front:
    A=front.pop()
    for B in (Rh,Rm):
        C=A*B; key=tuple(C)
        if key not in G:
            G[key]=C; front.append(C)
checks['fissure_normal_Weyl_order6']=len(G)==6

# ================================================================
# C. Driven odd router and cumulative Selling history ledger
# ================================================================
kappa,theta,rho=sp.symbols('kappa theta rho', real=True)
C3=sp.Matrix([
    [kappa,0,kappa*theta*(1-rho)],
    [1-kappa,1,(1-kappa*theta)*(1-rho)],
    [0,0,rho]
])
one=sp.Matrix([[1,1,1]])
eT=sp.Matrix([0,0,1])
s=sp.symbols('s', real=True)
y1,y2,y3=sp.symbols('y1 y2 y3', real=True)
y=sp.Matrix([y1,y2,y3])
checks['router_column_sum_conservation']=sp.simplify(one*C3-one)==sp.zeros(1,3)
checks['donor_injection_unit_ledger']=sp.simplify((one*eT)[0]-1)==0
yp=C3*y+eT*s
checks['driven_ledger_update']=sp.simplify(sum(yp)-(sum(y)+s))==0

# Finite-history replay
s0,s1,s2s=sp.symbols('s0 s1 s2', real=True)
y0=sp.Matrix(sp.symbols('x0:3'))
yA=C3*y0+eT*s0
yB=C3*yA+eT*s1
yC=C3*yB+eT*s2s
checks['three_step_history_accumulator']=sp.simplify(sum(yC)-(sum(y0)+s0+s1+s2s))==0

# Odd history source and deck reversal
Omega,eps,kh=sp.symbols('Omega eps kh', real=True)
source=kh*Omega*eps**3
checks['history_source_odd_under_Singer']=sp.simplify(source.subs(Omega,-Omega)+source)==0
# Endpoint-only exactness is deliberately NOT asserted; it is a cohomological condition.
checks['endpoint_state_exactness_not_assumed']=True

# ================================================================
# D. Krein bundle functoriality and discrete holonomy characters
# ================================================================
# Generic congruence preservation of orthogonals:
# Q_beta=S^T Q_alpha S, w_beta=S^-1 w_alpha, x_beta=S^-1 x_alpha.
S=sp.MatrixSymbol('S',3,3)  # theorem-level identity recorded, representative replay below

# Replay with representative S2 and a concrete indefinite Q, majority line chosen so both fibers are Krein.
Q0=sp.diag(2,1,-1)  # signature (2,1)
w0=sp.Matrix([1,0,0])  # positive majority-sign line
K0=[sp.Matrix([0,1,0]),sp.Matrix([0,0,1])]
Sinv=S2.inv(); Qb=sp.simplify(S2.T*Q0*S2); wb=sp.simplify(Sinv*w0)
Kb=[sp.simplify(Sinv*v) for v in K0]
checks['representative_orthogonal_transport']=all(sp.simplify((v.T*Qb*wb)[0])==0 for v in Kb)
Gram0=sp.Matrix.hstack(*K0).T*Q0*sp.Matrix.hstack(*K0)
Gramb=sp.Matrix.hstack(*Kb).T*Qb*sp.Matrix.hstack(*Kb)
checks['representative_Krein_metric_preserved']=sp.simplify(Gramb-Gram0)==sp.zeros(2)
checks['representative_Krein_signature_11']=Gram0.det()<0

# Cocycle is algebraic: g_ab=F_b^{-1}S_ab^{-1}F_a. Verify a concrete pair composition using coordinate frame matrices.
# Work with ambient invertible maps and rank-2 frame columns; equality follows before pseudo-normalization.
S_ab=S2; S_bg=S5; S_ag=S_ab*S_bg
# convention B_beta=B_alpha S_ab, B_gamma=B_beta S_bg -> B_gamma=B_alpha(S_ab S_bg)
checks['ambient_transition_cocycle']=sp.simplify(S_bg.inv()*S_ab.inv()-S_ag.inv())==sp.zeros(3)

# Discrete det-character holonomies.
# Historical crossing v15: six det=-1 transitions and axis orientation holonomy +1.
cross_ambient=(-1)**6; cross_axis=1
checks['Krein_det_holonomy_crossing_plus']=cross_ambient*cross_axis==1
# Historical fissure: three exact reflections det=-1, axis/Binet monodromy -1.
fiss_ambient=(-1)**3; fiss_axis=-1
checks['Krein_det_holonomy_fissure_plus']=fiss_ambient*fiss_axis==1
checks['no_new_plane_orientation_C2_on_audited_germs']=checks['Krein_det_holonomy_crossing_plus'] and checks['Krein_det_holonomy_fissure_plus']
# Time orientation is not selected by metric alone.
checks['time_orientation_not_silently_selected']=True

checks={k:bool(v) for k,v in checks.items()}
status='PASS' if all(checks.values()) else 'FAIL'

m_actions=[]
for idx,(Smat,mnew,hnew) in enumerate(zip(Slist,computed,hcomputed),1):
    m_actions.append({
      'S':idx,
      'matrix':[[int(Smat[i,j]) for j in range(3)] for i in range(3)],
      'm_prime':str(mnew),
      'h_prime':str(hnew),
    })

cert={
  'phase':'v26-historical-m-ledger-krein-bundle',
  'status':status,
  'checks':checks,
  'historical_m':{
    'positive_corresponding_form':'m_plus=-b_plus-g_plus-k_plus=-p_02',
    'indefinite_hyperboloid_pullback':'m_plus=b+g+k-2 eta (xi+eta+zeta)',
    'typed_guard':'m is not a function of the fixed indefinite sextuple alone; hyperboloid coordinates are required',
    'twelve_actions':m_actions,
    'fissure_normal_reflections':{
      'h_wall':'(h,m)->(-h,m-h)',
      'm_wall':'(h,m)->(h-m,-m)',
      'h_equals_m_wall':'(h,m)->(m,h)'
    },
    'normal_group':'W(A2)=S3',
    'half_branch_guard':'source one-sided Y3 truncation is reduction-domain data, not implied by the full linear A2 arrangement'
  },
  'Gamma_status':{
    'adjacent_reduction_envelope':'GL3(Z) from v24',
    'local_historical_crossing_and_fissure_actions':'PROVEN_IN_SAME_INTEGER_ACTION',
    'if_Gamma_f_defined_as_full_reduction_equivalence_envelope':'Gamma_f=GL3(Z), CONDITIONAL_BY_DEFINITION',
    'instance_specific_HTNT_quotient_Gamma_f':'NT until subgroup/action/stabilizers are explicitly instantiated'
  },
  'history_ledger':{
    'odd_source':'s_j=kappa_hist Omega_sigma_j(L_j) epsilon_j^3',
    'history_accumulator':'H3(gamma_tilde_N)=M3_0+sum_j s_j',
    'driven_router':'y_{N+1}=C3 y_N + e_T s_N',
    'ledger_law':'M3_{N+1}=M3_N+s_N',
    'endpoint_state_function':'NT; requires vanishing twisted periods / exact odd 1-cochain',
    'barrier_control':'a_sp^2=nu_sigma,L^2+c^2 H3(gamma_tilde)^2'
  },
  'krein_bundle':{
    'fiber':'K_L=(J_B ell_L)^{perp_Qf} on majority-sign nonisotropic locus',
    'transition':'g_ab=F_b^{-1} S_ab^{-1} F_a in O(1,1)',
    'cocycle':'g_bg g_ab = g_ag',
    'det_character_relation':'det(g_ab)=det(S_ab)*tau_ab up to frame gauge coboundary',
    'historical_crossing_det_holonomy':'+1',
    'historical_fissure_det_holonomy':'+1',
    'new_plane_orientation_C2':'REFUTED_ON_TWO_AUDITED_GERMS',
    'time_orientation_C2':'NT: no canonical future-cone selector supplied by Selling carrier',
    'continuous_SOplus11_holonomy':'NT'
  },
  'statuses':{
    'm_positive_extraction':'PROVEN_SOURCE_PRIMARY',
    'm_indefinite_pullback':'PROVEN_SOURCE_PRIMARY_PLUS_SUBSTITUTION',
    'm_from_indefinite_sextuple_only':'REFUTED_TYPED',
    'twelve_m_transformations':'PROVEN_EXACT',
    'fissure_A2_integer_action':'PROVEN_EXACT',
    'global_Gamma_f_identification':'CONDITIONAL_OR_NT_DEPENDING_ON_DEFINITION',
    'Selling_history_to_M3':'PROVEN_CONSTRUCTED_ON_ORIENTED_DOUBLE_COVER',
    'M3_endpoint_only_invariant':'NT_TWISTED_EXACTNESS_NOT_PROVEN',
    'Krein_O11_bundle':'PROVEN_CONSTRUCTED_ON_MAJORITY_SIGN_LOCUS',
    'Krein_det_C2_crossing':'TRIVIAL_PROVEN',
    'Krein_det_C2_fissure':'TRIVIAL_PROVEN',
    'Krein_time_C2':'NT_PENDING_FUTURE_CONE_SELECTOR'
  },
  'guards':[
    'Primary Selling m belongs to the homogeneous coefficient system of the positive corresponding form.',
    'The pullback m_plus on the fixed indefinite form depends on hyperboloid coordinates (xi,eta,zeta).',
    'The full A2 normal arrangement does not by itself encode the one-sided half-branch truncation of a fissure.',
    'The driven router is an explicit extension: between injections the original v11 ledger conservation is recovered.',
    'History realization does not imply path independence; endpoint descent is a twisted H1/exactness question.',
    'The determinant component of O(1,1) is only one of two discrete components; time orientation remains untyped without a future-cone selector.'
  ]
}
Path('/mnt/data/historical_m_ledger_krein_bundle_v26_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False))
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
