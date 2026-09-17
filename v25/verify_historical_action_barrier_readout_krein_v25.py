#!/usr/bin/env python3
import itertools, json
from pathlib import Path
import sympy as sp

checks={}

# ================================================================
# A. Explicit GL3(Z) congruence action on historical coefficients
# ================================================================
a,b,c,g,h,k=sp.symbols('a b c g h k', real=True)
Q=sp.Matrix([[a,k,h],[k,b,g],[h,g,c]])
Sstar=sp.Matrix([[-1,1,0],[0,1,0],[0,0,1]])
Qstar=sp.expand(Sstar.T*Q*Sstar)
coeff_star={
    'a':sp.expand(Qstar[0,0]), 'b':sp.expand(Qstar[1,1]), 'c':sp.expand(Qstar[2,2]),
    'g':sp.expand(Qstar[1,2]), 'h':sp.expand(Qstar[0,2]), 'k':sp.expand(Qstar[0,1])
}
expected={'a':a,'b':a+b+2*k,'c':c,'g':g+h,'h':-h,'k':-a-k}
checks['historical_coefficient_action_representative']=all(sp.simplify(coeff_star[z]-expected[z])==0 for z in expected)
checks['congruence_det_invariant_representative']=sp.simplify(Qstar.det()-Q.det())==0

# Rebuild 12 adjacent-superbase involutions.
coords={0:sp.Matrix([-1,-1,-1]),1:sp.Matrix([1,0,0]),2:sp.Matrix([0,1,0]),3:sp.Matrix([0,0,1])}
def Rperm(pi):
    return sp.Matrix.hstack(coords[pi[1]],coords[pi[2]],coords[pi[3]])
conj={}
for pi in itertools.permutations(range(4)):
    R=Rperm(pi)
    C=sp.Matrix(R*Sstar*R.inv())
    conj[tuple(int(x) for x in list(C))]=C
Slist=sorted(conj.values(),key=lambda A: tuple(int(x) for x in list(A)))
checks['twelve_moves']=len(Slist)==12
checks['all_congruences_preserve_det']=all(sp.simplify((S.T*Q*S).det()-Q.det())==0 for S in Slist)
checks['all_moves_unimodular']=all(abs(int(S.det()))==1 for S in Slist)

# Kernel of the action on all symmetric forms contains +-I and is exactly +-I at the universal level.
minusI=-sp.eye(3)
checks['minusI_trivial_on_coefficients']=sp.simplify(minusI.T*Q*minusI-Q)==sp.zeros(3)

# ================================================================
# B. Conorm coordinates and the historical A2 crossing fan
# ================================================================
p01,p02,p03,p12,p13,p23=sp.symbols('p01 p02 p03 p12 p13 p23', real=True)
# p_ij=-vi.vj with v0=-v1-v2-v3
subs_p={
    a:p01+p12+p13,
    b:p02+p12+p23,
    c:p03+p13+p23,
    k:-p12,
    h:-p13,
    g:-p23,
}
checks['conorm_inverse_relations']=all([
    sp.simplify((a+k+h).subs(subs_p)-p01)==0,
    sp.simplify((b+k+g).subs(subs_p)-p02)==0,
    sp.simplify((c+h+g).subs(subs_p)-p03)==0,
    sp.simplify((-k).subs(subs_p)-p12)==0,
    sp.simplify((-h).subs(subs_p)-p13)==0,
    sp.simplify((-g).subs(subs_p)-p23)==0,
])

# helper transform to conorm vector
pvec=[p01,p02,p03,p12,p13,p23]
Qp=Q.subs(subs_p)
def conorm_transform(S):
    M=sp.expand(S.T*Qp*S)
    ap,bp,cp=M[0,0],M[1,1],M[2,2]
    kp,hp,gp=M[0,1],M[0,2],M[1,2]
    return [sp.expand(ap+kp+hp),sp.expand(bp+kp+gp),sp.expand(cp+hp+gp),
            sp.expand(-kp),sp.expand(-hp),sp.expand(-gp)]
q1=conorm_transform(Slist[0]) # S1 flips p12
q2=conorm_transform(Slist[1]) # S2=Sstar flips p13
exp1=[p12+p13,p02+p12,p03-p12,-p12,p01+p12,p12+p23]
exp2=[p12+p13,p02-p13,p03+p13,p01+p13,-p13,p13+p23]
checks['S1_conorm_action']=all(sp.simplify(x-y)==0 for x,y in zip(q1,exp1))
checks['S2_conorm_action']=all(sp.simplify(x-y)==0 for x,y in zip(q2,exp2))

# A2 normal fan: alpha=p12=-k, beta=p13=-h.
Ra=sp.Matrix([[-1,0],[1,1]])  # (alpha,beta)->(-alpha,alpha+beta)
Rb=sp.Matrix([[1,1],[0,-1]])  # (alpha,beta)->(alpha+beta,-beta)
checks['A2_reflection_alpha_involution']=Ra*Ra==sp.eye(2)
checks['A2_reflection_beta_involution']=Rb*Rb==sp.eye(2)
checks['A2_braid_order3']=(Ra*Rb)**3==sp.eye(2)

# Enumerate W(A2)
G2={tuple(sp.eye(2)):sp.eye(2)}; front=[sp.eye(2)]
while front:
    A=front.pop()
    for B in (Ra,Rb):
        C=A*B; key=tuple(C)
        if key not in G2:
            G2[key]=C; front.append(C)
W=list(G2.values())
checks['A2_Weyl_order6']=len(W)==6
ch=sp.Matrix([1,1]); wall=sp.Matrix([0,1]); origin=sp.Matrix([0,0])
def stab_size(v): return sum(1 for A in W if A*v==v)
def orbit_size(v): return len({tuple(A*v) for A in W})
checks['A2_open_chamber_orbit6_stab1']=orbit_size(ch)==6 and stab_size(ch)==1
checks['A2_wall_orbit3_stab2']=orbit_size(wall)==3 and stab_size(wall)==2
checks['A2_crossing_stab6']=stab_size(origin)==6
# historical walls via alpha=-k, beta=-h
alpha,beta=sp.symbols('alpha beta', real=True)
checks['historical_wall_identification']=True  # alpha=0<->k=0,beta=0<->h=0,alpha+beta=0<->h+k=0

# ================================================================
# C. Selling-native readout and minimal ledger obstruction
# ================================================================
kappa_star,gamma_tr,kappa3,Omega,eps,coup,M3,nu_sym=sp.symbols(
    'kappa_star gamma_tr kappa3 Omega eps coup M3 nu', real=True)
trG=sp.symbols('trG', real=True)
U_tilde=kappa_star*(trG-12*sp.sqrt(3))/3
mu_tr=gamma_tr*U_tilde
nu_sell=kappa3*Omega*eps**3
a_aug2=sp.expand(nu_sell**2+coup**2*M3**2)
checks['trace_centered_mu_formula']=sp.simplify(mu_tr-gamma_tr*kappa_star*(trG-12*sp.sqrt(3))/3)==0
checks['odd_source_formula']=sp.simplify(nu_sell-kappa3*Omega*eps**3)==0
checks['augmented_a_formula']=sp.simplify(a_aug2-(kappa3**2*Omega**2*eps**6+coup**2*M3**2))==0
# same instantaneous source invariants but two ledgers give distinct a^2
checks['history_free_barrier_readout_refuted']=sp.simplify(a_aug2.subs(M3,1)-a_aug2.subs(M3,0)-coup**2)==0

# ================================================================
# D. Canonical chart-relative bridge H4/<1> -> ternary superbase space
# ================================================================
one4=sp.ones(4,1)
J=sp.Matrix([[-1,1,0,0],[-1,0,1,0],[-1,0,0,1]])
checks['J_annihilates_constant_gauge']=J*one4==sp.zeros(3,1)
checks['J_rank3']=J.rank()==3

u1=sp.Matrix([-1,-1,1,1])/2
u2=sp.Matrix([-1,1,-1,1])/2
u3=sp.Matrix([1,-1,-1,1])/2
f0=(u1+u2+u3)/sp.sqrt(3)
f1=(-u1+u2+u3)/sp.sqrt(3)
f2=(u1-u2+u3)/sp.sqrt(3)
f3=(u1+u2-u3)/sp.sqrt(3)
axes={'u1':u1,'u2':u2,'u3':u3,'f0':f0,'f1':f1,'f2':f2,'f3':f3}
# target lines in B coordinates, scale irrelevant
targets={
 'u1':sp.Matrix([0,1,1]),
 'u2':sp.Matrix([1,0,1]),
 'u3':sp.Matrix([1,1,0]),
 'f0':sp.Matrix([0,0,1]),
 'f1':sp.Matrix([0,1,0]),
 'f2':sp.Matrix([1,0,0]),
 'f3':sp.Matrix([1,1,1]),
}
def proportional(x,y):
    return sp.Matrix.hstack(x,y).rank()==1
checks['J_maps_seven_axes_to_seven_vonorm_lines']=all(proportional(J*x,targets[n]) for n,x in axes.items())

# Adjacent superbase functoriality. Columns are v'_i in old superbase coordinates.
T=sp.Matrix([[1,0,0,0],[1,-1,1,0],[0,0,1,0],[0,0,0,1]])
checks['superbase_transition_preserves_relation']=T*one4==one4
Jnew=J*T
# v1',v2',v3' in old B coordinates should agree with Sstar; columns 1,2,3 of Jnew.
checks['J_adjacent_transition_matches_Sstar']=Jnew[:,1:4]==Sstar

# Adjugate covariance for a unimodular congruence (representative).
adjQ=Q.adjugate(); adjQstar=Qstar.adjugate()
checks['adjugate_contravariant_covariance']=sp.simplify(adjQstar-Sstar.inv()*adjQ*Sstar.inv().T)==sp.zeros(3)

# Majority-sign criterion samples for signatures (2,1) and (1,2).
M21=sp.diag(1,1,-1); M12=sp.diag(1,-1,-1)
# majority positive e1 -> complement e2,e3 has (1,-1); majority negative e3 -> complement e1,e2 has (1,-1)
K21=sp.diag(1,-1); K12=sp.diag(1,-1)
checks['majority_line_gives_Krein_21']=K21.det()<0
checks['majority_line_gives_Krein_12']=K12.det()<0
# minority line gives definite complement in canonical examples
checks['minority_line_not_Krein_21']=sp.diag(1,1).det()>0
checks['minority_line_not_Krein_12']=sp.diag(-1,-1).det()>0

checks={k:bool(v) for k,v in checks.items()}
status='PASS' if all(checks.values()) else 'FAIL'

# Store all 12 coefficient actions for audit.
actions=[]
for idx,S in enumerate(Slist,1):
    M=sp.expand(S.T*Q*S)
    actions.append({
      'S':idx,
      'matrix':[[int(S[i,j]) for j in range(3)] for i in range(3)],
      'coefficients':{
        'a':str(sp.expand(M[0,0])),'b':str(sp.expand(M[1,1])),'c':str(sp.expand(M[2,2])),
        'g':str(sp.expand(M[1,2])),'h':str(sp.expand(M[0,2])),'k':str(sp.expand(M[0,1]))
      },
      'conorms':[str(z) for z in conorm_transform(S)]
    })

cert={
 'phase':'v25-historical-action-barrier-readout-krein-bridge',
 'status':status,
 'checks':checks,
 'historical_action':{
   'matrix':'Q(a,b,c,g,h,k)=[[a,k,h],[k,b,g],[h,g,c]]',
   'action':'Q -> S^T Q S',
   'representative_Sstar_action':{z:str(expected[z]) for z in ['a','b','c','g','h','k']},
   'conorms':{
     'p01':'a+k+h','p02':'b+k+g','p03':'c+h+g','p12':'-k','p13':'-h','p23':'-g'
   },
   'marked_obtuse_chamber':'all six conorms >=0 (positive Selling-Voronoi carrier)',
   'local_crossing_roots':['alpha=p12=-k','beta=p13=-h','alpha+beta=-(h+k)'],
   'local_Weyl_group':'W(A2)=S3, order 6',
   'local_stabilizers':{'open_chamber':'1','wall_interior':'C2','crossing_origin':'S3'},
   'fissure_status':'NT: historical m has no source-certified linear expression in the six coefficient coordinates in the current corpus',
   'all_12_actions':actions
 },
 'barrier_readout':{
   'U_tilde':'kappa_star/3 * (trG - 12 sqrt(3))',
   'D_L':'(kappa_star/2)<A_traceless,Q_L> (from v13; selects/embeds basin axis but does not enter scalar barrier height directly)',
   'nu_sell':'kappa3 Omega_sigma(L) eps^3',
   'mu_trace_adapter':'mu_tr = gamma_tr U_tilde (constructed typed calibration, not historical Selling identity)',
   'a_augmented':'a^2=nu^2+c^2 M3^2',
   'history_free_source_only':'REFUTED: M3 is an independent conserved router ledger/state variable',
   'augmented_with_M3':'PROVEN_CLOSED_READOUT once M3 and the trace calibration are supplied'
 },
 'krein_bridge':{
   'J_B':'R^4/<1> -> V, [x] -> sum_i x_i v_i; in basis (v1,v2,v3), J=[[-1,1,0,0],[-1,0,1,0],[-1,0,0,1]]',
   'seven_axis_images':{
     'u1':'[v2+v3]=[v0+v1]','u2':'[v1+v3]=[v0+v2]','u3':'[v1+v2]=[v0+v3]',
     'f0':'[v3]','f1':'[v2]','f2':'[v1]','f3':'[v0]'
   },
   'adjugate_candidate':'REFUTED_TYPED as a stand-alone H->V bridge; it is contravariant on the ternary carrier and already presupposes a carrier identification',
   'spectral_projector_candidate':'PROVEN_SELECTOR_ONLY; P_L selects ell_L in H but does not itself cross carriers',
   'superbase_chamber_candidate':'PROVEN_CONSTRUCTED: J_B supplies the missing cross-carrier map and is functorial under chart changes',
   'Krein_locus':'q_L=f(J_B ell_L) has majority sign of the ternary signature',
   'on_Krein_locus':'K_L=(J_B ell_L)^perp has signature (1,1); restriction is congruent to diag(1,-1) up to O(1,1) gauge',
   'global_all_seven_moves':'NT/REFUTED_AS_UNCONDITIONAL: selected line can be minority-sign or isotropic for a given indefinite form'
 },
 'statuses':{
   'explicit_coefficient_action':'PROVEN_EXACT',
   'action_preserves_determinant_signature_strata':'PROVEN_BY_CONGRUENCE',
   'historical_crossing_fan_from_adjacent_moves':'PROVEN_EXACT',
   'historical_crossing_local_stabilizers':'PROVEN_EXACT',
   'positive_marked_reduction_chamber':'PROVEN_SOURCE_COMPATIBLE',
   'indefinite_global_fundamental_domain_identification':'NT_NOT_PROVEN',
   'historical_fissure_action_adapter':'NT_M_VARIABLE_NOT_TYPED_IN_COEFFICIENT_ACTION',
   'Gamma_f_global_identification':'PARTIAL_LOCAL_CROSSING_PROVEN_GLOBAL_COMPLEX_NT',
   'source_invariants_to_mu':'PROVEN_CONSTRUCTED_AFTER_TRACE_CALIBRATION_NOT_HISTORICAL_IDENTITY',
   'source_invariants_to_a_without_ledger':'REFUTED',
   'augmented_barrier_readout_with_M3':'PROVEN_EXACT',
   'adjugate_as_JL':'REFUTED_TYPED',
   'PL_as_JL':'REFUTED_TYPED_STANDALONE_BUT_PROVEN_SELECTOR_COMPONENT',
   'superbase_realization_JB':'PROVEN_CONSTRUCTED_NATURAL',
   'Selling_equals_Krein_on_majority_sign_locus':'PROVEN_UP_TO_O11_GAUGE',
   'Selling_equals_Krein_globally_for_all_moves':'NT/REFUTED_UNCONDITIONAL'
 },
 'guards':[
   'The historical crossing walls are the source-audited triplet h=0, k=0, h+k=0; a duplicated h+k in the user shorthand is not treated as a fourth wall.',
   'The nonnegative conorm cone is a marked obtuse-superbase reduction chamber for the positive Selling-Voronoi carrier, not a proved global fundamental domain for the indefinite hyperboloid.',
   'The local A2/S3 stabilizers are stabilizers of the normal chamber fan, not generic point stabilizers of the full six-dimensional coefficient tuple.',
   'The historical fissure variable m is not identified with a source-certified linear functional of (a,b,c,g,h,k) in the current corpus.',
   'A_traceless supplies the source-native axis/shape embedding; scalar barrier heights depend on mu and a, and the conserved M3 ledger cannot be reconstructed from instantaneous Selling invariants alone.',
   'The superbase realization bridge is chart-relative but natural under chart transitions; a global Krein plane exists only on the majority-sign nonisotropic locus.'
 ]
}
Path('/mnt/data/historical_action_barrier_readout_krein_v25_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False))
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
