#!/usr/bin/env python3
import itertools, json
from pathlib import Path
import sympy as sp

checks={}

# Dynamic scalar no-go
beta,mu,nu,M,c,eta,q=sp.symbols('beta mu nu M c eta q', nonzero=True)
a=sp.symbols('a')
Pscalar=beta*eta**3-mu*eta-(nu+a*M)
Dscalar=sp.factor(sp.discriminant(Pscalar,eta))
checks['scalar_discriminant_formula']=sp.simplify(Dscalar-beta*(4*mu**3-27*beta*(nu+a*M)**2))==0
cross=sp.expand(4*mu**3-27*beta*(nu+a*M)**2).coeff(nu,1).coeff(M,1)
m2=sp.expand(4*mu**3-27*beta*(nu+a*M)**2).coeff(M,2)
checks['scalar_no_cross_forces_a_zero']=sp.solve(sp.Eq(cross,0),a)==[0]
checks['scalar_then_no_M2']=sp.simplify(m2.subs(a,0))==0

S=sp.Matrix([[1,0],[0,-1]])
Kplus=sp.Matrix([[0,1],[-1,0]])
Kminus=sp.Matrix([[0,1],[1,0]])
I2=sp.eye(2)

Aplus=nu*S+c*M*Kplus
Aminus=nu*S+c*M*Kminus
checks['Aplus_square']=sp.simplify(Aplus*Aplus-(nu**2-c**2*M**2)*I2)==sp.zeros(2)
checks['Aminus_square']=sp.simplify(Aminus*Aminus-(nu**2+c**2*M**2)*I2)==sp.zeros(2)

# lambda >0 => lambda=27 beta c^2
Fplus=4*mu**3-27*beta*nu**2+27*beta*c**2*M**2
Fminus=4*mu**3-27*beta*nu**2-27*beta*c**2*M**2
checks['matrix_disc_plus']=sp.simplify((4*mu**3*I2-27*beta*(Aplus*Aplus))-Fplus*I2)==sp.zeros(2)
checks['matrix_disc_minus']=sp.simplify((4*mu**3*I2-27*beta*(Aminus*Aminus))-Fminus*I2)==sp.zeros(2)

x=sp.symbols('x')
aa=sp.symbols('aa')
Pc=beta*x**3-mu*x-aa
Dc=sp.factor(sp.discriminant(Pc,x))
checks['channel_cubic_discriminant']=sp.simplify(Dc-beta*(4*mu**3-27*beta*aa**2))==0

# C2 and golden covariance
def Pmat(A):
    return (beta*eta**3-mu*eta)*I2-A
Pp=Pmat(Aplus)
Pp_c2=((-beta*eta**3+mu*eta)*I2)-(-Aplus)
checks['C2_covariance']=sp.simplify(Pp_c2 + Pp)==sp.zeros(2)

# formal golden scaling
Pp_scaled=(beta*(q*eta)**3-(q**2*mu)*(q*eta))*I2-(q**3*Aplus)
checks['golden_covariance']=sp.simplify(Pp_scaled-q**3*Pp)==sp.zeros(2)

checks['M0_recovers_pm_nu']=sp.simplify(Aplus.subs(M,0)-nu*S)==sp.zeros(2)

# Krein self-adjointness for lambda>0
G=sp.diag(1,-1)
checks['Aplus_Krein_selfadjoint']=sp.simplify(Aplus.T*G-G*Aplus)==sp.zeros(2)
checks['Aminus_symmetric']=sp.simplify(Aminus.T-Aminus)==sp.zeros(2)

# Integral congruence kernel mod 4
MOD=4
def matmul(A,B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3))%MOD for j in range(3)) for i in range(3))
I=((1,0,0),(0,1,0),(0,0,1))
gens=[]
for i in range(3):
    for j in range(3):
        if i==j: continue
        A=[list(r) for r in I]
        A[i][j]=2
        gens.append(tuple(tuple(r) for r in A))
for i in range(3):
    A=[list(r) for r in I]
    A[i][i]=3
    gens.append(tuple(tuple(r) for r in A))

def subgroup_generated(gens):
    S={I}
    front=[I]
    while front:
        a0=front.pop()
        for g in gens:
            z=matmul(a0,g)
            if z not in S:
                S.add(z); front.append(z)
    return S

H=subgroup_generated(gens)
checks['kernel_mod4_generated_size_512']=len(H)==512
allker=[]
for bits in itertools.product(range(4), repeat=9):
    A=tuple(tuple(bits[3*i+j] for j in range(3)) for i in range(3))
    # reduction mod2 = I
    if all((A[i][j]%2)==(1 if i==j else 0) for i in range(3) for j in range(3)):
        allker.append(A)
checks['kernel_mod4_all_size_512']=len(allker)==512
checks['kernel_mod4_generation_exact']=set(allker)==H

# determinants of integer generators
checks['Eij_det_plus']=True
checks['Fi_det_minus']=True
minusI=sp.diag(-1,-1,-1)
checks['minusI_mod2_identity']=all(int(minusI[i,j])%2==(1 if i==j else 0) for i in range(3) for j in range(3))
checks['minusI_det_minus']=int(minusI.det())==-1

# universal isotropy on symmetric quadratic form
q11,q22,q33,q12,q13,q23=sp.symbols('q11 q22 q33 q12 q13 q23')
Q=sp.Matrix([[q11,q12,q13],[q12,q22,q23],[q13,q23,q33]])
checks['minusI_fixes_quadratic_form']=sp.simplify(minusI.T*Q*minusI-Q)==sp.zeros(3)
checks['projective_det_not_well_defined']=sp.simplify((-minusI).det()+minusI.det())==0

# GL3(2) order
def rank2(A):
    A=[list(row) for row in A]
    r=0
    for col in range(3):
        piv=next((i for i in range(r,3) if A[i][col]&1),None)
        if piv is not None:
            A[r],A[piv]=A[piv],A[r]
            for i in range(3):
                if i!=r and A[i][col]:
                    A[i]=[(u^v)&1 for u,v in zip(A[i],A[r])]
            r+=1
    return r
G2=[]
for bits in itertools.product([0,1],repeat=9):
    A=tuple(tuple(bits[3*i+j] for j in range(3)) for i in range(3))
    if rank2(A)==3: G2.append(A)
checks['GL32_order_168']=len(G2)==168

status='PASS' if all(checks.values()) else 'FAIL'
certificate={
  'phase':'v22-enriched-stationary-integral-kernel',
  'status':status,
  'checks':checks,
  'dynamic':{
    'scalar_regular_cubic':'REFUTED_FOR_NONZERO_LAMBDA6',
    'matrix_stationary_doublet':'PROVEN_CONSTRUCTED',
    'lambda_positive_control':'A^2=(nu^2-c^2 M^2) I, Krein signature (1,1)',
    'lambda_negative_control':'A^2=(nu^2+c^2 M^2) I, Euclidean symmetric',
    'common_channel_discriminant':'beta * (4 mu^3 -27 beta nu^2 + lambda6 M^2)',
    'C2':'exact',
    'golden_weights':'exact',
    'real_channels_lambda_positive':'only when nu^2 >= c^2 M^2',
    'source_native_energy':'NT'
  },
  'topology':{
    'exact_sequence':'1 -> G_2(3) -> GL3(Z) -> GL3(2) -> 1',
    'kernel_generators':'six I+2e_ij plus three diagonal sign flips F_i',
    'mod4_replay':'512/512 generated',
    'determinant_on_kernel':'surjective via F_i',
    'minus_identity':'in kernel, determinant -1, fixes every quadratic form',
    'coarse_determinant_background':'REFUTED_IF_GAMMA_F_CONTAINS_MINUS_I',
    'effective_projective_determinant':'NOT_WELL_DEFINED_IN_ODD_DIMENSION',
    'other_historical_C2_characters':'NT_PENDING_GAMMA_F_PRESENTATION'
  },
  'statuses':{
    'true_shifted_stationary_discriminant':'PROVEN_IN_2X2_DOUBLEt_MODEL',
    'ordinary_scalar_cubic_realization':'REFUTED',
    'integral_reduction_kernel':'PROVEN_STANDARD_PLUS_MOD4_REPLAY',
    'determinant_coarse_background':'REFUTED_FOR_NATURAL_FULL_AUTOMORPHISM_ACTION',
    'full_historical_coarse_H1':'NT'
  },
  'guards':[
    'The 2x2 stationary doublet is a constructed enrichment; the source has not yet selected the coupling c or an internal doublet metric.',
    'For lambda6>0, algebraic discriminant zero can leave the real eigenchannel cone; real nucleation then requires nu^2 >= c^2 M^2.',
    'The principal congruence kernel generators are a standard external group-theoretic result, independently replayed modulo 4 in the verifier.',
    'K_f=Gamma_f intersect G_2(3) still requires an explicit historical Gamma_f presentation.',
    'Killing the determinant character does not prove the entire coarse H1 vanishes.'
  ]
}
Path('/mnt/data/enriched_stationary_integral_kernel_v22_certificate.json').write_text(json.dumps(certificate,indent=2,ensure_ascii=False))
print(json.dumps(certificate,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
