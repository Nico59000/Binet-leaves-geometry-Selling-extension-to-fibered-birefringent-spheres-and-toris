#!/usr/bin/env python3
import json, math, cmath, csv
from pathlib import Path
from collections import deque
import sympy as sp

OUT=Path('/mnt/data/v31_work')
checks={}

# Primary historical form and source-normalized Fig. 1 substitutions.
Q0=sp.diag(12,-1,-5)
P=sp.Matrix([[1,1,0],[3,4,0],[0,0,1]])
Q=sp.Matrix([[2,0,5],[0,1,0],[3,0,8]])
R=sp.Matrix([[-4,0,5],[0,-1,0],[-3,0,4]])
S=sp.Matrix([[-2,1,0],[-3,2,0],[0,0,-1]])

# Derived source transformations.
T=sp.simplify(P*R*P.inv())
U=sp.simplify(Q*S*Q.inv())
DP=sp.diag(1,-1,1)    # differs from Selling's displayed sign matrix at most by central -I
DQ=sp.diag(1,1,-1)    # same remark; commutator is unchanged by D -> -D
V=sp.simplify(P*DP*P.inv()*DP)
W=sp.simplify(Q*DQ*Q.inv()*DQ)
Y_v30=sp.Matrix([[7,2,0],[24,7,0],[0,0,1]])

expected={
 'T':sp.Matrix([[-13,3,5],[-36,8,15],[-12,3,4]]),
 'U':sp.Matrix([[-17,2,10],[-24,2,15],[-24,3,14]]),
 'V':sp.Matrix([[7,2,0],[24,7,0],[0,0,1]]),
 'W':sp.Matrix([[31,0,20],[0,1,0],[48,0,31]]),
}
checks['T_exact']=T==expected['T']
checks['U_exact']=U==expected['U']
checks['TU_commute']=T*U==U*T
checks['V_source_alias_equals_Y_v30']=V==Y_v30==expected['V']
checks['W_exact']=W==expected['W']

for name,A in [('P',P),('Q',Q),('R',R),('S',S),('T',T),('U',U),('V',V),('W',W)]:
    checks[f'det_{name}_one']=A.det()==1

# Primary edge target forms.
targets={name:sp.simplify(A.T*Q0*A) for name,A in [('P',P),('Q',Q),('R',R),('S',S)]}
checks['P_target']=targets['P']==sp.diag(3,-4,-5)
checks['Q_target']=targets['Q']==sp.diag(3,-1,-20)
checks['R_target']=targets['R']==sp.Matrix([[147,0,-180],[0,-1,0],[-180,0,220]])
checks['S_target']=targets['S']==sp.Matrix([[39,-18,0],[-18,8,0],[0,0,-5]])
for name,A in [('T',T),('U',U),('V',V),('W',W)]:
    checks[f'{name}_preserves_Q0']=sp.simplify(A.T*Q0*A)==Q0

# Mod-2 quotient.
def key(M):
    return tuple(int(M[i,j])%2 for i in range(3) for j in range(3))
def mul2(a,b):
    A=sp.Matrix(3,3,list(a)); B=sp.Matrix(3,3,list(b)); C=A*B
    return tuple(int(C[i,j])%2 for i in range(3) for j in range(3))
Ikey=key(sp.eye(3))
gens=[key(P),key(Q),key(R),key(S)]
seen={Ikey}; q=deque([Ikey])
while q:
    a=q.popleft()
    for g in gens:
        b=mul2(a,g)
        if b not in seen:
            seen.add(b); q.append(b)
checks['mod2_image_size_168']=len(seen)==168
checks['R_mod2_equals_Q_mod2']=key(R)==key(Q)
checks['V_mod2_identity']=key(V)==Ikey
checks['W_mod2_identity']=key(W)==Ikey

# Historical Krein boosts for the two kernel periods V and W.
J=sp.diag(1,-1)
F_P=sp.diag(1/sp.sqrt(12),1)
V2=sp.Matrix([[V[0,0],V[0,1]],[V[1,0],V[1,1]]])
gV=sp.simplify(F_P.inv()*V2.inv()*F_P)
checks['gV_exact']=gV==sp.Matrix([[7,-4*sp.sqrt(3)],[-4*sp.sqrt(3),7]])
checks['gV_SO11']=sp.simplify(gV.T*J*gV-J)==sp.zeros(2) and gV.det()==1

F_Q=sp.diag(1/sp.sqrt(12),1/sp.sqrt(5))
W2=sp.Matrix([[W[0,0],W[0,2]],[W[2,0],W[2,2]]])
gW=sp.simplify(F_Q.inv()*W2.inv()*F_Q)
checks['gW_exact']=gW==sp.Matrix([[31,-8*sp.sqrt(15)],[-8*sp.sqrt(15),31]])
checks['gW_SO11']=sp.simplify(gW.T*J*gW-J)==sp.zeros(2) and gW.det()==1
checks['chiV_nonzero']=sp.acosh(7)!=0
checks['chiW_nonzero']=sp.acosh(31)!=0

# Fano/Singer odd memory adapter, inherited from v30.
sigma=[0,4,1,5,2,6,3]
omega=cmath.exp(2j*math.pi/7)
lines={'u1':{0,1,3},'u2':{0,4,5},'u3':{0,2,6}}
def Omega(L):
    Z=sum((1 if sigma[n] in L else -1)*omega**n for n in range(7))
    return Z, (Z**3).imag/(16*math.sqrt(2))
Z1,Om1=Omega(lines['u1']); Z2,Om2=Omega(lines['u2'])
PV=Om1*(4**3-1**3)
PW=Om2*(20**3-5**3)
checks['Omega_u1_norm2_8']=abs(abs(Z1)**2-8)<1e-10
checks['Omega_u2_norm2_8']=abs(abs(Z2)**2-8)<1e-10
checks['V_memory_63Omega_nonzero']=abs(PV)>1e-10 and abs(PV-63*Om1)<1e-10
checks['W_memory_7875Omega_nonzero']=abs(PW)>1e-10 and abs(PW-7875*Om2)<1e-10

# The R/S moves are multi-step historical macros: a single epsilon is not source-canonical.
def conorms(M):
    a,b,c=M[0,0],M[1,1],M[2,2]; k,h,g=M[0,1],M[0,2],M[1,2]
    return [sp.expand(a+k+h),sp.expand(b+k+g),sp.expand(c+h+g),sp.expand(-k),sp.expand(-h),sp.expand(-g)]
con={name:conorms(M) for name,M in [('Q0',Q0),*targets.items()]}
checks['R_not_single_axis_conorm_transition']=sum(1 for x,y in zip(con['Q0'],con['R']) if x!=y)>1
checks['S_not_single_axis_conorm_transition']=sum(1 for x,y in zip(con['Q0'],con['S']) if x!=y)>1

checks={k:bool(v) for k,v in checks.items()}
status='PASS' if all(checks.values()) else 'FAIL'

rows=[]
for name,A,role in [
    ('P',P,'primary boundary species'),('Q',Q,'primary boundary species'),
    ('R',R,'primary doubled boundary species'),('S',S,'primary doubled boundary species'),
    ('T',T,'derived self-automorphism P R P^-1'),('U',U,'derived self-automorphism Q S Q^-1'),
    ('V',V,'primary/derived P-period self-automorphism; legacy alias Y_v30'),
    ('W',W,'primary/derived Q-period self-automorphism')]:
    rows.append({
      'name':name,'role':role,'matrix':str([[int(A[i,j]) for j in range(3)] for i in range(3)]),
      'det':int(A.det()),'mod2':str([[int(A[i,j])%2 for j in range(3)] for i in range(3)]),
      'preserves_Q0':bool(sp.simplify(A.T*Q0*A)==Q0),
      'target_form':str([[int((A.T*Q0*A)[i,j]) for j in range(3)] for i in range(3)])
    })
with (OUT/'historical_lift_generators_v31.csv').open('w',newline='',encoding='utf-8') as fh:
    w=csv.DictWriter(fh,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

cert={
  'phase':'v31-historical-lift-groupoid',
  'status':status,
  'checks':checks,
  'source_normalization':{
    'Q0':[[12,0,0],[0,-1,0],[0,0,-5]],
    'P':[[int(P[i,j]) for j in range(3)] for i in range(3)],
    'Q':[[int(Q[i,j]) for j in range(3)] for i in range(3)],
    'R':[[int(R[i,j]) for j in range(3)] for i in range(3)],
    'S':[[int(S[i,j]) for j in range(3)] for i in range(3)],
    'T':[[int(T[i,j]) for j in range(3)] for i in range(3)],
    'U':[[int(U[i,j]) for j in range(3)] for i in range(3)],
    'V':[[int(V[i,j]) for j in range(3)] for i in range(3)],
    'W':[[int(W[i,j]) for j in range(3)] for i in range(3)],
    'legacy_alias':'V_Selling = Y_v30',
    'relations_proved':['T=P R P^-1','U=Q S Q^-1','T U=U T','V=P DP P^-1 DP','W=Q DQ Q^-1 DQ']
  },
  'finite_vs_infinite':{
    'finite_statement':'Selling supports finite field development modulo repetition for a fixed form.',
    'infinite_statement':'The lifted line/automorphism system is ramified and infinite; its morphisms are alternating integer words in V,W with further R/S-linked systems.',
    'decision':'REFUTED_AS_FINITE_EDGE_CENSUS__PROMOTED_AS_FINITE_GENERATOR_PRESENTATION_WITH_INFINITE_MORPHISMS'
  },
  'mod2':{
    'image_order':len(seen),
    'image':'GL_3(2)',
    'R_equals_Q_in_mod2':key(R)==key(Q),
    'kernel_witnesses':['V','W'],
    'V_mod2':'I','W_mod2':'I'
  },
  'kernel_observables':{
    'V':{
      'legacy_name':'Y_v30','chi_abs_exact':'acosh(7)=log(7+4sqrt(3))',
      'chi_abs_numeric':float(sp.N(sp.acosh(7),16)),
      'memory_exact':'63 Omega_sigma_minus(u1)', 'memory_numeric':PV,
      'memory_status':'PROVEN_CONSTRUCTED_IN_DECLARED_CONORM_LIFT'
    },
    'W':{
      'chi_abs_exact':'acosh(31)=log(31+8sqrt(15))',
      'chi_abs_numeric':float(sp.N(sp.acosh(31),16)),
      'memory_exact':'7875 Omega_sigma_minus(u2)', 'memory_numeric':PW,
      'memory_status':'PROVEN_CONSTRUCTED_IN_DECLARED_CONORM_LIFT'
    }
  },
  'rank_tests':{
    'v29_finite_quotient':'INHERITED: rank Dc=168 < rank[Dc|s]=169; rank D0=rank[D0|chi]=167',
    'complete_lifted_census':'NT_NOT_A_FINITE_MATRIX_PROBLEM_WITHOUT_A_FINITE_FUNDAMENTAL_DOMAIN_OR_PRESENTATION_COMPLEX',
    'generator_level_kernel_test':'PASS_NONZERO_LIFT_OBSERVABLES'
  },
  'typing_guards':[
    'P,Q,R,S are source matrices. T,U,V,W are exact source-supported compositions/self-transformations.',
    'The letter V in the primary scan is retained; v30 legacy label Y is preserved as an alias, not silently rewritten.',
    'R and S are multi-step historical substitutions. No single epsilon_R or epsilon_S is promoted without an explicit primitive path decomposition.',
    'V and W boosts live on differently selected Krein (1,1) planes. Their separate nonzero rapidities prove failure of mod-2 descent, but a single global scalar chi cocycle across mixed fibers remains PENDING-MORPHISM until a common trivialization/transport is supplied.',
    'The W memory formula is the exact analogue of the v30 V/Y two-crossing construction inside the declared Singer/Binet conorm-lift; it is not claimed as a formula written by Selling.',
    'The v29 normalized quotient charge and v31 historical lift charges remain SEPARATED until an explicit cochain pullback/comparison map is constructed.'
  ]
}
(OUT/'historical_lift_groupoid_v31_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
