#!/usr/bin/env python3
import json, math, cmath
from pathlib import Path
import sympy as sp

checks={}

Q0=sp.diag(12,-1,-5)
P=sp.Matrix([[1,1,0],[3,4,0],[0,0,1]])
Qsub=sp.Matrix([[2,0,5],[0,1,0],[3,0,8]])
Y=sp.Matrix([[7,2,0],[24,7,0],[0,0,1]])
D=sp.diag(1,-1,1)
K=sp.eye(3); K[0,1]=2

QP=P.T*Q0*P
QQ=Qsub.T*Q0*Qsub
QK=K.T*Q0*K

checks['P_target_exact']=QP==sp.diag(3,-4,-5)
checks['Q_target_exact']=QQ==sp.diag(3,-1,-20)
checks['Y_preserves_Q0']=Y.T*Q0*Y==Q0
checks['Y_factorization']=Y==P*D*P.inv()*D
checks['det_P']=P.det()==1
checks['det_Q']=Qsub.det()==1
checks['det_Y']=Y.det()==1

def conorms(Q):
    a,b,c=Q[0,0],Q[1,1],Q[2,2]
    k,h,g=Q[0,1],Q[0,2],Q[1,2]
    return {
        "p01":sp.expand(a+k+h),
        "p02":sp.expand(b+k+g),
        "p03":sp.expand(c+h+g),
        "p12":sp.expand(-k),
        "p13":sp.expand(-h),
        "p23":sp.expand(-g),
    }
c0=conorms(Q0); cP=conorms(QP); cQ=conorms(QQ); cK=conorms(QK)
checks['conorms_Q0']=c0=={'p01':12,'p02':-1,'p03':-5,'p12':0,'p13':0,'p23':0}
checks['conorms_QP']=cP=={'p01':3,'p02':-4,'p03':-5,'p12':0,'p13':0,'p23':0}
checks['conorms_QQ']=cQ=={'p01':3,'p02':-1,'p03':-20,'p12':0,'p13':0,'p23':0}

# P and Q historical conorm amplitudes
epsP0,epsP1,epsQ0,epsQ1=1,4,5,20
checks['historical_eps_P_1_to_4']=(-c0['p02']==epsP0 and -cP['p02']==epsP1)
checks['historical_eps_Q_5_to_20']=(-c0['p03']==epsQ0 and -cQ['p03']==epsQ1)

# Krein transitions
J=sp.diag(1,-1)
F0P=sp.diag(1/sp.sqrt(12),1)
F1P=sp.diag(1/sp.sqrt(3),sp.Rational(1,2))
gP=sp.simplify(F1P.inv()*P[:2,:2].inv()*F0P)
checks['gP_exact']=gP==sp.Matrix([[2,-sp.sqrt(3)],[-sp.sqrt(3),2]])
checks['gP_O11']=sp.simplify(gP.T*J*gP-J)==sp.zeros(2)

F0Q=sp.diag(1/sp.sqrt(12),1/sp.sqrt(5))
F1Q=sp.diag(1/sp.sqrt(3),1/sp.sqrt(20))
Q2=sp.Matrix([[2,5],[3,8]])
gQ=sp.simplify(F1Q.inv()*Q2.inv()*F0Q)
checks['gQ_exact']=gQ==sp.Matrix([[4,-sp.sqrt(15)],[-sp.sqrt(15),4]])
checks['gQ_O11']=sp.simplify(gQ.T*J*gQ-J)==sp.zeros(2)

gY=sp.simplify(F0P.inv()*Y[:2,:2].inv()*F0P)
checks['gY_exact']=gY==sp.Matrix([[7,-4*sp.sqrt(3)],[-4*sp.sqrt(3),7]])
checks['gY_SOplus11']=sp.simplify(gY.T*J*gY-J)==sp.zeros(2) and gY.det()==1 and gY[0,0]>0
checks['chiY_double_chiP']=abs(float(sp.N(sp.acosh(7)-2*sp.acosh(2),30)))<1e-20
checks['sinh_chiY']=sp.simplify(sp.sinh(sp.acosh(7))-4*sp.sqrt(3))==0

# Historical Singer observable and P-period memory.
sigma=[0,4,1,5,2,6,3]
L={0,1,3} # u1
omega=cmath.exp(2j*math.pi/7)
Z=sum((1 if sigma[n] in L else -1)*omega**n for n in range(7))
Omega=(Z**3).imag/(16*math.sqrt(2))
period=Omega*(4**3-1**3)
checks['Omega_norm_u1']=abs(abs(Z)**2-8)<1e-10
checks['P_period_memory_nonzero']=abs(period)>1e-8
checks['P_period_formula_63Omega']=abs(period-63*Omega)<1e-12

# Descent obstruction to v29 mod2 quotient.
checks['K_mod2_identity']=all(int(K[i,j])%2==(1 if i==j else 0) for i in range(3) for j in range(3))
checks['K_changes_form']=QK!=Q0
checks['K_changes_conorm_amplitude_data']=cK!=c0
checks['Y_mod2_identity']=all(int(Y[i,j])%2==(1 if i==j else 0) for i in range(3) for j in range(3))
checks['Y_has_nonzero_historical_boost']=float(sp.N(sp.acosh(7),16))>0

# Strong consequence: neither epsilon nor chi is a function of the GL3(2) vertex alone.
checks['epsilon_not_descend_GL32']=checks['K_mod2_identity'] and checks['K_changes_conorm_amplitude_data']
checks['chi_not_descend_GL32']=checks['Y_mod2_identity'] and checks['Y_has_nonzero_historical_boost']

checks={k: bool(v) for k,v in checks.items()}
status='PASS' if all(checks.values()) else 'FAIL'
cert={
  "phase":"v30-explicit-historical-selling-form",
  "status":status,
  "checks":checks,
  "historical_form":{
    "Q0":[[12,0,0],[0,-1,0],[0,0,-5]],
    "P":[[1,1,0],[3,4,0],[0,0,1]],
    "Q":[[2,0,5],[0,1,0],[3,0,8]],
    "Y":[[7,2,0],[24,7,0],[0,0,1]]
  },
  "actual_amplitudes":{
    "P_side":[1,4],
    "Q_side":[5,20],
    "P_period_memory":"63 Omega_sigma_minus(u1)",
    "P_period_memory_numeric":period
  },
  "actual_krein":{
    "gP":"[[2,-sqrt(3)],[-sqrt(3),2]]",
    "chiP":"-acosh(2)",
    "gQ":"[[4,-sqrt(15)],[-sqrt(15),4]]",
    "chiQ":"-acosh(4)",
    "gY":"[[7,-4sqrt(3)],[-4sqrt(3),7]]",
    "chiY":"-acosh(7)",
    "abs_chiY_numeric":float(sp.N(sp.acosh(7),16))
  },
  "decisions":{
    "v29_unit_epsilon_can_be_replaced_on_same_GL32_quotient":"REFUTED",
    "v29_flat_krein_can_be_replaced_on_same_GL32_quotient":"REFUTED",
    "reason":"both data live in congruence-kernel/lift fibers",
    "historical_Y_boost":"PROVEN_NONZERO",
    "historical_P_period_memory_in_declared_conorm_lift":"PROVEN_NONZERO",
    "full_fig1_chamber_census":"NT_NOT_EXHAUSTED"
  },
  "guards":[
    "Selling 1877 supplies Q0, P, Q and Y; the Binet/Singer tau and odd cochain are later project adapters.",
    "The period 63 Omega(u1) is proved in the declared two-crossing conorm-lift model Y=P D P^{-1} D, not claimed as a formula written by Selling.",
    "The explicit kernel matrices show that mod-2 congruence forgets exactly the data needed for physical amplitudes and continuous Krein rapidity.",
    "No claim is made that the three-edge skeleton exhausts figure 1."
  ]
}
Path('/mnt/data/verify_explicit_historical_form_v30_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False))
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
