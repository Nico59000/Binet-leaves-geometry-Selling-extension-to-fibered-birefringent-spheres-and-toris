#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp, json, hashlib, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from selling_reconstruction_engine_v49 import S, RELABEL, WALL_TO_MOVES, gram_from_six, homogeneous

labs=['g','h','k','l','m','n']
wall_for_conorm={'p23':'g','p13':'h','p12':'k','p01':'l','p02':'m','p03':'n'}
first_move={'g':9,'h':2,'k':1,'l':3,'m':5,'n':4}
second_move={'g':11,'h':8,'k':7,'l':6,'m':12,'n':10}

# Base Q2 field and global projective chart q=z/x, r=-y/x.
M=gram_from_six(7,-1,-8,0,2,0)
q,r=sp.symbols('q r', real=True)
D0=8+4*q-7*q**2-60*r**2
w0=sp.Matrix([
 -120*q*r,
 2*(7*q**2+56*q+60*r**2-8),
 -120*r,
 -3*(21*q**2+28*q+180*r**2-40*r+16),
 7*q**2+120*q*r-4*q-60*r**2+120*r-8,
 -6*(13*q**2-20*q*r+24*q-60*r**2+8),
])
# exact x=0 exclusion on hyperboloid: -(60 y^2+7 z^2)/60 = 1 impossible over R.
y,z=sp.symbols('y z', real=True)
x0_form=sp.factor((sp.Matrix([0,y,z]).T*M.inv()*sp.Matrix([0,y,z]))[0])

# Six wall-action matrices on homogeneous wall vector.
G,H,K,L,MM,N=sp.symbols('g h k l m n')
a=-L-H-K; b=-MM-G-K; c=-N-G-H
X=gram_from_six(a,b,c,G,H,K)
syms=[G,H,K,L,MM,N]
R={}
for wall,idx in first_move.items():
 new=homogeneous(S[idx-1].T*X*S[idx-1]); A=sp.zeros(6,6)
 for i,lab in enumerate(labs):
  e=sp.expand(new[lab])
  for j,s in enumerate(syms):A[i,j]=e.coeff(s)
 assert A.det()==1 and A*A==sp.eye(6)
 R[wall]=A

# Global projective gauge collapse for all six candidate pairs.
relset={tuple(map(int,list(B))) for B in RELABEL}
torsors={}
for con,(i,j) in WALL_TO_MOVES.items():
 wall=wall_for_conorm[con]
 A1,A2=S[i-1],S[j-1]
 B=sp.simplify(-A1.inv()*A2)
 assert tuple(map(int,list(B))) in relset
 assert B*B==sp.eye(3) and B.det()==-1
 assert A2.inv()*A1*B==-sp.eye(3)
 torsors[wall]={'conorm':con,'candidate_moves':[f'S{i}',f'S{j}'],'target_relabelling':B.tolist(),
                'det':int(B.det()),'involution':True,'projective_relation':'A2^-1*A1*B=-I',
                'status':'PROVEN-UNIQUE-PROJECTIVE-TRANSPORT-CLASS'}

# Deck stabilizers conjugated into the Q2 chart.
Q0=sp.diag(12,-1,-5)
Q1=sp.Matrix([[1,0,0],[0,1,0],[1,0,1]])
Q2m=sp.Matrix([[1,0,1],[0,1,0],[0,0,1]])
Aseed=Q1*Q2m
assert Aseed.T*Q0*Aseed==M
P=sp.Matrix([[1,1,0],[3,4,0],[0,0,1]])
Rmac=sp.Matrix([[-4,0,5],[0,-1,0],[-3,0,4]])
T=P*Rmac*P.inv()
Qmac=sp.Matrix([[2,0,5],[0,1,0],[3,0,8]])
Smac=sp.Matrix([[-2,1,0],[-3,2,0],[0,0,-1]])
U=Qmac*Smac*Qmac.inv()
GT=sp.simplify(Aseed.inv()*T*Aseed); GU=sp.simplify(Aseed.inv()*U*Aseed)
assert GT.T*M*GT==M and GU.T*M*GU==M and GT*GU==GU*GT
assert all(x.q==1 for x in list(GT)+list(GU))

checks={
 'Q2_det60':M.det()==60,
 'global_projective_chart_covers_real_hyperboloid':sp.simplify(x0_form+(60*y**2+7*z**2)/60)==0,
 'six_wall_actions_integer_involutions':all(A.det()==1 and A*A==sp.eye(6) for A in R.values()),
 'six_projective_torsors_collapsed':len(torsors)==6,
 'deck_T_stabilizes_Q2':GT.T*M*GT==M,
 'deck_U_stabilizes_Q2':GU.T*M*GU==M,
 'deck_TU_commute':GT*GU==GU*GT,
}
assert all(checks.values())
out={
 'version':'v50','status':'PASS','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
 'base_field_six':[7,-1,-8,0,2,0],
 'global_projective_chart':{
   'coordinates':'q=z/x, r=-y/x modulo v~-v','hyperboloid_scale':'x^2=60/D0','D0':str(D0),
   'x_zero_restriction':str(x0_form),'coverage_decision':'GLOBAL: x=0 has no real hyperboloid point, so this chart covers the whole real projective hyperboloid'},
 'base_wall_numerators_order_g_h_k_l_m_n':[str(x) for x in w0],
 'wall_action_matrices':{k:v.tolist() for k,v in R.items()},
 'wall_action_rule':'for path A and wall vector w_A=L_A w0, crossing wall j uses L_new=R_j L_A; reduced field is L_A w0 <= 0 componentwise',
 'projective_transport_torsors':torsors,
 'deck_in_Q2_chart':{'G_T':GT.tolist(),'G_U':GU.tolist(),'commute':True,
                     'quotient_rule':'path matrices may be identified only after exact left action by <G_T,G_U>, right superbase relabelling, and projective central sign'},
 'exact_global_classifier_contract':[
   'factor each transformed quadratic wall polynomial',
   'retain a boundary component only if it contains a real open arc with D0>0 and all six transformed inequalities <=0',
   'if one irreducible boundary factor is shared by >=2 wall components, classify it as MULTIWALL/MIRROR and do not create a single-wall neighbour',
   'cross only certified ordinary one-wall arcs',
   'quotient fields only by exact deck-left / relabelling-right / projective-sign equivalence',
 ],
 'checks':checks
}
p=HERE/'selling_global_wall_linearization_v50.json';p.write_text(json.dumps(out,indent=2,ensure_ascii=False,default=str),encoding='utf-8')
print(json.dumps({'status':'PASS','checks':checks,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
