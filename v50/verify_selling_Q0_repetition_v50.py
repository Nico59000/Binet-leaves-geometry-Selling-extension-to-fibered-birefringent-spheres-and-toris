#!/usr/bin/env python3
from __future__ import annotations
import json, csv
from pathlib import Path
import sympy as sp

OUT=Path(__file__).resolve().parent
x,y,z=sp.symbols('xi eta zeta', real=True)

Q0=sp.diag(12,-1,-5)

# Selling 1877, Fig.1, pp. 189--191: route substitutions and their printed factor order.
P1=sp.Matrix([[1,0,0],[3,1,0],[0,0,1]])
P2=sp.Matrix([[1,1,0],[0,1,0],[0,0,1]])
P=P1*P2

Q1=sp.Matrix([[1,0,0],[0,1,0],[1,0,1]])
Q2=sp.Matrix([[1,0,1],[0,1,0],[0,0,1]])
Q3=sp.Matrix([[1,0,0],[0,1,0],[1,0,1]])
Q4=sp.Matrix([[1,0,2],[0,1,0],[0,0,1]])
Q=Q1*Q2*Q3*Q4

R1=sp.Matrix([[1,0,1],[0,1,0],[0,0,1]])
R2=sp.Matrix([[1,0,0],[0,1,0],[1,0,1]])
R3=sp.Matrix([[-1,0,0],[0,-1,0],[-1,0,1]])
R4=sp.Matrix([[1,0,0],[0,1,0],[-1,0,1]])
R5=sp.Matrix([[1,0,-1],[0,1,0],[0,0,1]])
R=R1*R2*R3*R4*R5

S1=sp.Matrix([[1,0,0],[1,1,0],[0,0,1]])
S2=sp.Matrix([[-1,1,0],[0,1,0],[0,0,-1]])
S3=sp.Matrix([[1,0,0],[-1,1,0],[0,0,1]])
SMAC=S1*S2*S3

P_EXPECT=sp.Matrix([[1,1,0],[3,4,0],[0,0,1]])
Q_EXPECT=sp.Matrix([[2,0,5],[0,1,0],[3,0,8]])
R_EXPECT=sp.Matrix([[-4,0,5],[0,-1,0],[-3,0,4]])
S_EXPECT=sp.Matrix([[-2,1,0],[-3,2,0],[0,0,-1]])

# Superbase wall combinatorics.
WALL_EDGE={'g':(2,3),'h':(1,3),'k':(1,2),'l':(0,1),'m':(0,2),'n':(0,3)}
WALLS=tuple(WALL_EDGE)
def opposite(u,v): return set(WALL_EDGE[u]).isdisjoint(WALL_EDGE[v])
OPPOSITE={u:next(v for v in WALLS if v!=u and opposite(u,v)) for u in WALLS}

def wall_pair_type(u,v):
    if u==v: return 'SAME'
    return 'FISSURE' if opposite(u,v) else 'CROSSING'

def wall_polynomials(M):
    a,b,c=M[0,0],M[1,1],M[2,2]
    k,h,g=M[0,1],M[0,2],M[1,2]
    s=x+y+z
    return {
      'g':sp.expand(2*y*z-g),
      'h':sp.expand(2*x*z-h),
      'k':sp.expand(2*x*y-k),
      'l':sp.expand(a+h+k-2*x*s),
      'm':sp.expand(b+k+g-2*y*s),
      'n':sp.expand(c+h+g-2*z*s),
    }

def hyperboloid_poly(M):
    v=sp.Matrix([x,y,z])
    return sp.factor((v.T*M.inv()*v)[0]-1)

def fissure_witness_check(M,u,v,point):
    assert wall_pair_type(u,v)=='FISSURE'
    subs={x:point[0],y:point[1],z:point[2]}
    ww=wall_polynomials(M)
    hp=sp.simplify(hyperboloid_poly(M).subs(subs))
    vals={w:sp.simplify(p.subs(subs)) for w,p in ww.items()}
    return hp==0 and vals[u]==0 and vals[v]==0 and all(bool(sp.Le(q,0)) for q in vals.values()), vals

# Canonical numerical form under the 24 permutations of a superbase.
def relabelling_matrices():
    from itertools import permutations
    V=[sp.Matrix([-1,-1,-1]),sp.eye(3)[:,0],sp.eye(3)[:,1],sp.eye(3)[:,2]]
    mats=[]; seen=set()
    for pi in permutations(range(4)):
        A=sp.Matrix.hstack(V[pi[1]],V[pi[2]],V[pi[3]])
        if abs(int(A.det()))!=1: continue
        key=tuple(map(int,list(A)))
        if key not in seen: seen.add(key);mats.append(A)
    return mats
RELABEL=relabelling_matrices()
def numerical_key(M):
    keys=[]
    for A in RELABEL:
        N=A.T*M*A
        keys.append(tuple(map(int,list(N))))
    return min(keys)

def six(M):
    return [int(M[0,0]),int(M[1,1]),int(M[2,2]),int(M[1,2]),int(M[0,2]),int(M[0,1])]

def trace_route(name,factors):
    A=sp.eye(3); first={numerical_key(Q0):0}; rows=[]
    rows.append({'route':name,'step':0,'factor':'START','a':12,'b':-1,'c':-5,'g':0,'h':0,'k':0,
                 'exact_Q0':True,'first_numerical_repeat':None,'central_exact_symmetry':False})
    prev=Q0
    for i,(label,F) in enumerate(factors,1):
        A=A*F; M=sp.simplify(A.T*Q0*A); nk=numerical_key(M)
        rep=first.get(nk)
        central=(M==prev)
        co=six(M)
        rows.append({'route':name,'step':i,'factor':label,'a':co[0],'b':co[1],'c':co[2],
                     'g':co[3],'h':co[4],'k':co[5],'exact_Q0':bool(M==Q0),
                     'first_numerical_repeat':rep,'central_exact_symmetry':central})
        first.setdefault(nk,i); prev=M
    return A,rows

T_FACT=[('P1',P1),('P2',P2),('R1',R1),('R2',R2),('R3',R3),('R4',R4),('R5',R5),
        ('P2^-1',P2.inv()),('P1^-1',P1.inv())]
U_FACT=[('Q1',Q1),('Q2',Q2),('Q3',Q3),('Q4',Q4),('S1',S1),('S2',S2),('S3',S3),
        ('Q4^-1',Q4.inv()),('Q3^-1',Q3.inv()),('Q2^-1',Q2.inv()),('Q1^-1',Q1.inv())]
T,trows=trace_route('T=P R P^-1',T_FACT)
U,urows=trace_route('U=Q S Q^-1',U_FACT)

# Source-side local symmetry continuations from p.191.
JR=sp.diag(-1,1,-1)
JQ=sp.diag(-1,-1,1)
V0=P*JR*P.inv()
V=V0*JR
W=Q*JQ*Q.inv()*JQ

# An exact fissure witness on the first nontrivial P-route marking.
M_P1=P1.T*Q0*P1
FW=sp.sqrt(3)
fiss_ok_pos,fiss_vals_pos=fissure_witness_check(M_P1,'g','l',(FW,-FW,0))
fiss_ok_neg,fiss_vals_neg=fissure_witness_check(M_P1,'g','l',(-FW,FW,0))

checks={
 'wall_pair_partition_3_fissure_12_crossing': sum(wall_pair_type(a,b)=='FISSURE' for i,a in enumerate(WALLS) for b in WALLS[i+1:])==3,
 'opposite_pairs_exact': {tuple(sorted((u,v))) for u,v in OPPOSITE.items()}=={('g','l'),('h','m'),('k','n')},
 'P_factorization_exact':P==P_EXPECT,
 'Q_factorization_exact':Q==Q_EXPECT,
 'R_factorization_exact':R==R_EXPECT,
 'S_factorization_exact':SMAC==S_EXPECT,
 'P_target_source_form':P.T*Q0*P==sp.diag(3,-4,-5),
 'Q_target_source_form':Q.T*Q0*Q==sp.diag(3,-1,-20),
 'R3_involution':R3*R3==sp.eye(3),
 'S2_involution':S2*S2==sp.eye(3),
 'T_closes_Q0':T.T*Q0*T==Q0,
 'U_closes_Q0':U.T*Q0*U==Q0,
 'TU_commute':T*U==U*T,
 'T_trace_returns_exact_Q0':trows[-1]['exact_Q0'] and trows[-1]['first_numerical_repeat']==0,
 'U_trace_returns_exact_Q0':urows[-1]['exact_Q0'] and urows[-1]['first_numerical_repeat']==0,
 'T_central_reversion_exact':trows[5]['central_exact_symmetry'] and [r['first_numerical_repeat'] for r in trows[5:]]==[4,3,2,1,0],
 'U_central_reversion_exact':urows[6]['central_exact_symmetry'] and [r['first_numerical_repeat'] for r in urows[6:]]==[5,4,3,2,1,0],
 'JR_stabilizes_P_target':JR.T*(P.T*Q0*P)*JR==P.T*Q0*P,
 'JQ_stabilizes_Q_target':JQ.T*(Q.T*Q0*Q)*JQ==Q.T*Q0*Q,
 'V0_matches_source_reciprocal':V0==sp.Matrix([[-7,2,0],[-24,7,0],[0,0,-1]]),
 'V_matches_source_total':V==sp.Matrix([[7,2,0],[24,7,0],[0,0,1]]),
 'V_stabilizes_Q0':V.T*Q0*V==Q0,
 'W_stabilizes_Q0':W.T*Q0*W==Q0,
 'fissure_gl_witness_positive':fiss_ok_pos,
 'fissure_gl_witness_negative':fiss_ok_neg,
}
status='PASS' if all(checks.values()) else 'FAIL'

branch={
 'version':'v50','status':status,
 'source_anchor':'Selling 1877, IV.b--c, esp. pp.156--170 (PDF pp.6--18)',
 'wall_edge_model':{k:list(v) for k,v in WALL_EDGE.items()},
 'opposite_wall':OPPOSITE,
 'pair_partition':{
   'fissure_pairs':[['g','l'],['h','m'],['k','n']],
   'crossing_pairs':[[a,b] for i,a in enumerate(WALLS) for b in WALLS[i+1:] if wall_pair_type(a,b)=='CROSSING']},
 'local_models':{
   'fissure':{
      'criterion':'two opposite superbase-edge wall coefficients vanish at a reduced hyperboloid point',
      'standard_example':'h=m=0',
      'three_local_lines':['h=0','h-m=0','m=0'],
      'field_boundary_continuation':'if no reduced crossing occurs on the current branch, continue from wall x through the fissure on opposite wall opp(x)',
   },
   'crossing':{
      'criterion':'two adjacent superbase-edge wall coefficients vanish at a reduced hyperboloid point',
      'standard_example':'h=k=0',
      'six_sector_lines':['h=0','h+k=0','k=0','retrograde(h=0)','retrograde(h+k=0)','retrograde(k=0)'],
   }
 },
 'branch_order_rule':{
   'from_crossing_xy_following_x':'test the three other walls adjacent to x (excluding y) for reduced crossings on the same branch; the consecutive one is followed first',
   'fallback':'if none exists on that branch, x terminates at fissure (x,opp(x)) and the boundary continues on opp(x)',
   'two_crossing_simplification':'when exactly two crossing points lie on wall x, Selling permits treating them as if they belonged to the same branch for field grouping',
   'closure_guard':'the final boundary line of a field must reconnect to its first; otherwise the field is not serialized as historical',
 },
 'standard_h_branch_discriminator':{
   'case_B_pos_a_c_neg':'two branches; distinguish by the signs of xi and zeta',
   'case_B_a_c_neg':'two parts; distinguish by the sign of adjoint H',
   'case_a_c_opposite_signs':'one branch',
   'case_a_c_positive':'wall impossible',
   'other_cases':'NT/FAIL-CLOSED until transformed to a certified standard marking',
 },
 'exact_fissure_test':{
   'hyperboloid':'v^T M^{-1} v = 1',
   'walls':'H_w(v)<=0 for w in {g,h,k,l,m,n}',
   'event':'H_x(v)=H_opp(x)(v)=0',
   'accept':'all six reduction inequalities hold exactly',
 },
 'certified_witness':{
   'form_after_P1':[[int(M_P1[i,j]) for j in range(3)] for i in range(3)],
   'pair':['g','l'],
   'points':['(sqrt(3),-sqrt(3),0)','(-sqrt(3),sqrt(3),0)'],
   'wall_values_at_either':{k:str(v) for k,v in fiss_vals_pos.items()},
 },
 'checks':checks,
}

closure={
 'version':'v50','status':status,
 'source_anchor':'Selling 1877, Fig.1 repetition construction, pp.189--191 (PDF pp.37--39)',
 'Q0':[[12,0,0],[0,-1,0],[0,0,-5]],
 'macros':{
   'P':{'matrix':P.tolist(),'factor_order':['P1','P2'],'target_form':(P.T*Q0*P).tolist()},
   'Q':{'matrix':Q.tolist(),'factor_order':['Q1','Q2','Q3','Q4'],'target_form':(Q.T*Q0*Q).tolist()},
   'R':{'matrix':R.tolist(),'factor_order':['R1','R2','R3','R4','R5'],'central_symmetry':'R3'},
   'S':{'matrix':SMAC.tolist(),'factor_order':['S1','S2','S3'],'central_symmetry':'S2'},
 },
 'deck_repetition':{
   'T':'P R P^-1','T_matrix':T.tolist(),'T_stabilizes_Q0':True,
   'U':'Q S Q^-1','U_matrix':U.tolist(),'U_stabilizes_Q0':True,
   'TU_equals_UT':True,
 },
 'additional_source_stabilizers':{
   'V0_P_reciprocal':V0.tolist(),
   'V_after_final_Q_symmetry':V.tolist(),
   'W_Q_analogue':W.tolist(),
 },
 'route_reversion':{
   'T_numerical_repeat_signature':[r['first_numerical_repeat'] for r in trows],
   'U_numerical_repeat_signature':[r['first_numerical_repeat'] for r in urows],
   'T_form_sequence':[six((sp.eye(3) if i==0 else sp.eye(3))) for i in []],
   'interpretation':'both source routes reach an exact central marking symmetry and then retrace the previously seen numerical forms in reverse order, ending at Q0 exactly',
 },
 'closure_decision':'PROVEN-SOURCE-REPETITION-CLOSURE at the P/Q/R/S route level',
 'microcell_decision':'NOT-YET-MATERIALIZED: macro repetition closure is now certified, but C0/C1/C2 require expansion of every internal field boundary event and transport',
 'historical_H1_H2':'GUARDED/NOT-RUN',
 'checks':checks,
}

(OUT/'selling_fissure_branch_automaton_v50.json').write_text(json.dumps(branch,indent=2,ensure_ascii=False,default=str),encoding='utf-8')
(OUT/'selling_Q0_repetition_closure_v50.json').write_text(json.dumps(closure,indent=2,ensure_ascii=False,default=str),encoding='utf-8')

rows=trows+urows
with (OUT/'selling_Q0_route_trace_v50.csv').open('w',newline='',encoding='utf-8') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows[0].keys()));wr.writeheader();wr.writerows(rows)

report={'version':'v50','status':status,'check_count':len(checks),'pass_count':sum(checks.values()),'checks':checks,
        'outputs':['selling_fissure_branch_automaton_v50.json','selling_Q0_repetition_closure_v50.json','selling_Q0_route_trace_v50.csv']}
(OUT/'selling_Q0_repetition_verification_v50.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(report,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
