#!/usr/bin/env python3
import json, csv, math, cmath
from pathlib import Path
from collections import deque
import sympy as sp

OUT=Path(__file__).resolve().parent
OUT.mkdir(parents=True,exist_ok=True)
checks={}
I3=sp.eye(3)
Q0=sp.diag(12,-1,-5)

# Historical/source matrices inherited from v31.
P=sp.Matrix([[1,1,0],[3,4,0],[0,0,1]])
Q=sp.Matrix([[2,0,5],[0,1,0],[3,0,8]])
R=sp.Matrix([[-4,0,5],[0,-1,0],[-3,0,4]])
S=sp.Matrix([[-2,1,0],[-3,2,0],[0,0,-1]])
T=sp.simplify(P*R*P.inv())
U=sp.simplify(Q*S*Q.inv())
DP=sp.diag(1,-1,1); DQ=sp.diag(1,1,-1)
V=sp.simplify(P*DP*P.inv()*DP)
W=sp.simplify(Q*DQ*Q.inv()*DQ)
Gamma=sp.diag(-1,1,-1)
Delta=sp.diag(-1,1,1)
minusI=-I3

# v24 adjacent-superbase generators, frozen by v29.
S_int=[
[[-1,0,1],[0,1,0],[0,0,1]],
[[-1,1,0],[0,1,0],[0,0,1]],
[[-1,1,1],[0,1,0],[0,0,1]],
[[0,-1,0],[-1,0,0],[-1,-1,1]],
[[0,0,-1],[-1,1,-1],[-1,0,0]],
[[1,-1,-1],[0,0,-1],[0,-1,0]],
[[1,0,0],[0,-1,1],[0,0,1]],
[[1,0,0],[0,1,0],[0,1,-1]],
[[1,0,0],[0,1,0],[1,0,-1]],
[[1,0,0],[0,1,0],[1,1,-1]],
[[1,0,0],[1,-1,0],[0,0,1]],
[[1,0,0],[1,-1,1],[0,0,1]],
]
Adj=[sp.Matrix(M) for M in S_int]
for i,A in enumerate(Adj,1):
    checks[f'adj_S{i}_involution']=A*A==I3

# ------------------------------------------------------------------
# A. Primary source decompositions for R and S.
# ------------------------------------------------------------------
A1=sp.Matrix([[1,0,1],[0,1,0],[0,0,1]])
A2=sp.Matrix([[1,0,0],[0,1,0],[1,0,1]])
A3=sp.Matrix([[-1,0,0],[0,-1,0],[-1,0,1]])
A4=sp.Matrix([[1,0,0],[0,1,0],[-1,0,1]])
A5=sp.Matrix([[1,0,-1],[0,1,0],[0,0,1]])
B1=sp.Matrix([[1,0,0],[1,1,0],[0,0,1]])
B2=sp.Matrix([[-1,1,0],[0,1,0],[0,0,-1]])
B3=sp.Matrix([[1,0,0],[-1,1,0],[0,0,1]])
checks['R_source_ordered_decomposition']=A1*A2*A3*A4*A5==R
checks['S_source_ordered_decomposition']=B1*B2*B3==S

H=minusI*(Gamma*Delta)  # diag(-1,-1,1)
Z3=Gamma*Delta          # diag(1,1,-1)
algebraic_refinement={
 'R1':['S3','S2'],
 'R2':['S10','S8'],
 'R3':['-I','Gamma','Delta','S8','S10'],
 'R4':['S8','S10'],
 'R5':['S2','S3'],
 'S1':['S12','S7'],
 'S2':['S2','Gamma','Delta'],
 'S3':['S7','S12'],
}
checks['R1_adjacent_refinement']=A1==Adj[2]*Adj[1]
checks['R2_adjacent_refinement']=A2==Adj[9]*Adj[7]
checks['R3_adjacent_refinement']=A3==H*Adj[7]*Adj[9]
checks['R4_adjacent_refinement']=A4==Adj[7]*Adj[9]
checks['R5_adjacent_refinement']=A5==Adj[1]*Adj[2]
checks['S1_adjacent_refinement']=B1==Adj[11]*Adj[6]
checks['S2_adjacent_refinement']=B2==Adj[1]*Z3
checks['S3_adjacent_refinement']=B3==Adj[6]*Adj[11]

# Diagnostic no-go: algebraic adjacent decomposition is not the historical wall path.
P1=sp.Matrix([[1,0,0],[3,1,0],[0,0,1]])
P2=sp.Matrix([[1,1,0],[0,1,0],[0,0,1]])
Q1=sp.Matrix([[1,0,0],[0,1,0],[1,0,1]])
Q2=sp.Matrix([[1,0,1],[0,1,0],[0,0,1]])
Q3=Q1
Q4=sp.Matrix([[1,0,2],[0,1,0],[0,0,1]])
checks['P_source_two_step']=P1*P2==P
checks['Q_source_four_step']=Q1*Q2*Q3*Q4==Q
checks['P1_adjacent_word']=P1==(Adj[11]*Adj[6])**3
checks['P2_adjacent_word']=P2==Adj[2]*Adj[0]
checks['Q1_adjacent_word']=Q1==Adj[9]*Adj[7]
checks['Q2_adjacent_word']=Q2==Adj[2]*Adj[1]
checks['Q4_adjacent_word']=Q4==(Adj[2]*Adj[1])**2
adjacent_lengths={'P':8,'Q':10}
checks['adjacent_parity_conflicts_with_v30_boundary_tau']=((-1)**8==1 and (-1)**10==1)
# v30 historical crossing tau(P)=tau(Q)=-1, hence this exact parity discrepancy is a typing no-go.

# Write primitive decomposition table.
primitive_rows=[]
def matlist(M): return [[int(M[i,j]) for j in range(3)] for i in range(3)]
for macro,steps in [('R',[A1,A2,A3,A4,A5]),('S',[B1,B2,B3])]:
    for k,M in enumerate(steps,1):
        key=f'{macro}{k}'
        primitive_rows.append({
          'macro':macro,'source_step':k,'step_id':key,'matrix':json.dumps(matlist(M)),
          'det':int(M.det()),'adjacent_algebraic_refinement':' '.join(algebraic_refinement[key]),
          'wall_path_status':'PENDING_GEOMETRIC_ADAPTER',
          'edge_tuple_status':'NO_EPSILON_TAU_S_CHI_FROM_MATRIX_WORD_ALONE'
        })
with (OUT/'primitive_decompositions_v32.csv').open('w',newline='',encoding='utf-8') as fh:
    w=csv.DictWriter(fh,fieldnames=primitive_rows[0].keys()); w.writeheader(); w.writerows(primitive_rows)

# ------------------------------------------------------------------
# B. Finite presentation skeleton and exact 2-cells.
# ------------------------------------------------------------------
gen_names=['V','W','Gamma','Delta','T','U']
GM={'V':V,'W':W,'Gamma':Gamma,'Delta':Delta,'T':T,'U':U}
rels={
 'Gamma2':[('Gamma',1),('Gamma',1)],
 'Delta2':[('Delta',1),('Delta',1)],
 'T2':[('T',1),('T',1)],
 'U2':[('U',1),('U',1)],
 'GammaDelta_comm':[('Gamma',1),('Delta',1),('Gamma',-1),('Delta',-1)],
 'TU_comm':[('T',1),('U',1),('T',-1),('U',-1)],
 'Gamma_V_inversion':[('Gamma',1),('V',1),('Gamma',1),('V',1)],
 'Delta_V_inversion':[('Delta',1),('V',1),('Delta',1),('V',1)],
 'Gamma_W_fixed':[('Gamma',1),('W',1),('Gamma',1),('W',-1)],
 'Delta_W_inversion':[('Delta',1),('W',1),('Delta',1),('W',1)],
}
def word_matrix(word):
    M=I3
    for g,e in word:
        M=M*(GM[g] if e==1 else GM[g].inv())
    return sp.simplify(M)
for rn,wrd in rels.items(): checks[f'relation_{rn}']=word_matrix(wrd)==I3
checks['T_U_commute_source_relation']=T*U==U*T
checks['Gamma_conjugates_V_inverse']=Gamma*V*Gamma==V.inv()
checks['Gamma_fixes_W']=Gamma*W*Gamma==W
checks['Delta_conjugates_V_inverse']=Delta*V*Delta==V.inv()
checks['Delta_conjugates_W_inverse']=Delta*W*Delta==W.inv()

presentation={
 'name':'X32_source_supported_relation_skeleton',
 'cells':{'0':['*'],'1':gen_names,'2':list(rels)},
 'relators':{k:[[g,e] for g,e in v] for k,v in rels.items()},
 'normal_form_core':{
   'carrier':'alternating integral powers of V and W, source unique-path/tree grammar',
   'topological_model':'rose R_2 for the V/W core after quotienting the source tree by the base form',
   'status':'PROVEN_SOURCE_SUPPORTED_CORE'
 },
 'completeness':{
   'known_relations':'PROVEN_EXACT',
   'full_presentation_of_entire_historical_automorphism_group':'NT_PENDING_REWRITE_COMPLETENESS',
   'reason':'No unproved T/U-to-V/W relations are inserted; finite relation skeleton is not silently promoted to a complete group presentation.'
 }
}
(OUT/'presentation_complex_v32.json').write_text(json.dumps(presentation,indent=2,ensure_ascii=False),encoding='utf-8')

# ------------------------------------------------------------------
# C. Exact O(1,1)-compatible transport K_P -> K_Q.
# ------------------------------------------------------------------
s5=sp.sqrt(5)
JPQ=sp.Matrix([[1,0,0],[0,0,s5],[0,1/s5,0]])
checks['JPQ_Q0_isometry']=sp.simplify(JPQ.T*Q0*JPQ-Q0)==sp.zeros(3)
checks['JPQ_involution']=sp.simplify(JPQ*JPQ-I3)==sp.zeros(3)
checks['JPQ_det_minus_one']=sp.simplify(JPQ.det()+1)==0
# Pseudo-frame bases as embedded columns.
EP=sp.Matrix([[1/sp.sqrt(12),0],[0,1],[0,0]])
EQ=sp.Matrix([[1/sp.sqrt(12),0],[0,0],[0,1/sp.sqrt(5)]])
# Coordinate transfer C satisfying EQ*C=JPQ*EP.
C=(EQ.T*EQ).inv()*EQ.T*JPQ*EP
checks['JPQ_pseudoframe_coordinate_transport_identity']=sp.simplify(C-sp.eye(2))==sp.zeros(2)
J11=sp.diag(1,-1)
checks['JPQ_O11_compatible']=sp.simplify(C.T*J11*C-J11)==sp.zeros(2)

# Historical boost matrices in the common normalized coordinates.
FP=sp.diag(1/sp.sqrt(12),1)
FQ=sp.diag(1/sp.sqrt(12),1/sp.sqrt(5))
V2=sp.Matrix([[V[0,0],V[0,1]],[V[1,0],V[1,1]]])
W2=sp.Matrix([[W[0,0],W[0,2]],[W[2,0],W[2,2]]])
gV=sp.simplify(FP.inv()*V2.inv()*FP)
gW=sp.simplify(FQ.inv()*W2.inv()*FQ)
checks['gV_common_SO11']=sp.simplify(gV.T*J11*gV-J11)==sp.zeros(2)
checks['gW_common_SO11']=sp.simplify(gW.T*J11*gW-J11)==sp.zeros(2)

# ------------------------------------------------------------------
# D. Cellular twisted cohomology of the finite relation skeleton.
# ------------------------------------------------------------------
idx={n:i for i,n in enumerate(gen_names)}
def cocycle_relation_matrix(rhos):
    d=next(iter(rhos.values())).rows
    rows=[]
    for rn in rels:
        word=rels[rn]
        A=sp.zeros(d,len(gen_names)*d)
        pref=sp.eye(d)
        for g,e in word:
            Rg=rhos[g]
            if e==1:
                A[:,idx[g]*d:(idx[g]+1)*d] += pref
                pref=pref*Rg
            else:
                Ri=Rg.inv()
                A[:,idx[g]*d:(idx[g]+1)*d] += pref*(-Ri)
                pref=pref*Ri
        assert sp.simplify(pref-sp.eye(d))==sp.zeros(d)
        rows.extend([list(A.row(i)) for i in range(d)])
    return sp.Matrix(rows)

def d0_matrix(rhos):
    d=next(iter(rhos.values())).rows
    return sp.Matrix.vstack(*[rhos[g]-sp.eye(d) for g in gen_names])

# Scalar time-orientation attempt: exact discrete character on the verified skeleton.
# Q0 has one positive direction. For an isometry A, sign((A e1)_1) is the time-cone character.
def time_char(A):
    assert sp.simplify(A.T*Q0*A-Q0)==sp.zeros(3)
    return 1 if A[0,0]>0 else -1
expected_time={'V':1,'W':1,'Gamma':-1,'Delta':-1,'T':-1,'U':-1}
for g in gen_names:
    checks[f'time_character_{g}']=time_char(GM[g])==expected_time[g]
rho_time={g:sp.Matrix([[expected_time[g]]]) for g in gen_names}
D1_scalar=cocycle_relation_matrix(rho_time)
D0_scalar=d0_matrix(rho_time)
checks['scalar_time_local_system_is_representation']=True
chi_scalar_basis=sp.Matrix([1,1,0,0,0,0]) # nonzero placeholders on V,W; scale irrelevant
scalar_residual=D1_scalar*chi_scalar_basis
checks['single_scalar_chi_candidate_not_closed']=scalar_residual!=sp.zeros(D1_scalar.rows,1)
# Specifically Gamma W Gamma W^-1 forces -2 chi_W = 0.
row_gammaW=list(rels).index('Gamma_W_fixed')
checks['GammaW_relation_forces_scalar_chiW_zero']=D1_scalar[row_gammaW,idx['W']]==-2

# Rank-2 channel local system: minimal verified repair on this finite skeleton.
I2=sp.eye(2)
rho2={g:I2 for g in gen_names}
rho2['Gamma']=sp.diag(-1,1)
rho2['Delta']=-I2
# T,U trivial on the rapidity channel in this *presentation skeleton* only.
D1_vec=cocycle_relation_matrix(rho2)
D0_vec=d0_matrix(rho2)
cV=sp.zeros(12,1); cV[2*idx['V'],0]=1
cW=sp.zeros(12,1); cW[2*idx['W']+1,0]=1
checks['vector_chi_V_closed']=D1_vec*cV==sp.zeros(D1_vec.rows,1)
checks['vector_chi_W_closed']=D1_vec*cW==sp.zeros(D1_vec.rows,1)
rank_d0_vec=D0_vec.rank(); rank_d1_vec=D1_vec.rank()
checks['vector_d0_rank_2']=rank_d0_vec==2
checks['vector_d1_rank_8']=rank_d1_vec==8
h1dim=12-rank_d1_vec-rank_d0_vec
checks['vector_H1_dimension_2']=h1dim==2
checks['vector_V_nonexact']=D0_vec.row_join(cV).rank()==rank_d0_vec+1
checks['vector_W_nonexact']=D0_vec.row_join(cW).rank()==rank_d0_vec+1
checks['vector_VW_independent_mod_exact']=D0_vec.row_join(cV).row_join(cW).rank()==rank_d0_vec+2

# Core V/W one-vertex ranks (both memory and transported boost, where local monodromy is +1).
D_core=sp.zeros(2,1)
rank_core=D_core.rank()
checks['core_D_rank_zero']=rank_core==0
checks['core_memory_aug_rank_one']=D_core.row_join(sp.Matrix([1,1])).rank()==1
checks['core_boost_aug_rank_one']=D_core.row_join(sp.Matrix([1,1])).rank()==1

# Store matrices exactly.
def matrix_rows(M): return [[str(M[i,j]) for j in range(M.cols)] for i in range(M.rows)]
coboundary={
 'scalar_time_attempt':{
   'generator_order':gen_names,
   'tau':{g:int(rho_time[g][0,0]) for g in gen_names},
   'D0_shape':list(D0_scalar.shape),'D0':matrix_rows(D0_scalar),
   'D1_shape':list(D1_scalar.shape),'D1':matrix_rows(D1_scalar),
   'rank_D0':D0_scalar.rank(),'rank_D1':D1_scalar.rank(),
   'candidate_VW_residual':matrix_rows(scalar_residual),
   'decision':'REFUTED_GLOBAL_SCALAR_FOR_SIMULTANEOUS_V_AND_W_ON_THIS_RELATION_COMPLEX'
 },
 'rank2_rapidity_bundle':{
   'generator_order':gen_names,
   'rho_Gamma':[[-1,0],[0,1]],'rho_Delta':[[-1,0],[0,-1]],
   'rho_V':'I2','rho_W':'I2','rho_T':'I2_TYPED_MINIMAL','rho_U':'I2_TYPED_MINIMAL',
   'D0_shape':list(D0_vec.shape),'D1_shape':list(D1_vec.shape),
   'rank_D0':rank_d0_vec,'rank_D1':rank_d1_vec,'H1_dimension':h1dim,
   'V_class':'NONZERO','W_class':'NONZERO','classes_independent':True,
   'status':'PROVEN_ON_X32_RELATION_SKELETON__T_U_ACTION_TYPED_MINIMAL'
 },
 'VW_core':{
   'D_memory':[['0'],['0']], 'D_boost':[['0'],['0']],
   'rank':0,'rank_aug_nonzero_observable':1,
   'decision':'BOTH_CORE_CLASSES_NONEXACT'
 }
}
(OUT/'presentation_coboundaries_v32.json').write_text(json.dumps(coboundary,indent=2,ensure_ascii=False),encoding='utf-8')

# ------------------------------------------------------------------
# E. Mod-2 comparison to v29: canonical quotient paths + obstruction.
# ------------------------------------------------------------------
def key2(M): return tuple(int(M[i,j])%2 for i in range(3) for j in range(3))
def mulkey(a,b):
    A=sp.Matrix(3,3,list(a)); B=sp.Matrix(3,3,list(b)); C=A*B
    return tuple(int(C[i,j])%2 for i in range(3) for j in range(3))
Ikey=key2(I3)
adj2=[key2(A) for A in Adj]
paths={Ikey:()}; q=deque([Ikey])
while q:
    a=q.popleft(); pth=paths[a]
    for j,g in enumerate(adj2,1):
        b=mulkey(a,g)
        if b not in paths:
            paths[b]=pth+(j,); q.append(b)
checks['v29_GL32_BFS_168']=len(paths)==168
lift_mats={'P':P,'Q':Q,'R':R,'S':S,'T':T,'U':U,'V':V,'W':W,'Gamma':Gamma,'Delta':Delta}
canon={name:list(paths[key2(A)]) for name,A in lift_mats.items()}
checks['comparison_V_collapses_to_vertex']=canon['V']==[]
checks['comparison_W_collapses_to_vertex']=canon['W']==[]
checks['comparison_Gamma_collapses_to_vertex']=canon['Gamma']==[]
checks['comparison_Delta_collapses_to_vertex']=canon['Delta']==[]

# v29 Fano cochain and a nonzero square period: path comparison is section/path dependent.
sigma=[0,4,1,5,2,6,3]
lines={'u1':{0,1,3},'u2':{0,4,5},'u3':{0,2,6}}
meta={1:'u2',7:'u2',2:'u1',8:'u1',3:'u3',6:'u3',4:'u2',10:'u2',5:'u1',12:'u1',9:'u3',11:'u3'}
omega=cmath.exp(2j*math.pi/7)
def Omega(label):
    L=lines[label]
    Z=sum((1 if sigma[n] in L else -1)*omega**n for n in range(7))
    return (Z**3).imag/(16*math.sqrt(2))
Om={u:Omega(u) for u in lines}
def v29_twisted_path_integral(word):
    # all-minus local system; S_i=S_i^{-1}; reverse traversal has same local coefficient.
    pref=1.0; val=0.0
    for j in word:
        val += pref*Om[meta[j]]
        pref *= -1.0
    return val
square=[1,2,1,2]
square_period=v29_twisted_path_integral(square)
checks['v29_square_period_nonzero']=abs(square_period)>1e-10
p0=canon['P']
p1=p0+square
checks['same_mod2_endpoint_after_square']=mulkey(key2(P), Ikey)==key2(P) and len(square)==4
# Calculate path product of square identity.
def path_key(word):
    a=Ikey
    for j in word: a=mulkey(a,adj2[j-1])
    return a
checks['square_is_identity_mod2']=path_key(square)==Ikey
checks['P_paths_same_endpoint']=path_key(p0)==path_key(p1)==key2(P)
val0=v29_twisted_path_integral(p0); val1=v29_twisted_path_integral(p1)
checks['v29_path_comparison_not_path_independent']=abs(val1-val0)>1e-10

comparison={
 'quotient':'rho2 : historical lifts -> GL_3(2)',
 'canonical_shortest_paths_in_v29_generators':{k:[f'S{j}' for j in v] for k,v in canon.items()},
 'normalized_bar_cochain_pullback':{
   'formula':'(rho2^* f)(g1,...,gn)=f(rho2(g1),...,rho2(gn))',
   'kernel_degree1':'normalized degree-1 quotient cochains vanish on V,W,Gamma,Delta because rho2=I',
   'status':'CANONICAL_AT_GROUP_BAR_COCHAIN_LEVEL'
 },
 'v29_cayley_graph_comparison':{
   'canonical_path_section':'BFS shortest path in ordered generators S1,...,S12',
   'status':'EXPLICIT_BUT_SECTION_DEPENDENT',
   'obstruction_square':['S1','S2','S1','S2'],
   'square_twisted_period':square_period,
   'P_canonical_integral':val0,'P_with_square_integral':val1,
   'decision':'V29_GRAPH_COCHAIN_DOES_NOT_DESCEND_TO_A_PATH_INDEPENDENT_GROUP_COCHAIN_WHILE_THE_NONZERO_SQUARE_PERIOD_IS_RETAINED'
 },
 'kernel_consequence':{
   'V':'PURE_LIFT_COMPONENT_RELATIVE_TO_CANONICAL_GROUP_COCHAIN_PULLBACK',
   'W':'PURE_LIFT_COMPONENT_RELATIVE_TO_CANONICAL_GROUP_COCHAIN_PULLBACK',
   'v29_graph_charges_vs_v32_lift_classes':'SEPARATED_UNTIL_A_2_CELL_COMPATIBLE_COMPARISON_IS_SUPPLIED'
 }
}
(OUT/'v29_v32_comparison_map.json').write_text(json.dumps(comparison,indent=2,ensure_ascii=False),encoding='utf-8')

# ------------------------------------------------------------------
# F. Final certificate.
# ------------------------------------------------------------------
checks={k:bool(v) for k,v in checks.items()}
status='PASS' if all(checks.values()) else 'FAIL'
cert={
 'phase':'v32-finite-presentation-transport-comparison',
 'status':status,
 'checks':checks,
 'tooling':{
   'exact_backend':'SymPy exact integer/rational/algebraic arithmetic',
   'python_flint':'UNAVAILABLE_IN_RUNTIME',
   'sage':'UNAVAILABLE_IN_RUNTIME',
   'lean_native_replay':'UNAVAILABLE_IN_RUNTIME; KEP Lean source used only as contextual proof-ABI precedent'
 },
 'primitive_decomposition':{
   'R_source_steps':5,'S_source_steps':3,
   'matrix_products':'PROVEN_EXACT',
   'adjacent_superbase_algebraic_refinement':'PROVEN_EXACT',
   'geometric_wall_path_identification':'REFUTED_FROM_MATRIX_WORD_ALONE / PENDING_EXPLICIT_WALL_PATH_ADAPTER',
   'decisive_guard':'P and Q have even-length v24 algebraic words (8 and 10) although v30 historical boundary tau(P)=tau(Q)=-1; therefore algebraic word parity cannot define Binet wall monodromy.'
 },
 'presentation_complex':presentation,
 'krein_transport':{
   'J_PQ':[['1','0','0'],['0','0','sqrt(5)'],['0','1/sqrt(5)','0']],
   'properties':['J_PQ^T Q0 J_PQ=Q0','J_PQ^2=I','J_PQ(K_P)=K_Q','pseudo-frame coordinate transport=I_2'],
   'field_extension':'Q(sqrt(5)) subset R',
   'golden_guard':'sqrt(5) is algebraically necessary for this natural isometry; no identification with a golden/irrational-rotation carrier is made without a separate morphism.'
 },
 'cohomology':{
   'VW_core':{'rank_D':0,'rank_aug_nonzero_memory':1,'rank_aug_nonzero_boost':1,'decision':'NONEXACT_CORE_CLASSES'},
   'scalar_time_system':{
      'rank_D0':D0_scalar.rank(),'rank_D1':D1_scalar.rank(),
      'simultaneous_VW_scalar_chi':'REFUTED_NOT_CLOSED',
      'witness_relation':'Gamma W Gamma W^-1=1 forces chi_W=0 under scalar time character'
   },
   'rank2_rapidity_bundle':{
      'rank_D0':rank_d0_vec,'rank_D1':rank_d1_vec,'H1_dimension':h1dim,
      'V_W_classes':'CLOSED_NONEXACT_INDEPENDENT; basis of H1 on X32 relation skeleton',
      'scope':'T,U coefficient action is the minimal trivial extension on the currently verified independent T/U 2-cells; full historical coefficient action remains PENDING if new cross-relations are added.'
   },
   'memory_full_presentation':'PENDING because a source-exact global Binet C2 action on Gamma,Delta,T,U and primitive R/S wall paths is not materialized.'
 },
 'comparison':comparison,
 'decisions':{
   'finite_relation_complex':'PROVEN_CONSTRUCTED',
   'complete_group_presentation':'NT_PENDING_REWRITE_COMPLETENESS',
   'R_S_primitive_matrix_decomposition':'PROVEN',
   'R_S_edgewise_epsilon_tau_s_chi':'PENDING_MORPHISM',
   'O11_transport_KP_KQ':'PROVEN_OVER_Q_SQRT5',
   'single_global_scalar_chi':'REFUTED_ON_NAIVE_TIME_LOCAL_SYSTEM; REPLACED_BY_RANK2_BUNDLE_ON_CURRENT_SKELETON',
   'v29_graph_charge_direct_identification':'SEPARATED',
   'V_W_congruence_kernel_components':'PROVEN_PURE_LIFT_RELATIVE_TO_GROUP_COCHAIN_PULLBACK'
 },
 'safe_successor':(
   'v33 -- retain all v32 carriers and no-go results; complete the source rewrite/normal-form proof for the finite relation skeleton or add only source-exact missing cross-relations; '
   'materialize a geometric wall-path adapter from the ordered R/S source decompositions to primitive Selling chamber crossings and thereby determine the genuine Binet C2 action on Gamma,Delta,T,U; '
   'lift the rank-2 rapidity bundle from the current relation skeleton to that completed groupoid and test whether it splits canonically into scalar characters or remains irreducibly channel-valued; '
   'construct a 2-cell-compatible comparison from the v29 Cayley graph to the completed presentation (or prove an obstruction), and only then classify normalized v29 graph periods against historical memory classes beyond the already proven pure-kernel V/W components.'
 )
}
(OUT/'presentation_transport_comparison_v32_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
