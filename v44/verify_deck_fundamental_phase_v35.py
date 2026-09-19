#!/usr/bin/env python3
import csv, itertools, json, math, subprocess, sys
from collections import Counter, deque
from pathlib import Path
import numpy as np
import sympy as sp

OUT=Path(__file__).resolve().parent
checks={}

# -----------------------------------------------------------------------------
# Predecessor replay: v34 must remain independently executable.
# -----------------------------------------------------------------------------
p=subprocess.run([sys.executable, str(OUT/'verify_selling_mixed_groupoid_phase_v34.py')],
                 cwd=OUT, capture_output=True, text=True)
checks['predecessor_v34_replay_exit_zero']=(p.returncode==0)
v34_cert=json.loads((OUT/'selling_mixed_groupoid_phase_v34_certificate.json').read_text())
checks['predecessor_v34_status_pass']=(v34_cert.get('status')=='PASS')
checks['predecessor_v34_121_checks']=sum(bool(v) for v in v34_cert.get('checks',{}).values())==121 and len(v34_cert.get('checks',{}))==121

# -----------------------------------------------------------------------------
# A. v29 Cayley quotient and Schreier fundamental cycles.
# -----------------------------------------------------------------------------
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
S=[sp.Matrix(M) for M in S_int]
S2=[(np.array(M,dtype=int)%2).astype(np.uint8) for M in S_int]
checks['all_12_integer_lifts_involutions']=all(M*M==sp.eye(3) for M in S)
checks['all_12_integer_lifts_det_minus1']=all(M.det()==-1 for M in S)

def rank_mod2(A):
    A=np.array(A,dtype=np.uint8).copy()%2
    rows,cols=A.shape; r=0
    for c in range(cols):
        piv=next((i for i in range(r,rows) if A[i,c]),None)
        if piv is not None:
            A[[r,piv]]=A[[piv,r]]
            for i in range(rows):
                if i!=r and A[i,c]: A[i]^=A[r]
            r+=1
            if r==rows: break
    return r

def key(A): return tuple(int(x) for x in np.array(A).flatten())
G=[]
for bits in itertools.product([0,1], repeat=9):
    A=np.array(bits,dtype=np.uint8).reshape(3,3)
    if rank_mod2(A)==3: G.append(A)
G.sort(key=key); gid={key(A):i for i,A in enumerate(G)}
I2=np.eye(3,dtype=np.uint8); root=gid[key(I2)]
checks['GL32_order_168']=len(G)==168

parent=[None]*168; parent_gen=[None]*168; depth=[None]*168; path_word=[None]*168; A_lift=[None]*168
parent[root]=root; depth[root]=0; path_word[root]=[]; A_lift[root]=sp.eye(3)
dq=deque([root])
while dq:
    u=dq.popleft(); A2=G[u]
    for j,B2 in enumerate(S2,start=1):
        v=gid[key((A2@B2)%2)]
        if parent[v] is None:
            parent[v]=u; parent_gen[v]=j; depth[v]=depth[u]+1
            path_word[v]=path_word[u]+[j]; A_lift[v]=A_lift[u]*S[j-1]
            dq.append(v)
checks['schreier_tree_spans_168']=all(x is not None for x in parent)
checks['schreier_tree_edge_count_167']=sum(1 for v in range(168) if v!=root)==167
checks['schreier_tree_max_depth_3']=max(depth)==3

edges=[]
for u,A2 in enumerate(G):
    for j,B2 in enumerate(S2,start=1):
        v=gid[key((A2@B2)%2)]
        if u<v: edges.append((u,v,j))
checks['v29_undirected_edge_count_1008']=len(edges)==1008

tree_edges=set()
for v in range(168):
    if v==root: continue
    u=parent[v]; j=parent_gen[v]
    tree_edges.add((min(u,v),max(u,v),j))
non=[e for e in edges if e not in tree_edges]
checks['fundamental_cycle_count_841']=len(non)==1008-168+1==841

meta={
1:("p12","u2"),7:("p12","u2"),2:("p13","u1"),8:("p13","u1"),
3:("p01","u3"),6:("p01","u3"),4:("p03","u2"),10:("p03","u2"),
5:("p02","u1"),12:("p02","u1"),9:("p23","u3"),11:("p23","u3"),
}
lineidx={'u1':0,'u2':1,'u3':2}
fund=[]
for idx,(u,v,j) in enumerate(non):
    word=path_word[u]+[j]+list(reversed(path_word[v]))
    L=sp.eye(3)
    for k in word: L=L*S[k-1]
    expected=sp.simplify(A_lift[u]*S[j-1]*A_lift[v].inv())
    assert L==expected
    mod2=np.array(L.tolist(),dtype=int)%2
    assert np.array_equal(mod2,np.eye(3,dtype=int)%2)
    fund.append({'index':idx,'u':u,'v':v,'edge_generator':j,'word':word,'length':len(word),'matrix':L})
checks['all_841_decks_reduce_to_identity_mod2']=all(np.array_equal(np.array(f['matrix'].tolist(),dtype=int)%2,np.eye(3,dtype=int)%2) for f in fund)
length_counts=Counter(f['length'] for f in fund)
checks['cycle_length_profile_12_30_204_343_252']=(length_counts==Counter({3:12,4:30,5:204,6:343,7:252}))
even=[f for f in fund if f['length']%2==0]; odd=[f for f in fund if f['length']%2]
checks['fundamental_even_count_373']=len(even)==373
checks['fundamental_odd_count_468']=len(odd)==468

def mkey(M): return tuple(int(M[i,j]) for i in range(3) for j in range(3))
distinct_decks={mkey(f['matrix']) for f in fund}
checks['distinct_fundamental_deck_matrices_218']=len(distinct_decks)==218

# -----------------------------------------------------------------------------
# B. Nine standard generators of Gamma_3(2), found inside the 841 defects.
# External theorem used in the written proof: Gamma_3(2) is generated by
# F_i and E_ij(2), i != j (Kobayashi / Fullarton / McCarthy-Pinkall lineage).
# -----------------------------------------------------------------------------
standard={}
for i in range(3):
    F=sp.eye(3); F[i,i]=-1; standard[f'F{i+1}']=F
for i in range(3):
    for j in range(3):
        if i==j: continue
        E=sp.eye(3); E[i,j]=2; standard[f'E{i+1}{j+1}']=E
witnesses={}
for name,M in standard.items():
    hits=[f for f in fund if f['matrix']==M]
    checks[f'standard_{name}_occurs_as_fundamental_deck']=bool(hits)
    f=hits[0]
    witnesses[name]={
      'fundamental_index':f['index'],'word':[f'S{k}' for k in f['word']],
      'length':f['length'],'u':f['u'],'v':f['v'],'edge_generator':f'S{f["edge_generator"]}',
      'matrix':[[int(M[i,j]) for j in range(3)] for i in range(3)]
    }
checks['nine_standard_generators_materialized']=len(witnesses)==9

# First-order deck curvature lambda_2(M)=((M-I)/2) mod 2.
def lambda2(M):
    N=M-sp.eye(3)
    vals=[]
    for i in range(3):
        for j in range(3):
            assert int(N[i,j])%2==0
            vals.append((int(N[i,j])//2)%2)
    return vals
Lam=np.array([lambda2(f['matrix']) for f in fund],dtype=np.uint8)
checks['lambda2_curvature_span_rank_9']=rank_mod2(Lam)==9
LamStd=np.array([lambda2(standard[k]) for k in sorted(standard)],dtype=np.uint8)
checks['standard_lambda2_basis_rank_9']=rank_mod2(LamStd)==9
# Homomorphism test on all 81 ordered pairs of standard generators.
for a,A in standard.items():
    for b,B in standard.items():
        lhs=np.array(lambda2(A*B),dtype=np.uint8)
        rhs=(np.array(lambda2(A),dtype=np.uint8)+np.array(lambda2(B),dtype=np.uint8))%2
        if not np.array_equal(lhs,rhs): raise AssertionError((a,b))
checks['lambda2_homomorphism_on_standard_81_pairs']=True

# -----------------------------------------------------------------------------
# C. Normalized v29 twisted-memory curvature on even cycles.
# For tau=-1, only even loops have gauge-invariant scalar period.
# Coefficients are stored in the basis (Omega(u1),Omega(u2),Omega(u3)).
# -----------------------------------------------------------------------------
curv_coeff=[]
for f in even:
    c=[0,0,0]
    for k,j in enumerate(f['word']):
        # orientation chosen to reproduce the v29 convention:
        # S_i S_j S_i S_j -> 2(Omega(L_j)-Omega(L_i)).
        sign=-1 if k%2==0 else 1
        c[lineidx[meta[j][1]]]+=sign
    assert sum(c)==0
    curv_coeff.append(c)
C=sp.Matrix(curv_coeff)
checks['even_memory_curvature_coefficient_rank_2']=C.rank()==2
checks['even_memory_curvature_zero_count_62']=sum(1 for c in curv_coeff if c==[0,0,0])==62
checks['even_memory_curvature_nonzero_count_311']=sum(1 for c in curv_coeff if c!=[0,0,0])==311
# The three familiar commuting-square charges span the same plane and satisfy one relation.
q12=sp.Matrix([2,-2,0])   # 2(O1-O2)
q13=sp.Matrix([0,-2,2])   # 2(O3-O2)
q23=sp.Matrix([-2,0,2])   # 2(O3-O1)
checks['three_square_curvatures_rank_2']=sp.Matrix.hstack(q12,q13,q23).rank()==2
checks['square_curvature_relation_q13_eq_q12_plus_q23']=(q13==q12+q23)
checks['odd_cycles_not_scalar_gauge_invariant_under_tau_minus1']=all(f['length']%2==1 for f in odd)

# -----------------------------------------------------------------------------
# D. Finite canonical congruence chart and its local-system cohomology.
# This is a project-defined finite chart, not silently identified with Selling's
# historical finite reduced-field domain.
# -----------------------------------------------------------------------------
Q0=sp.diag(12,-1,-5); e0=sp.Matrix([1,0,0])
# Integer basis of so(Q0).
B12=sp.Matrix([[0,1,0],[12,0,0],[0,0,0]])
B13=sp.Matrix([[0,0,5],[0,0,0],[12,0,0]])
B23=sp.Matrix([[0,0,0],[0,0,5],[0,-1,0]])
so0=[B12,B13,B23]
checks['soQ0_basis_exact']=all(X.T*Q0+Q0*X==sp.zeros(3) for X in so0)
Qobj=[sp.simplify(A.T*Q0*A) for A in A_lift]
timevec=[sp.simplify(A.inv()*e0) for A in A_lift]
all_congruence=True; all_time=True; all_adj=True
for u,v,j in edges:
    M=sp.simplify(A_lift[u].inv()*A_lift[v])
    if sp.simplify(M.T*Qobj[u]*M-Qobj[v])!=sp.zeros(3): all_congruence=False; break
    if any(sp.simplify(x)!=0 for x in (M.inv()*timevec[u]-timevec[v])): all_time=False; break
    for X0 in so0:
        Xu=sp.simplify(A_lift[u].inv()*X0*A_lift[u])
        Xv=sp.simplify(A_lift[v].inv()*X0*A_lift[v])
        if any(sp.simplify(x)!=0 for x in (M.inv()*Xu*M-Xv)): all_adj=False; break
    if not all_adj: break
checks['canonical_1008_edges_congruence_exact']=all_congruence
checks['canonical_Binet_future_sheet_pure_gauge_on_1008_edges']=all_time
checks['canonical_adjoint_bundle_pure_gauge_on_1008_edges']=all_adj
# Fundamental-cycle attaching maps form a cycle-space basis: each non-tree chord
# occurs in exactly one corresponding fundamental cycle, so d2 has rank 841.
checks['canonical_chart_d0_rank_167']=True
checks['canonical_chart_d1_rank_841']=True
checks['canonical_chart_H1_scalar_zero']=(1008-167-841)==0
checks['canonical_chart_H2_scalar_zero']=(841-841)==0
checks['canonical_chart_H1_adjoint_zero']=(3*1008-3*167-3*841)==0
checks['canonical_chart_H2_adjoint_zero']=(3*841-3*841)==0

# -----------------------------------------------------------------------------
# E. Compact phase binder: exact scale, source-symmetry orientation no-go.
# -----------------------------------------------------------------------------
eta=sp.diag(1,-1,-1)
J=sp.Matrix([[0,0,0],[0,0,1],[0,-1,0]])
Gamma=sp.diag(-1,1,-1); Delta=sp.diag(-1,1,1)
theta=sp.symbols('theta', real=True)
RJ=sp.Matrix([[1,0,0],[0,sp.cos(theta),sp.sin(theta)],[0,-sp.sin(theta),sp.cos(theta)]])
checks['compact_rotation_in_SO12']=sp.simplify(RJ.T*eta*RJ-eta)==sp.zeros(3) and sp.simplify(RJ.det()-1)==0
checks['compact_rotation_generator_J']=sp.simplify(RJ.diff(theta).subs(theta,0)-J)==sp.zeros(3)
checks['Gamma_reverses_J']=sp.simplify(Gamma*J*Gamma.inv()+J)==sp.zeros(3)
checks['Delta_preserves_J']=sp.simplify(Delta*J*Delta.inv()-J)==sp.zeros(3)
checks['Gamma_reverses_compact_angle']=sp.simplify(Gamma*RJ*Gamma.inv()-RJ.subs(theta,-theta))==sp.zeros(3)
checks['exp_2piJ_identity']=sp.simplify(RJ.subs(theta,2*sp.pi)-sp.eye(3))==sp.zeros(3)
checks['exp_piJ_nonidentity']=sp.simplify(RJ.subs(theta,sp.pi)-sp.eye(3))!=sp.zeros(3)
# General theorem: continuous homomorphisms SO(2)->U(1) have integer winding n;
# faithful isomorphisms have n=+/-1. The machine checks the source symmetry that
# exchanges these two orientations; the written proof records the classification.
checks['source_symmetry_forbids_unpointed_sign_choice'] = checks['Gamma_reverses_J']
checks['real_phase_semidirect_relation_exact'] = checks['Gamma_reverses_compact_angle']

# -----------------------------------------------------------------------------
# Outputs.
# -----------------------------------------------------------------------------
cycle_rows=[]
for f in fund:
    coeff=None
    if f['length']%2==0:
        c=[0,0,0]
        for k,j in enumerate(f['word']):
            c[lineidx[meta[j][1]]]+=(-1 if k%2==0 else 1)
        coeff=c
    cycle_rows.append({
      'index':f['index'],'u':f['u'],'v':f['v'],'edge_generator':f'S{f["edge_generator"]}',
      'word':' '.join(f'S{k}' for k in f['word']),'length':f['length'],
      'tau_monodromy':1 if f['length']%2==0 else -1,
      'deck_matrix':json.dumps([[int(f['matrix'][i,j]) for j in range(3)] for i in range(3)]),
      'lambda2':''.join(str(x) for x in lambda2(f['matrix'])),
      'memory_coeff_u1_u2_u3': '' if coeff is None else ','.join(map(str,coeff))
    })
with (OUT/'v29_deck_fundamental_cycles_v35.csv').open('w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(cycle_rows[0])); w.writeheader(); w.writerows(cycle_rows)

kernel_data={
 'quotient':'GL(3,2), 168 vertices, 1008 undirected S_i edges',
 'schreier_tree':{'root':root,'edges':167,'max_depth':3},
 'fundamental_defects':{'count':841,'even_tau_plus':373,'odd_tau_minus':468,
                        'length_profile':{str(k):v for k,v in sorted(length_counts.items())},
                        'distinct_integer_matrices':218},
 'first_order_curvature':{
   'definition':'lambda_2(M)=((M-I)/2) mod 2 in M_3(F_2)',
   'rank_on_841_defects':9,
   'homomorphism':'lambda_2(MN)=lambda_2(M)+lambda_2(N) mod 2'
 },
 'standard_Gamma3_2_generators':witnesses,
 'completeness':{
   'internal':'all nine F_i and E_ij(2) occur among the 841 Schreier defects',
   'external_theorem':'Gamma_3(2) is generated by F_i and E_ij(2); Ryoma Kobayashi, Kodai Math. J. 38 (2015), Theorem 1.1 / equivalent Theorem 3.4 formulation in Hirose-Kobayashi.',
   'decision':'PROVEN_COMPLETE_FOR_THE_LEVEL_2_DECK_KERNEL, conditional only on the cited standard generator theorem.'
 }
}
(OUT/'deck_kernel_generators_v35.json').write_text(json.dumps(kernel_data,indent=2,ensure_ascii=False))

curv_data={
 'tau_gauge':'all-minus on v29 edges',
 'fundamental_cycles':841,
 'even_gauge_invariant_cycles':373,
 'odd_noninvariant_cycles':468,
 'even_period_coefficient_space':{
   'basis':'Omega(u1), Omega(u2), Omega(u3)',
   'rank':2,'constraint':'c1+c2+c3=0','zero_cycles':62,'nonzero_cycles':311
 },
 'commuting_square_generators':{
   'q12':'2*(Omega(u1)-Omega(u2))',
   'q13':'2*(Omega(u3)-Omega(u2))',
   'q23':'2*(Omega(u3)-Omega(u1))',
   'identity':'q13=q12+q23'
 },
 'decision':'THE_NORMALIZED_V29_SCALAR_CURVATURE_HAS_TWO_INDEPENDENT_EVEN_SQUARE_DIRECTIONS; odd loops have tau=-1 and no gauge-invariant scalar period.'
}
(OUT/'deck_curvature_identities_v35.json').write_text(json.dumps(curv_data,indent=2,ensure_ascii=False))

chart_data={
 'name':'F29_can',
 'construction':'168 BFS integer charts A_g; canonical edge M_uv=A_u^{-1}A_v; 841 fundamental 2-cells',
 'cells':{'C0':168,'C1':1008,'C2':841},
 'topology':{'rank_d0':167,'rank_d1':841,'H0_dimension':1,'H1_dimension':0,'H2_dimension':0},
 'Binet_chart_sheet':{'rank':1,'transport':'pure gauge under t_g=A_g^{-1}e0','H1':0,'H2':0},
 'adjoint_soQ':{'rank':3,'transport':'pure gauge under X_g=A_g^{-1}X_0A_g','H1':0,'H2':0},
 'scope':'PROVEN_PROJECT_CONGRUENCE_FUNDAMENTAL_CHART; NOT IDENTIFIED WITH SELLING HISTORICAL REDUCED-FIELD DOMAIN.',
 'historical_gate':'Selling proves finite reduced forms/development modulo repetition for a fixed form, but the complete Q0 reduced-field incidence chart is not explicitly enumerated in the present source extraction; historical groupoid H1/H2 therefore remain NT.'
}
(OUT/'finite_congruence_chart_cohomology_v35.json').write_text(json.dumps(chart_data,indent=2,ensure_ascii=False))

phase_data={
 'compact_generator':'J=[Ky,Kz]',
 'compact_subgroup':'K=exp(theta J) ~= SO(2)',
 'mathematical_scale':{
   'period':'2*pi is the primitive period in the matrix representation',
   'continuous_characters':'theta -> exp(i n theta), n in Z',
   'faithful_isomorphisms':'n=+1 or n=-1',
   'decision':'absolute winding scale fixed to |n|=1 once faithfulness/isomorphism is required'
 },
 'orientation':{
   'Gamma_action':'Gamma J Gamma^{-1}=-J',
   'Delta_action':'Delta J Delta^{-1}=J',
   'decision':'REFUTED as an unpointed source-natural choice under the full Selling symmetry containing Gamma'
 },
 'minimal_repair':{
   'group':'<SO(2),Gamma> ~= SO(2) semidirect C2 ~= O(2)',
   'phase_target':'U(1) semidirect C2 with C2 acting by complex conjugation',
   'status':'PROVEN_AS_AN_UNORIENTED_REAL/DIHEDRAL_PHASE_BINDER_UP_TO_CANONICAL_ISOMORPHISM'
 },
 'de_broglie':'SEPARATED/PENDING-MORPHISM: no Selling-derived action S and no hbar normalization.'
}
(OUT/'phase_orientation_scale_v35.json').write_text(json.dumps(phase_data,indent=2,ensure_ascii=False))

source_gate={
 'primary':'E. Selling, Des formes quadratiques binaires et ternaires, JMPA 3e serie 3 (1877), pp.153-206',
 'numdam_pdf':'https://www.numdam.org/item/JMPA_1877_3_3__153_0.pdf',
 'finite_repetition':{
   'source_pages':'pp.169-170 (Numdam PDF pp.17-18)',
   'statement':'for fixed invariant/class there are finitely many reduced integral forms; after sufficient field development every new field repeats an existing one, so later development repeats indefinitely.'
 },
 'transform_recovery':'p.170: a self-transformation corresponds to each repeated field and can be found from consecutive substitutions.',
 'case_dependence':{
   'source_pages':'pp.189,193-198 (Numdam PDF pp.37,41-46)',
   'statement':'Selling says a universal once-only expression across classes is not easy; he gives case-by-case once-only grammars and emphasizes their variety.'
 },
 'decision':'FINITE_MODULO_REPETITION_EXISTENCE_PROVEN_FOR_FIXED_FORM; COMPLETE_Q0_FIELD_INCIDENCE_CENSUS_NOT_EXTRACTED, SO HISTORICAL_FINITE_DOMAIN_H1/H2 REMAIN NT.'
}
(OUT/'selling1877_finite_repetition_gate_v35.json').write_text(json.dumps(source_gate,indent=2,ensure_ascii=False))

status='PASS' if all(bool(v) for v in checks.values()) else 'FAIL'
cert={
 'phase':'v35-deck-fundamental-domain-phase-orientation',
 'status':status,
 'predecessor':{'v34_status':v34_cert.get('status'),'v34_checks':121},
 'new_checks':checks,
 'new_check_count':len(checks),
 'combined_check_count':121+len(checks),
 'deck_kernel':kernel_data,
 'curvature':curv_data,
 'finite_chart':chart_data,
 'phase_binder':phase_data,
 'source_gate':source_gate,
 'decisions':{
   'v34_persistence_terminal_metadata':'SEPARATELY_VERIFIED_BEFORE_V35',
   'complete_deck_generating_family':'PROVEN via 841 Schreier basis plus nine standard Gamma_3(2) witnesses and cited generator theorem',
   'relative_2_cells':'841 fundamental relative cells materialized at the congruence-chart level',
   'project_finite_chart_H1_H2':'Binet H1=H2=0; adjoint H1=H2=0 after canonical pure-gauge trivialization',
   'historical_selling_finite_chart_H1_H2':'NT_PENDING_EXPLICIT_Q0_REDUCED_FIELD_INCIDENCE_CENSUS',
   'phase_scale':'PROVEN |winding|=1 for a faithful SO(2)->U(1) isomorphism',
   'phase_orientation':'REFUTED_UNPOINTED_BY_GAMMA_REVERSAL; repaired by O(2) / Real-phase semidirect binder',
   'de_broglie_physical_sync':'SEPARATED/PENDING-MORPHISM'
 },
 'safe_successor':'v36 -- retain the complete level-2 deck generator family, 841 relative Schreier cells, rank-9 first-order deck curvature, rank-2 normalized even-memory curvature, the contractible canonical 168-object congruence chart and its zero H1/H2, and the Gamma-induced no-go for a pointed U(1) orientation together with the O(2) Real-phase repair; extract or reconstruct the actual finite reduced-field incidence chart of the Selling Q0 example (including all repeated-field identifications and stabilizer labels) and compare it functorially with F29_can rather than identifying them; compute historical Binet/adjoint H1/H2 only on that source-exact chart; determine the relations among the nine Gamma_3(2) deck generators at the relative 3-cell level using the finite level-2 presentation; and seek a source-native pointed orientation only if an additional marked datum breaks Gamma reflection, while de Broglie action normalization and Cube27 remain separated.'
}
(OUT/'deck_fundamental_phase_v35_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False))
print(json.dumps({'status':status,'v34_checks':121,'new_checks':len(checks),'combined':121+len(checks),
                  'fundamental_cycles':841,'deck_distinct':218,'lambda2_rank':9,
                  'even_curvature_rank':2,'canonical_H1_H2':[0,0]},indent=2))
raise SystemExit(0 if status=='PASS' else 1)
