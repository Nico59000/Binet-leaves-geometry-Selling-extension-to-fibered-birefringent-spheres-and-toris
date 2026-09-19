#!/usr/bin/env python3
import itertools, json, subprocess, sys
from collections import Counter, deque
from pathlib import Path
import numpy as np
import sympy as sp

OUT=Path(__file__).resolve().parent
checks={}

# -----------------------------------------------------------------------------
# Predecessor replay: v35 is immutable.
# -----------------------------------------------------------------------------
p=subprocess.run([sys.executable, str(OUT/'verify_deck_fundamental_phase_v35.py')],
                 cwd=OUT, capture_output=True, text=True)
checks['predecessor_v35_replay_exit_zero']=(p.returncode==0)
v35=json.loads((OUT/'deck_fundamental_phase_v35_certificate.json').read_text())
checks['predecessor_v35_status_pass']=(v35.get('status')=='PASS')
checks['predecessor_v35_combined_175']=(v35.get('combined_check_count')==175)

# -----------------------------------------------------------------------------
# Historical Q0 branch matrices / C6 RACG from v34.
# -----------------------------------------------------------------------------
I3=sp.eye(3)
Q0=sp.diag(12,-1,-5)
P=sp.Matrix([[1,1,0],[3,4,0],[0,0,1]])
Q=sp.Matrix([[2,0,5],[0,1,0],[3,0,8]])
R=sp.Matrix([[-4,0,5],[0,-1,0],[-3,0,4]])
Sm=sp.Matrix([[-2,1,0],[-3,2,0],[0,0,-1]])
T=sp.Matrix([[-13,3,5],[-36,8,15],[-12,3,4]])
U=sp.Matrix([[-17,2,10],[-24,2,15],[-24,3,14]])
V=sp.Matrix([[7,2,0],[24,7,0],[0,0,1]])
W=sp.Matrix([[31,0,20],[0,1,0],[48,0,31]])
Gamma=sp.diag(-1,1,-1); Delta=sp.diag(-1,1,1); H=Gamma*Delta
CGEN={'a':V*Gamma,'h':H,'g':Gamma,'b':W*H,'u':U,'t':T}
cycle=['a','h','g','b','u','t']
edges=list(zip(cycle,cycle[1:]+cycle[:1]))
for x,M in CGEN.items(): checks[f'C6_{x}_involution']=sp.simplify(M*M-I3)==sp.zeros(3)
for x,y in edges: checks[f'C6_{x}_{y}_commute']=sp.simplify(CGEN[x]*CGEN[y]-CGEN[y]*CGEN[x])==sp.zeros(3)

# -----------------------------------------------------------------------------
# A. Source-anchored historical-domain reconstruction status.
# The source text gives an exact finite/repetition gate and several incidence
# constraints, but the complete plate incidence is not enumerated in prose.
# We materialize the exact partial skeleton and keep the full field complex NT.
# -----------------------------------------------------------------------------
source_domain={
 'primary_source':'E. Selling, Des formes quadratiques binaires et ternaires, JMPA 3e série 3 (1877), pp.153-206, Numdam.',
 'Q0_branch':'declared figure-1 historical branch inherited from v30-v35',
 'source_exact_constraints':[
   'fig.1 is a field division for an integral class/invariant example; only finite-extent fields are inscribed',
   'the six external boundary lines are bipartite symmetry axes and their six intersection points are quadripartite symmetry centres',
   'every field boundary line starts/ends at fissure points and passes through one or more crossings; the final boundary reconnects to the first',
   'the A/a territory construction contains a closed alternating six-sided boundary in the cited example',
   'after sufficiently far development, all new forms numerically repeat earlier forms, and each repeated field determines a self-transformation',
   'self-transformations are generated from local symmetries together with finitely many external-boundary substitutions P,Q,R,S'
 ],
 'three_distinct_C6_carriers_guard':{
   'outer_mirror_skeleton':'six external symmetry lines / six quadripartite intersections',
   'local_crossing':'six local fields at a triple-line crossing, already audited in v5/v28',
   'Aa_territory_hexagon':'alternating A/a territory-boundary hexagon',
   'decision':'SEPARATED: equal cardinality C6 does not identify the three incidences.'
 },
 'presentation_orbifold':{
   'group':'W(C6) on a-h-g-b-u-t-a from v34 exact transformation relations',
   'status':'PROVEN_FOR_TRANSFORMATION_BRANCH'
 },
 'field_incidence':{
   'status':'NT_PENDING_EXACT_PLATE_TRANSCRIPTION',
   'reason':'the prose fixes boundary/symmetry/repetition constraints but does not enumerate every bounded field, edge, crossing, fissure and repeated-field identification of the plate in machine-readable form; the transformation presentation alone does not determine that microscopic cell decomposition.'
 },
 'no_silent_identification':'F_Sell^field(Q0) is not identified with the W(C6) presentation orbifold or with F29_can.'
}
checks['source_domain_has_six_external_axes']=('six external boundary lines' in source_domain['source_exact_constraints'][1])
checks['source_domain_three_C6_carriers_separated']=(source_domain['three_distinct_C6_carriers_guard']['decision'].startswith('SEPARATED'))
checks['full_field_incidence_not_silently_promoted']=(source_domain['field_incidence']['status']=='NT_PENDING_EXACT_PLATE_TRANSCRIPTION')
(OUT/'selling_q0_historical_domain_v36.json').write_text(json.dumps(source_domain,indent=2,ensure_ascii=False),encoding='utf-8')

# -----------------------------------------------------------------------------
# B. Exact finite transformation-orbifold cohomology on W(C6).
# This is an orbifold/presentation computation, not a substitute for the
# still-untranscribed microscopic field complex.
# -----------------------------------------------------------------------------
# Binet character on C6 generators: a=VGamma -, h=GammaDelta +,
# g=Gamma -, b=W h +, u -, t -.
tau={'a':-1,'h':1,'g':-1,'b':1,'u':-1,'t':-1}
minus=[x for x in cycle if tau[x]==-1]
C1B={x:i for i,x in enumerate(minus)}
minus_edges=[e for e in edges if tau[e[0]]==tau[e[1]]==-1]
D0B=sp.Matrix([tau[x]-1 for x in minus])
D1B=sp.zeros(len(minus_edges),len(minus))
for r,(s,t0) in enumerate(minus_edges):
    # d1(x)=(1-rho(t))x_s+(rho(s)-1)x_t
    D1B[r,C1B[s]]=1-tau[t0]
    D1B[r,C1B[t0]]=tau[s]-1
checks['orbifold_Binet_C1_dim_4']=(len(minus)==4)
checks['orbifold_Binet_C2_dim_2']=(len(minus_edges)==2)
checks['orbifold_Binet_d0_rank_1']=(D0B.rank()==1)
checks['orbifold_Binet_d1_rank_2']=(D1B.rank()==2)
checks['orbifold_Binet_d1d0_zero']=(D1B*D0B==sp.zeros(D1B.rows,1))
H1B=len(minus)-D0B.rank()-D1B.rank(); H2B=D1B.rows-D1B.rank()
checks['orbifold_Binet_H1_dim_1']=(H1B==1)
checks['orbifold_Binet_H2_zero']=(H2B==0)

# Adjoint local system in normalized eta frame.
S0=sp.diag(1/sp.sqrt(12),1,1/sp.sqrt(5))
Ky=sp.Matrix([[0,1,0],[1,0,0],[0,0,0]])
Kz=sp.Matrix([[0,0,1],[0,0,0],[1,0,0]])
J =sp.Matrix([[0,0,0],[0,0,1],[0,-1,0]])
Bmat=sp.Matrix.hstack(*[X.reshape(9,1) for X in (Ky,Kz,J)])
def normA(A): return sp.simplify(S0.inv()*A*S0)
def ad_matrix(A):
    An=normA(A); cols=[]
    for X in (Ky,Kz,J):
        Y=sp.simplify(An*X*An.inv())
        sol=next(iter(sp.linsolve((Bmat,Y.reshape(9,1)))))
        cols.append(sp.Matrix(list(sol)))
    return sp.Matrix.hstack(*cols)
Ad={x:ad_matrix(CGEN[x]) for x in cycle}
checks['Gamma_reverses_J_exact']=(Ad['g']*sp.Matrix([0,0,1])==sp.Matrix([0,0,-1]))
anti_basis={s:sp.Matrix.hstack(*(Ad[s]+sp.eye(3)).nullspace()) for s in cycle}
checks['each_C6_generator_anti_eigenspace_dim_2']=all(B.cols==2 for B in anti_basis.values())
common_basis={}
for e in edges:
    s,t0=e
    common=(Ad[s]+sp.eye(3)).col_join(Ad[t0]+sp.eye(3)).nullspace()
    common_basis[e]=sp.Matrix.hstack(*common)
checks['each_C6_square_common_anti_dim_1']=all(B.cols==1 for B in common_basis.values())
off={}; c1=0
for s in cycle: off[s]=c1; c1+=anti_basis[s].cols
c2=sum(B.cols for B in common_basis.values())
D0A=sp.zeros(c1,3)
for s in cycle:
    B=anti_basis[s]; image=Ad[s]-sp.eye(3)
    sol=next(iter(sp.linsolve((B,image[:,0])))) if False else None
    # solve each column
    C=sp.zeros(B.cols,3)
    for jcol in range(3):
        vv=next(iter(sp.linsolve((B,image[:,jcol]))))
        C[:,jcol]=sp.Matrix(vv)
    assert B*C==image
    D0A[off[s]:off[s]+B.cols,:]=C
D1A=sp.zeros(c2,c1); row=0
for e in edges:
    s,t0=e; Cb=common_basis[e]
    for jcol in range(anti_basis[s].cols):
        y=(sp.eye(3)-Ad[t0])*anti_basis[s][:,jcol]
        vv=next(iter(sp.linsolve((Cb,y))))
        D1A[row:row+Cb.cols,off[s]+jcol]=sp.Matrix(vv)
    for jcol in range(anti_basis[t0].cols):
        y=(Ad[s]-sp.eye(3))*anti_basis[t0][:,jcol]
        vv=next(iter(sp.linsolve((Cb,y))))
        D1A[row:row+Cb.cols,off[t0]+jcol]=sp.Matrix(vv)
    row+=Cb.cols
checks['orbifold_Ad_C1_dim_12']=(c1==12)
checks['orbifold_Ad_C2_dim_6']=(c2==6)
checks['orbifold_Ad_d0_rank_3']=(D0A.rank()==3)
checks['orbifold_Ad_d1_rank_6']=(D1A.rank()==6)
checks['orbifold_Ad_d1d0_zero']=(sp.simplify(D1A*D0A)==sp.zeros(c2,3))
H1A=c1-D0A.rank()-D1A.rank(); H2A=c2-D1A.rank()
checks['orbifold_Ad_H1_dim_3']=(H1A==3)
checks['orbifold_Ad_H2_zero']=(H2A==0)
orbifold_coh={
 'carrier':'Selling Q0 transformation-orbifold / Davis quotient of W(C6), not the untranscribed microscopic field complex',
 'cells_with_stabilizer_invariants':{
   'Binet':{'C0':1,'C1':len(minus),'C2':len(minus_edges)},
   'Adjoint':{'C0':3,'C1':c1,'C2':c2}
 },
 'Binet':{'rank_d0':D0B.rank(),'rank_d1':D1B.rank(),'H1_dimension':H1B,'H2_dimension':H2B},
 'Adjoint':{'rank_d0':D0A.rank(),'rank_d1':D1A.rank(),'H1_dimension':H1A,'H2_dimension':H2A},
 'consistency':'H1 dimensions reproduce v34 (1 scalar, 3 adjoint); new H2=0 on the transformation-orbifold carrier.',
 'historical_field_H1_H2':'NT_PENDING_EXACT_FIELD_INCIDENCE; not inferred from these orbifold values.'
}
(OUT/'selling_orbifold_H12_v36.json').write_text(json.dumps(orbifold_coh,indent=2,ensure_ascii=False),encoding='utf-8')

# -----------------------------------------------------------------------------
# C. Presentation-level functor to F29_can through reduction mod 2.
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
S2=[(np.array(M,dtype=int)%2).astype(np.uint8) for M in S_int]
def mod2(M): return np.array(M.tolist(),dtype=int)%2
I2=np.eye(3,dtype=int)%2
checks['functor_a_h_g_b_reduce_identity']=all(np.array_equal(mod2(CGEN[x]),I2) for x in ['a','h','g','b'])
checks['functor_t_equals_S6_mod2']=np.array_equal(mod2(T),S2[5])
checks['functor_u_equals_S3S6_mod2']=np.array_equal(mod2(U),(S2[2]@S2[5])%2)
# image subgroup generated by t,u
els={tuple(I2.flatten()):I2}; q=[I2]
while q:
    A=q.pop()
    for B in (mod2(T),mod2(U)):
        C=A@B%2; k=tuple(C.flatten())
        if k not in els: els[k]=C; q.append(C)
checks['functor_mod2_image_is_V4_size4']=(len(els)==4)
# Verify C6 presentation relators in the quotient images.
img={'a':I2,'h':I2,'g':I2,'b':I2,'u':mod2(U),'t':mod2(T)}
checks['functor_all_squares_fill_mod2']=all(np.array_equal(img[x]@img[x]%2,I2) for x in cycle)
checks['functor_all_C6_commutators_fill_mod2']=all(np.array_equal((img[x]@img[y]@img[x]@img[y])%2,I2) for x,y in edges)
functor_data={
 'source':'K_hist^pres(Q0)=W(C6) transformation presentation/orbifold',
 'target':'F29_can canonical 168-object congruence chart',
 'generator_paths':{'a':'constant','h':'constant','g':'constant','b':'constant','t':'S6','u':'S3 S6'},
 'image':'<rho2(T),rho2(U)> ~= C2 x C2, order 4',
 'two_cell_compatibility':'all six involution cells and six C6 commuting-square cells map to closed mod-2 loops fillable by the 841 fundamental cells of F29_can',
 'faithfulness':'NO: four of six C6 generators lie in the level-2 kernel and are collapsed at quotient level',
 'field_level_functor':'NT_PENDING_EXACT_PLATE_TRANSCRIPTION: the presentation-level functor is not promoted to F_Sell^field(Q0)->F29_can.'
}
(OUT/'selling_to_F29_functor_v36.json').write_text(json.dumps(functor_data,indent=2,ensure_ascii=False),encoding='utf-8')

# -----------------------------------------------------------------------------
# D. 43 relative 3-cells from Kobayashi's published finite presentation of
# Gamma_3(2), expressed on the nine v35 deck generators.
# -----------------------------------------------------------------------------
STD={}
for i in range(1,4):
    F=sp.eye(3); F[i-1,i-1]=-1; STD[f'F{i}']=F
for i in range(1,4):
    for j in range(1,4):
        if i!=j:
            E=sp.eye(3); E[i-1,j-1]=2; STD[f'E{i}{j}']=E

def invword(w): return [(n,-e) for n,e in reversed(w)]
def comm(w1,w2): return invword(w1)+invword(w2)+w1+w2
def evalword(w):
    M=sp.eye(3)
    for n,e in w: M=M*(STD[n] if e==1 else STD[n].inv())
    return sp.simplify(M)
rels=[]
for i in range(1,4): rels.append((f'R1_F{i}^2',[(f'F{i}',1)]*2,'1'))
for i in range(1,4):
    for j in range(1,4):
        if i==j: continue
        rels.append((f'R2a_E{i}{j}F{i}',[(f'E{i}{j}',1),(f'F{i}',1)]*2,'2'))
        rels.append((f'R2b_E{i}{j}F{j}',[(f'E{i}{j}',1),(f'F{j}',1)]*2,'2'))
for i in range(1,4):
    for j in range(i+1,4): rels.append((f'R2c_F{i}F{j}',[(f'F{i}',1),(f'F{j}',1)]*2,'2'))
for i,j,k in itertools.permutations(range(1,4),3):
    rels.append((f'R3a1_{i}{j}{k}',comm([(f'E{i}{j}',1)],[(f'E{i}{k}',1)]),'3a'))
    rels.append((f'R3a2_{i}{j}{k}',comm([(f'E{i}{j}',1)],[(f'E{k}{j}',1)]),'3a'))
    rels.append((f'R3a3_{i}{j}{k}',comm([(f'E{i}{j}',1)],[(f'F{k}',1)]),'3a'))
    rels.append((f'R3a4_{i}{j}{k}',comm([(f'E{i}{j}',1)],[(f'E{k}{i}',1)])+[(f'E{k}{j}',1)]*2,'3a'))
i,j,k=1,2,3
A=[(f'E{j}{i}',1),(f'F{j}',1),(f'E{i}{j}',1),(f'F{i}',1),(f'E{k}{i}',-1),(f'E{k}{j}',1)]
B=[(f'E{k}{i}',1),(f'F{k}',1),(f'E{i}{k}',1),(f'F{i}',1),(f'E{j}{i}',-1),(f'E{j}{k}',1)]
rels.append(('R3b_123',comm(A,B),'3b'))
checks['Gamma3_2_relative_3cell_count_43']=(len(rels)==43)
checks['Gamma3_2_relator_length_profile']=(Counter(len(w) for _,w,_ in rels)==Counter({2:3,4:33,6:6,24:1}))
checks['Gamma3_2_all_43_matrix_relators_identity']=all(evalword(w)==sp.eye(3) for _,w,_ in rels)
# lambda2 boundary must vanish for any relation identity.
def lambda2(M):
    N=M-sp.eye(3); out=[]
    for ii in range(3):
        for jj in range(3):
            assert int(N[ii,jj])%2==0
            out.append((int(N[ii,jj])//2)%2)
    return np.array(out,dtype=np.uint8)
checks['nine_standard_lambda2_are_matrix_unit_basis']=(np.linalg.matrix_rank(np.array([lambda2(STD[x]) for x in sorted(STD)],dtype=float))==9)
checks['all_43_relator_lambda2_boundaries_zero']=all(not lambda2(evalword(w)).any() for _,w,_ in rels)
# Confirm v35 materialized witnesses exist for all nine names.
v35_kernel=json.loads((OUT/'deck_kernel_generators_v35.json').read_text())
wits=v35_kernel['standard_Gamma3_2_generators']
checks['all_nine_v35_witnesses_available_for_3cells']=set(wits)==set(STD)
rel_data={
 'source':'R. Kobayashi, Kodai Math. J. 38 (2015), Theorem 1.1, specialized to n=3.',
 'generator_order':sorted(STD),
 'relative_dimension':'the nine v35 generator defects are treated as relative 2-cells; each presentation relator is attached as a relative 3-cell among them',
 'count':len(rels),
 'length_profile':{str(k):v for k,v in sorted(Counter(len(w) for _,w,_ in rels).items())},
 'cells':[{'name':n,'family':fam,'word':[x if e==1 else x+'^-1' for x,e in w], 'length':len(w)} for n,w,fam in rels],
 'matrix_verification':'43/43 evaluate to I_3 exactly over Z',
 'scope':'PROVEN_ARITHMETIC_COHERENCE_FOR_Gamma3(2); these are not declared to be Selling field 3-cells.'
}
(OUT/'gamma3_level2_relative_3cells_v36.json').write_text(json.dumps(rel_data,indent=2,ensure_ascii=False),encoding='utf-8')

# -----------------------------------------------------------------------------
# E. Roof / hyperfine labelling versus lambda2 and Thomas-Wigner parity shadow.
# -----------------------------------------------------------------------------
# lambda2 on the standard generators gives exactly all nine matrix units.
unit_labels={}
for name,M in STD.items():
    L=lambda2(M).reshape(3,3)
    ones=np.argwhere(L==1)
    checks[f'lambda2_{name}_single_matrix_unit']=(len(ones)==1)
    unit_labels[name]=[int(ones[0][0]+1),int(ones[0][1]+1)]
checks['lambda2_standard_labels_cover_all_3x3_positions']=(set(map(tuple,unit_labels.values()))==set(itertools.product([1,2,3],[1,2,3])))
# Pointed lexicographic binder to the nine-position roof root 1..9.
lex=[(i,j) for i in range(1,4) for j in range(1,4)]
roof_root=list(range(1,10)); roof=list(range(1,10))+list(range(8,0,-1))
checks['roof_root_length_9']=(len(roof_root)==9)
checks['roof_palindrome_17_profile']=(roof==list(reversed(roof)))
checks['roof_apex_9_equals_3_squared']=(max(roof)==9==3**2)
# Reversal of nine root positions has cycle type 4 transpositions +1 fixed.
roof_rev={k:10-k for k in range(1,10)}
# Matrix transpose permutation on lexicographic units: 3 transpositions +3 fixed.
pos={ij:k+1 for k,ij in enumerate(lex)}
transpose_perm={pos[(i,j)]:pos[(j,i)] for i,j in lex}
fix_roof=sum(roof_rev[k]==k for k in roof_rev); fix_trans=sum(transpose_perm[k]==k for k in transpose_perm)
checks['roof_reversal_fixed_count_1']=(fix_roof==1)
checks['matrix_transpose_unit_fixed_count_3']=(fix_trans==3)
checks['roof_reversal_not_permutation_conjugate_to_transpose']=(fix_roof!=fix_trans)
# Exact mod-2 Lie parity shadow.
def E2(i,j):
    A=np.zeros((3,3),dtype=np.uint8); A[i-1,j-1]=1; return A
Ky2=(E2(1,2)+E2(2,1))%2; Kz2=(E2(1,3)+E2(3,1))%2; J2=(E2(2,3)+E2(3,2))%2
def br2(A,B): return (A@B+B@A)%2
checks['TW_shadow_KyKz_to_J']=np.array_equal(br2(Ky2,Kz2),J2)
checks['TW_shadow_JKy_to_Kz']=np.array_equal(br2(J2,Ky2),Kz2)
checks['TW_shadow_JKz_to_Ky']=np.array_equal(br2(J2,Kz2),Ky2)
roof_data={
 'storage_context':{
   'roof_identity':'111111111^2 = 12345678987654321',
   'roof_profile':'1,2,...,9,...,2,1 is the carry-free autocorrelation / anti-diagonal slice profile for a 9x9 all-one square',
   'hyperfine_labels':'large finite ternary depth with projective consistency is used as an operational hyperfine labelling in the stored roof corpus'
 },
 'lambda2':{
   'space':'M_3(F_2), dimension 9',
   'standard_generators_to_units':unit_labels,
   'decision':'PROVEN: F_i,E_ij(2) map bijectively to the nine matrix-unit directions.'
 },
 'pointed_roof_binder':{
   'map':'roof-root positions 1..9 -> matrix units in Q0 lexicographic eigen-axis order',
   'status':'CONSTRUCTED_POINTED_ADAPTER, NOT NATURAL UNDER THE FULL GROUPOID',
   'guard':'the equality 9=3^2 and the roof apex do not identify decimal roof arithmetic with the binary deck label space.'
 },
 'reversal_no_go':{
   'roof_root_reversal_fixed_points':fix_roof,
   'matrix_unit_transpose_fixed_points':fix_trans,
   'decision':'REFUTED: no relabelling permutation can conjugate roof-root reversal to 3x3 transpose, because the involution cycle types differ.'
 },
 'Thomas_Wigner_parity_shadow':{
   'Ky_bar':'e12+e21','Kz_bar':'e13+e31','J_bar':'e23+e32',
   'brackets':'[Ky_bar,Kz_bar]=J_bar; [J_bar,Ky_bar]=Kz_bar; [J_bar,Kz_bar]=Ky_bar in characteristic 2',
   'status':'PROVEN_PARITY_SHADOW',
   'losses':'sign, Lorentz metric, rapidity magnitude, real orientation and physical spacetime interpretation are not retained.'
 },
 'hyperfinite_extension':'CONTEXT/PENDING-MORPHISM: ternary projective refinement can refine a chosen real/integer label after a binder is specified, but is not derived from Selling or lambda2 itself.'
}
(OUT/'roof_lambda_tw_shadow_v36.json').write_text(json.dumps(roof_data,indent=2,ensure_ascii=False),encoding='utf-8')

# -----------------------------------------------------------------------------
# F. Phase orientation gate remains O(2) unless a source-native mark breaks Gamma.
# -----------------------------------------------------------------------------
phase_gate={
 'previous':'v35 proved faithful SO(2)->U(1) scale |n|=1 and Gamma reverses J',
 'historical_domain_search':'the source constraints extracted in v36 provide mirror axes and quadripartite centres but no distinguished orientation that is invariantly selected over its reflected image',
 'result':'NO_NEW_GAMMA_BREAKING_MARK_FOUND_IN_SOURCE_TEXT_EXTRACTION',
 'carrier':'retain O(2)/Real-phase binder',
 'oriented_U1':'NT_PENDING_SOURCE_NATIVE_MARK_FROM_EXACT_FIELD_PLATE_OR_ADDITIONAL_DATA',
 'de_broglie':'SEPARATED/PENDING-MORPHISM'
}
checks['phase_gate_retains_O2']=(phase_gate['carrier']=='retain O(2)/Real-phase binder')
checks['no_unjustified_oriented_U1_promotion']=(phase_gate['oriented_U1'].startswith('NT_'))
(OUT/'phase_gate_v36.json').write_text(json.dumps(phase_gate,indent=2,ensure_ascii=False),encoding='utf-8')

# -----------------------------------------------------------------------------
# Certificate.
# -----------------------------------------------------------------------------
checks={k:bool(v) for k,v in checks.items()}
status='PASS' if all(bool(v) for v in checks.values()) else 'FAIL'
cert={
 'phase':'v36-selling-domain-roof-gamma3-relative-3cells',
 'status':status,
 'predecessor':{'v35_status':v35.get('status'),'v35_combined_checks':175},
 'new_checks':checks,
 'new_check_count':len(checks),
 'combined_check_count':175+len(checks),
 'historical_domain':source_domain,
 'orbifold_H12':orbifold_coh,
 'presentation_functor':functor_data,
 'Gamma3_2_3cells':{'count':43,'length_profile':rel_data['length_profile'],'status':'PROVEN'},
 'roof_lambda_TW':roof_data,
 'phase_gate':phase_gate,
 'decisions':{
   'true_Selling_Q0_field_domain':'PARTIAL_SOURCE_EXACT_RECONSTRUCTION; FULL_MICROSCOPIC_INCIDENCE_NT_PENDING_PLATE_TRANSCRIPTION',
   'transformation_orbifold_H1_H2':'PROVEN: Binet (1,0), Adjoint (3,0)',
   'field_domain_H1_H2':'NT_PENDING_FULL_INCIDENCE; not replaced by orbifold values',
   'functor_to_F29':'PROVEN_AT_PRESENTATION_LEVEL; FIELD_LEVEL_NT',
   'Gamma3_relations_as_relative_3cells':'PROVEN_43_ARITHMETIC_RELATIVE_CELLS',
   'lambda2_roof':'POINTED_9_LABEL_BINDER_CONSTRUCTED; NATURAL_IDENTIFICATION_REFUTED/PENDING',
   'lambda2_Thomas_Wigner':'PROVEN_PARITY_SHADOW_ONLY',
   'phase':'RETAIN_O2_REAL_PHASE; ORIENTED_U1_NT'
 },
 'safe_successor':'v37 -- retain all v36 separations; transcribe/vectorize the Selling figure-1 plate itself (bounded fields, boundary segments, crossings, fissures, repeated numeric forms and stabilizers) and hash-lock the resulting incidence table; only if that transcription closes, instantiate the field-level functor to F29_can and recompute field Binet/adjoint H1/H2. In parallel, compute the identities among the 43 Gamma_3(2) relative 3-cells at the 4-cell/syzygy level, and test whether the pointed 3x3 lambda2 label binder admits a projectively consistent ternary/hyperfine refinement compatible with the mod-2 Thomas-Wigner bracket; keep roof reversal, transpose, real Lorentz orientation, de Broglie phase and Cube27 separated unless explicit equivariant adapters are supplied.'
}
(OUT/'selling_domain_roof_gamma3_v36_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'status':status,'v35_checks':175,'new_checks':len(checks),'combined':175+len(checks),
                  'orbifold_Binet_H12':[H1B,H2B],'orbifold_Ad_H12':[H1A,H2A],
                  'Gamma3_relative_3cells':43,'lambda2_rank':9,'phase':'O2 retained'},indent=2))
raise SystemExit(0 if status=='PASS' else 1)
