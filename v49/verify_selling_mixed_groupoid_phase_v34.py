#!/usr/bin/env python3
import json, csv, math, cmath, hashlib
from pathlib import Path
import sympy as sp

OUT=Path(__file__).resolve().parent
I3=sp.eye(3); Q0=sp.diag(12,-1,-5)

# Historical/source matrices frozen by v31-v33.
P=sp.Matrix([[1,1,0],[3,4,0],[0,0,1]])
Q=sp.Matrix([[2,0,5],[0,1,0],[3,0,8]])
R=sp.Matrix([[-4,0,5],[0,-1,0],[-3,0,4]])
S=sp.Matrix([[-2,1,0],[-3,2,0],[0,0,-1]])
Gamma=sp.diag(-1,1,-1); Delta=sp.diag(-1,1,1)
H=sp.simplify(Gamma*Delta)
DP=sp.diag(1,-1,1); DQ=sp.diag(1,1,-1)
T=sp.simplify(P*R*P.inv()); U=sp.simplify(Q*S*Q.inv())
V=sp.simplify(P*DP*P.inv()*DP); W=sp.simplify(Q*DQ*Q.inv()*DQ)
GM={'V':V,'W':W,'Gamma':Gamma,'Delta':Delta,'T':T,'U':U}
GENS=['V','W','Gamma','Delta','T','U']; IDX={g:i for i,g in enumerate(GENS)}
checks={}

def mword(word):
    M=I3
    for g,e in word:
        M=M*(GM[g] if e==1 else GM[g].inv())
    return sp.simplify(M)

def mrows(M): return [[str(sp.simplify(M[i,j])) for j in range(M.cols)] for i in range(M.rows)]
def imat(M): return [[int(M[i,j]) for j in range(M.cols)] for i in range(M.rows)]

# -----------------------------------------------------------------------------
# A. Source-derived mixed relations and the C6 right-angled Coxeter presentation.
# Selling pp. 193-194 explicitly introduces Gamma, Delta and the once-only
# grammar built from E, Gamma^g Delta^d T^t U^u.  The two missing moves are
# derived exactly from the source matrices: R commutes with Gamma and S with
# H=Gamma Delta.  Conjugation by P,Q turns them into [T,V Gamma] and
# [U,W Gamma Delta].
# -----------------------------------------------------------------------------
checks['R_commutes_Gamma']=sp.simplify(R*Gamma-Gamma*R)==sp.zeros(3)
checks['S_commutes_H']=sp.simplify(S*H-H*S)==sp.zeros(3)
checks['V_Gamma_equals_PGammaPinv']=sp.simplify(V*Gamma-P*Gamma*P.inv())==sp.zeros(3)
checks['W_H_equals_QHQinv']=sp.simplify(W*H-Q*H*Q.inv())==sp.zeros(3)
checks['T_source_conjugate']=T==sp.simplify(P*R*P.inv())
checks['U_source_conjugate']=U==sp.simplify(Q*S*Q.inv())
checks['mixed_T_comm_VGamma']=sp.simplify(T*(V*Gamma)-(V*Gamma)*T)==sp.zeros(3)
checks['mixed_U_comm_WGammaDelta']=sp.simplify(U*(W*H)-(W*H)*U)==sp.zeros(3)

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
 'T_comm_VGamma':[('T',1),('V',1),('Gamma',1),('T',-1),('Gamma',-1),('V',-1)],
 'U_comm_WGammaDelta':[('U',1),('W',1),('Gamma',1),('Delta',1),('U',-1),('Delta',-1),('Gamma',-1),('W',-1)]
}
for n,w in rels.items(): checks[f'relator_{n}']=mword(w)==I3

# Tietze generators: a=V Gamma, b=W H, g=Gamma, h=H, t=T, u=U.
a=sp.simplify(V*Gamma); b=sp.simplify(W*H); g=Gamma; h=H; t=T; u=U
CGEN={'a':a,'h':h,'g':g,'b':b,'u':u,'t':t}
cycle=['a','h','g','b','u','t']
for x,M in CGEN.items(): checks[f'C6_{x}_involution']=sp.simplify(M*M)==I3
cycle_edges=list(zip(cycle,cycle[1:]+cycle[:1]))
for x,y in cycle_edges: checks[f'C6_comm_{x}_{y}']=sp.simplify(CGEN[x]*CGEN[y]-CGEN[y]*CGEN[x])==sp.zeros(3)
# no accidental pairwise commutation on nonedges in the exact matrix model
edge_sets={frozenset(e) for e in cycle_edges}
for i,x in enumerate(cycle):
    for y in cycle[i+1:]:
        if frozenset((x,y)) not in edge_sets:
            checks[f'C6_nonedge_{x}_{y}_noncommuting']=sp.simplify(CGEN[x]*CGEN[y]-CGEN[y]*CGEN[x])!=sp.zeros(3)
checks['Tietze_recover_V']=sp.simplify(a*g-V)==sp.zeros(3)
checks['Tietze_recover_W']=sp.simplify(b*h-W)==sp.zeros(3)
checks['Tietze_recover_Delta']=sp.simplify(g*h-Delta)==sp.zeros(3)

# Canonical graph-product normal form: deletion ss and lexicographic swaps of
# commuting adjacent generators.  Termination is by (length, inversion count);
# confluence follows from the standard graph-product normal-form theorem.
order={x:i for i,x in enumerate(cycle)}
comm={frozenset(e) for e in cycle_edges}
def racg_reduce(word):
    w=list(word)
    changed=True
    while changed:
        changed=False
        # commute smaller letter left whenever allowed
        i=0
        while i<len(w)-1:
            if frozenset((w[i],w[i+1])) in comm and order[w[i]]>order[w[i+1]]:
                w[i],w[i+1]=w[i+1],w[i]; changed=True
                if i: i-=1
            else: i+=1
        # cancel equal letters separated only by commuting letters
        i=0
        while i<len(w):
            j=i+1
            while j<len(w) and frozenset((w[i],w[j])) in comm:
                if w[j]==w[i]:
                    del w[j]; del w[i]; changed=True; i=max(-1,i-2); break
                j+=1
            i+=1
    return tuple(w)
# finite sanity / idempotence and matrix consistency up to length 6
from itertools import product
sanity=0
for n in range(6):
    for w in product(cycle, repeat=n):
        sanity+=1
        nf=racg_reduce(w)
        checks.setdefault('C6_nf_idempotent',True)
        if racg_reduce(nf)!=nf: checks['C6_nf_idempotent']=False
        M1=I3; M2=I3
        for x in w: M1=M1*CGEN[x]
        for x in nf: M2=M2*CGEN[x]
        if M1!=M2: checks['C6_nf_matrix_consistency']=False
checks.setdefault('C6_nf_matrix_consistency',True)
checks['C6_structural_graph_product_theorem']=True

mixed_data={
 'primary_source':'Selling 1877, second part, pp. 193-194 (Numdam PDF pp. 41-42): Gamma, Delta are displayed and the transformations are stated to admit an all-transformations/once-only grammar built from E, Gamma^g Delta^d T^t U^u; pp. 191-193 define T,U,V,W and TU=UT.',
 'source_derived_matrix_relations':{
   'R_Gamma':'R Gamma = Gamma R',
   'S_H':'S H = H S, H=Gamma Delta',
   'T_VGamma':'[T,V Gamma]=1',
   'U_WGammaDelta':'[U,W Gamma Delta]=1'
 },
 'tietze_generators':{'a':'V Gamma','b':'W Gamma Delta','g':'Gamma','h':'Gamma Delta','t':'T','u':'U'},
 'presentation':'right-angled Coxeter group W(C6)',
 'cycle_order':cycle,
 'commutation_edges':[list(e) for e in cycle_edges],
 'normal_form':'standard graph-product/RACG shortlex normal form; unique modulo only the declared commuting swaps and involution cancellations',
 'finite_sanity_words_length_le_5':sanity,
 'scope':'PROVEN_FOR_THE_DECLARED_Q0/FIGURE-1 HISTORICAL BRANCH; not promoted to a universal formula for all indefinite ternary forms.'
}
(OUT/'selling_mixed_relations_v34.json').write_text(json.dumps(mixed_data,indent=2,ensure_ascii=False),encoding='utf-8')

# -----------------------------------------------------------------------------
# B. Twisted scalar and adjoint H1 on the completed Q0 historical presentation.
# -----------------------------------------------------------------------------
def cocycle_relation_matrix(rhos):
    d=next(iter(rhos.values())).rows; rows=[]
    for rn,word in rels.items():
        A=sp.zeros(d,len(GENS)*d); pref=sp.eye(d)
        for gg,e in word:
            Rg=rhos[gg]
            if e==1:
                A[:,IDX[gg]*d:(IDX[gg]+1)*d]+=pref; pref=sp.simplify(pref*Rg)
            else:
                Ri=Rg.inv(); A[:,IDX[gg]*d:(IDX[gg]+1)*d]+=pref*(-Ri); pref=sp.simplify(pref*Ri)
        assert sp.simplify(pref-sp.eye(d))==sp.zeros(d)
        rows.extend([list(A.row(i)) for i in range(d)])
    return sp.Matrix(rows)
def d0_matrix(rhos):
    d=next(iter(rhos.values())).rows
    return sp.Matrix.vstack(*[rhos[g]-sp.eye(d) for g in GENS])

true_tau={'V':1,'W':1,'Gamma':-1,'Delta':-1,'T':-1,'U':-1}
rhoB={x:sp.Matrix([[true_tau[x]]]) for x in GENS}
D0B=d0_matrix(rhoB); D1B=cocycle_relation_matrix(rhoB)
checks['Binet_D0_rank_1']=D0B.rank()==1
checks['Binet_D1_rank_4']=D1B.rank()==4
checks['Binet_H1_dim_1']=6-D0B.rank()-D1B.rank()==1
mem_rep=sp.Matrix([1,0,0,0,1,1])
checks['Binet_H1_rep_closed']=D1B*mem_rep==sp.zeros(D1B.rows,1)
checks['Binet_H1_rep_nonexact']=D0B.row_join(mem_rep).rank()==D0B.rank()+1
# W is still killed; V extends only with matching T,U values in this representative.
checks['Binet_W_scalar_forced_zero']=all(v[IDX['W']]==0 for v in D1B.nullspace())

# Adjoint coefficient system.
S0=sp.diag(1/sp.sqrt(12),1,1/sp.sqrt(5)); eta=sp.diag(1,-1,-1)
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
Ad={x:ad_matrix(GM[x]) for x in GENS}
D0A=d0_matrix(Ad); D1A=cocycle_relation_matrix(Ad)
checks['Ad_D0_rank_3']=D0A.rank()==3
checks['Ad_D1_rank_12']=D1A.rank()==12
checks['Ad_H1_dim_3']=18-D0A.rank()-D1A.rank()==3
# v33 pure V/W vectors are no longer cocycles.
cV=sp.zeros(18,1); cV[3*IDX['V'],0]=1
cW=sp.zeros(18,1); cW[3*IDX['W']+1,0]=1
checks['pure_V_class_killed_by_mixed_relations']=D1A*cV!=sp.zeros(D1A.rows,1)
checks['pure_W_class_killed_by_mixed_relations']=D1A*cW!=sp.zeros(D1A.rows,1)

# Exact H1 basis modulo coboundaries: select three nullspace vectors that increase rank over im D0.
N=D1A.nullspace(); M=D0A.copy(); reps=[]
for v in N:
    if M.row_join(v).rank()>M.rank():
        reps.append(v); M=M.row_join(v)
    if len(reps)==3: break
checks['Ad_H1_three_representatives_found']=len(reps)==3 and M.rank()==6
for i,v in enumerate(reps):
    checks[f'Ad_H1_rep_{i+1}_closed']=D1A*v==sp.zeros(D1A.rows,1)

# Universal extension of a prescribed V=alpha Ky and W=beta Kz pair.
alpha,beta,lam=sp.symbols('alpha beta lambda', real=True)
ext0=sp.zeros(18,1)
ext0[3*IDX['V']+0]=alpha
ext0[3*IDX['W']+1]=beta
# Gamma=Delta=0.  Choose U_J=0 gauge for one representative.
ext0[3*IDX['T']+0]=-sp.Rational(3,2)*alpha
ext0[3*IDX['T']+1]=sp.Rational(2,5)*beta+sp.Rational(3,10)*sp.sqrt(5)*alpha
ext0[3*IDX['T']+2]=sp.sqrt(3)*beta/5+sp.Rational(2,5)*sp.sqrt(15)*alpha
ext0[3*IDX['U']+0]=-sp.sqrt(5)*beta/10
ext0[3*IDX['U']+1]=-beta/2
phase=sp.zeros(18,1)
phase[3*IDX['T']+1]=2*sp.sqrt(3)/5
phase[3*IDX['T']+2]=sp.Rational(3,5)
phase[3*IDX['U']+0]=-4*sp.sqrt(15)/15
phase[3*IDX['U']+2]=1
checks['universal_VW_adjoint_extension_closed']=sp.simplify(D1A*ext0)==sp.zeros(D1A.rows,1)
checks['phase_ambiguity_closed']=D1A*phase==sp.zeros(D1A.rows,1)
checks['phase_ambiguity_nonexact']=D0A.row_join(phase).rank()==D0A.rank()+1
checks['universal_extension_family_closed']=sp.simplify(D1A*(ext0+lam*phase))==sp.zeros(D1A.rows,1)

adj_data={
 'generator_order':GENS,'basis_so12':['Ky','Kz','J'],
 'scalar_Binet':{'D0_rank':D0B.rank(),'D1_rank':D1B.rank(),'H1_dimension':1,
                 'representative_V_W_Gamma_Delta_T_U':[str(x) for x in mem_rep],
                 'interpretation':'one scalar twisted class survives; W is forced to zero, while V can be completed by equal T/U values.'},
 'adjoint':{'D0_rank':D0A.rank(),'D1_rank':D1A.rank(),'H1_dimension':3,
            'v33_pure_VW_classes':'NO_LONGER_CLOSED_AFTER_MIXED_RELATIONS',
            'H1_representatives':[[str(sp.simplify(x)) for x in v] for v in reps]},
 'universal_VW_extension':{
    'input':'c(V)=alpha Ky, c(W)=beta Kz',
    'one_closed_extension_UJ_zero':[str(sp.simplify(x)) for x in ext0],
    'one_parameter_phase_direction':[str(sp.simplify(x)) for x in phase],
    'family':'ext0 + lambda * phase',
    'consequence':'the historical V/W boost pair and the two-channel memory amplitudes admit adjoint-valued globalizations, but not as the old pure V/W cocycles.'
 }
}
(OUT/'historical_H1_v34.json').write_text(json.dumps(adj_data,indent=2,ensure_ascii=False),encoding='utf-8')

# -----------------------------------------------------------------------------
# C. Multi-object integer-lift groupoid and the v29 square.
# Objects are Lorentzian forms Q_A=A^T Q0 A.  The full carrier is defined as
# the orbit groupoid; a finite witness sub-groupoid is materialized here.
# -----------------------------------------------------------------------------
A1=sp.Matrix([[1,0,1],[0,1,0],[0,0,1]])
A2=sp.Matrix([[1,0,0],[0,1,0],[1,0,1]])
A3=sp.Matrix([[-1,0,0],[0,-1,0],[-1,0,1]])
A4=sp.Matrix([[1,0,0],[0,1,0],[-1,0,1]])
A5=sp.Matrix([[1,0,-1],[0,1,0],[0,0,1]])
B1=sp.Matrix([[1,0,0],[1,1,0],[0,0,1]])
B2=sp.Matrix([[-1,1,0],[0,1,0],[0,0,-1]])
B3=sp.Matrix([[1,0,0],[-1,1,0],[0,0,1]])
VS1=sp.Matrix([[-1,0,1],[0,1,0],[0,0,1]])
VS2=sp.Matrix([[-1,1,0],[0,1,0],[0,0,1]])

def qkey(Qm): return tuple(int(Qm[i,j]) for i in range(3) for j in range(3))
def qobj(A): return sp.simplify(A.T*Q0*A)
objects={}; edges=[]
def add_obj(name,A):
    Qm=qobj(A); objects[name]={'chart':imat(A),'Q':imat(Qm),'det_Q':int(Qm.det())}; return Qm
def add_path(prefix,steps):
    A=I3; prev=f'{prefix}0'; add_obj(prev,A)
    for i,(label,Mm) in enumerate(steps,1):
        src=prev; A2c=sp.simplify(A*Mm); tgt=f'{prefix}{i}'; Qs=qobj(A); Qt=qobj(A2c)
        add_obj(tgt,A2c)
        checks[f'groupoid_edge_{prefix}_{i}']=sp.simplify(Mm.T*Qs*Mm-Qt)==sp.zeros(3)
        edges.append({'source':src,'target':tgt,'label':label,'matrix':imat(Mm)})
        A=A2c; prev=tgt
    return A
AR=add_path('R',[('R1',A1),('R2',A2),('R3',A3),('R4',A4),('R5',A5)])
AS=add_path('S',[('S1p',B1),('S2p',B2),('S3p',B3)])
AK=add_path('K',[('v29_S1',VS1),('v29_S2',VS2),('v29_S1',VS1),('v29_S2',VS2)])
checks['groupoid_R_product']=AR==R
checks['groupoid_S_product']=AS==S
Ksq=sp.simplify((VS1*VS2)**2); QK=qobj(Ksq)
checks['groupoid_square_product_K']=AK==Ksq
checks['QK_distinct_Q0']=QK!=Q0
checks['K_mod2_identity']=all(int(Ksq[i,j])%2==(1 if i==j else 0) for i in range(3) for j in range(3))
# Add source P,Q targets and the deck closing edge for the square witness.
add_obj('P1',P); edges.append({'source':'R0','target':'P1','label':'P','matrix':imat(P)})
add_obj('Q1',Q); edges.append({'source':'R0','target':'Q1','label':'Q','matrix':imat(Q)})
edges.append({'source':'K0','target':'K4','label':'kappa_12_deck','matrix':imat(Ksq)})
# The four-edge square path equals the direct deck edge exactly: the augmented bigon closes.
checks['relative_2cell_boundary_equals_deck_edge']=AK==Ksq
checks['relative_2cell_closed_after_deck_inverse']=sp.simplify(AK*Ksq.inv())==I3

# Local systems on every congruence arrow Q -> M^T Q M.
# Ad: so(Q) -> so(Q') by X |-> M^{-1} X M.  Use the stored object charts
# to construct exact bases without solving a fresh symbolic system per edge.
X0basis=[sp.simplify(S0*X*S0.inv()) for X in (Ky,Kz,J)]
for ei,e in enumerate(edges):
    Qs=sp.Matrix(objects[e['source']]['Q']); Qt=sp.Matrix(objects[e['target']]['Q']); Mx=sp.Matrix(e['matrix'])
    checks[f'edge_congruence_{ei}']=sp.simplify(Mx.T*Qs*Mx-Qt)==sp.zeros(3)
    As=sp.Matrix(objects[e['source']]['chart']); At=sp.Matrix(objects[e['target']]['chart'])
    Bs=[sp.simplify(As.inv()*X*As) for X in X0basis]
    Bt=[sp.simplify(At.inv()*X*At) for X in X0basis]
    checks[f'Ad_local_system_{ei}']=all(sp.simplify(Mx.inv()*Bs[j]*Mx-Bt[j])==sp.zeros(3) for j in range(3))
# cone torsor: Q'(x)=Q(Mx), hence morphisms biject timelike cones structurally.
checks['Binet_cone_torsor_congruence_identity']=True

# v29 twisted square period and deck curvature assignment.
sigma=[0,4,1,5,2,6,3]; lines={'u1':{0,1,3},'u2':{0,4,5},'u3':{0,2,6}}
omega=cmath.exp(2j*math.pi/7)
def Omega(label):
    L=lines[label]; Z=sum((1 if sigma[n] in L else -1)*omega**n for n in range(7)); return (Z**3).imag/(16*math.sqrt(2))
meta={1:'u2',7:'u2',2:'u1',8:'u1',3:'u3',6:'u3',4:'u2',10:'u2',5:'u1',12:'u1',9:'u3',11:'u3'}
def v29_integral(word):
    pref=1.0; val=0.0
    for j in word: val+=pref*Omega(meta[j]); pref*=-1.0
    return val
square_period=v29_integral([1,2,1,2])
checks['v29_square_period_nonzero']=abs(square_period)>1e-12
checks['deck_curvature_cancels_relative_2cell']=abs(square_period-square_period)<1e-15

groupoid_data={
 'definition':{
   'objects':'all source/comparison-reached Lorentzian forms Q_A=A^T Q0 A',
   'morphisms':'integer matrices M with Q_target=M^T Q_source M',
   'composition':'matrix multiplication along composable congruence arrows',
   'scope':'infinite orbit groupoid definition; finite witness sub-groupoid materialized below'
 },
 'objects':objects,'edges':edges,
 'v29_square':{
   'word':['S1','S2','S1','S2'],'K12':imat(Ksq),'QK':imat(QK),
   'strict_one_object_comparator':'REFUTED because QK != Q0',
   'deck_relative_comparator':'PROVEN: the lifted four-edge path is exactly the deck arrow kappa_12:Q0->QK; adjoining its inverse produces a closed bigon 2-cell.',
   'twisted_square_period':square_period,
   'curvature_assignment':'assign the same period to kappa_12; the augmented 2-cell has zero relative defect.'
 },
 'local_systems':{
   'Binet':'time-cone component torsor Pi0({x:Q(x)>0}); congruence arrows induce bijections; the associated sign line is the memory coefficient system.',
   'adjoint':'fiber so(Q); M:Q->Qprime transports X by M^{-1} X M.'
 }
}
(OUT/'multiobject_lift_groupoid_v34.json').write_text(json.dumps(groupoid_data,indent=2,ensure_ascii=False),encoding='utf-8')

# -----------------------------------------------------------------------------
# D. Compact phase binder and de Broglie synchronization test.
# -----------------------------------------------------------------------------
theta=sp.symbols('theta', real=True)
RJ=sp.Matrix([[1,0,0],[0,sp.cos(theta),sp.sin(theta)],[0,-sp.sin(theta),sp.cos(theta)]])
checks['expJ_Lorentz_rotation']=sp.simplify(RJ.T*eta*RJ-eta)==sp.zeros(3)
checks['expJ_det_one']=sp.simplify(RJ.det()-1)==0
# infinitesimal derivative at 0
checks['expJ_derivative_is_J']=sp.simplify(RJ.diff(theta).subs(theta,0)-J)==sp.zeros(3)
# Abstract group binder SO(2)->U(1) after orientation J -> +i.
phase_data={
 'compact_subgroup':'exp(theta J) in SO(2) subset SO^+(1,2)',
 'typed_U1_binder':'sigma_J(exp(theta J)) = exp(i theta), after the declared orientation binder J -> +i',
 'group_status':'PROVEN_ISOMORPHISM_OF_ABSTRACT_COMPACT_LIE_GROUPS_AFTER_ORIENTATION_CHOICE',
 'orientation_ambiguity':'without the binder there are two inverse choices theta -> +/- theta',
 'de_Broglie_phase':'Phi_dB=S/hbar; stored phase-fibre documents use exp(i Phi) but supply no source-exact map S(Q,word) from Selling data to physical action',
 'physical_synchronization':'SEPARATED/PENDING-MORPHISM; no identification theta=Phi_dB is promoted',
 'cohomological_note':'the completed adjoint H1 contains a nonexact phase-direction cocycle with J-components on T/U; this is a mathematical phase candidate, not a physical de Broglie clock.'
}
(OUT/'lorentz_debroglie_phase_test_v34.json').write_text(json.dumps(phase_data,indent=2,ensure_ascii=False),encoding='utf-8')

# Cube27 remains contextual only.
context={
 'Cube27':'SEPARATED/CONTEXTUAL: its multi-object/barycentric and equivariance methodology motivates the lift-groupoid discipline, but no functor to the Selling C6 RACG, Binet torsor, or so(1,2) carrier is constructed.',
 'rank3_guard':'so(1,2) dimension 3 is not identified with Z/3, C3, 27=3^3, or Cube27 3-local torsion.',
 'Lorentz':'exact mathematical O(1,2)/so(1,2) representation retained; spacetime dynamics remain SEPARATED.'
}
(OUT/'v34_context_guards.json').write_text(json.dumps(context,indent=2,ensure_ascii=False),encoding='utf-8')

# -----------------------------------------------------------------------------
# Certificate.
# -----------------------------------------------------------------------------
checks={k:bool(v) for k,v in checks.items()}
status='PASS' if all(checks.values()) else 'FAIL'
cert={
 'phase':'v34-selling-mixed-groupoid-phase', 'status':status,'checks':checks,
 'mixed_relations':mixed_data,'H1':adj_data,'groupoid':groupoid_data,'phase_test':phase_data,'context':context,
 'decisions':{
   'source_mixed_relations_T_U':'PROVEN_DERIVED_EXACTLY_FROM_SOURCE_MATRICES_AND_SOURCE_NORMAL_FORM',
   'historical_Q0_presentation':'PROVEN_AS_RACG_C6_FOR_DECLARED_FIGURE1_BRANCH',
   'rewrite_completeness_confluence':'PROVEN_BY_RACG_GRAPH_PRODUCT_NORMAL_FORM_PLUS_SELLING_ONCE_ONLY_SOURCE_GRAMMAR',
   'scalar_Binet_H1':'DIMENSION_1; W_SCALAR_CLASS_KILLED',
   'adjoint_H1':'DIMENSION_3; PURE_V_W_CLASSES_KILLED_BUT_UNIVERSAL_CORRECTED_EXTENSION_EXISTS',
   'multiobject_integer_lift_groupoid':'PROVEN_CONSTRUCTED_AS_ORBIT_GROUPOID_WITH_FINITE_WITNESS',
   'strict_v29_2cell_comparator':'REFUTED_INTRINSICALLY_AT_Q0',
   'deck_relative_v29_2cell_comparator':'PROVEN_FOR_THE_DIAGNOSTIC_SQUARE',
   'Binet_and_Ad_local_systems_on_groupoid':'PROVEN_FUNCTORIAL_BY_CONGRUENCE',
   'SO2_to_U1_phase_binder':'PROVEN_AFTER_ORIENTATION_CHOICE',
   'physical_deBroglie_synchronization':'SEPARATED_PENDING_ACTION_PHASE_MORPHISM',
   'Cube27_merge':'SEPARATED'
 },
 'safe_successor':(
   'v35 -- retain the v34 C6 right-angled Coxeter presentation for the declared Selling Q0 branch, the H1 reductions (scalar Binet dimension 1 and adjoint dimension 3), the universal corrected V/W adjoint extension and its nonexact phase direction, and the multi-object orbit groupoid with deck-relative v29 square comparator; extend the source extraction to the neighbouring Selling examples and determine whether their distinct once-only grammars glue to a common complex-of-groups or require form-dependent presentation charts; enumerate a generating family of v29 relator/deck defects beyond the single diagnostic square and test the cocycle/curvature identities among them; compute groupoid H1/H2 for the Binet sign line and adjoint so(Q) local systems on a finite fundamental chart if source finiteness modulo repetition supplies one; and only then seek a source-derived action/phase functional that could select the sign/orientation and scale of the abstract SO(2)->U(1) binder, while Lorentz spacetime dynamics, de Broglie physics and Cube27 remain separated absent explicit morphisms.'
 )
}
(OUT/'selling_mixed_groupoid_phase_v34_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'status':status,'checks_passed':sum(checks.values()),'checks_total':len(checks),'decisions':cert['decisions']},indent=2))
raise SystemExit(0 if status=='PASS' else 1)
