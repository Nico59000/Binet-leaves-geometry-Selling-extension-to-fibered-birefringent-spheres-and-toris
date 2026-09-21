#!/usr/bin/env python3
"""Source-typed Selling field-construction primitives for v50.

This is deliberately fail-closed.  It implements only algebraic steps already
source-anchored in the v24--v31 lineage: homogeneous coefficients/conorms,
Selling reduction inequalities, the twelve adjacent superbase moves, exact
congruence transport, and a finite-cell output schema.  It does NOT infer a
historical figure from raster morphology or invent missing figure seeds.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from fractions import Fraction
from itertools import permutations
from pathlib import Path
import argparse, json
import sympy as sp

# ---------------------------------------------------------------------------
# Exact historical/adjacent matrices already certified in the lineage.
# ---------------------------------------------------------------------------
S = [
sp.Matrix([[-1,0,1],[0,1,0],[0,0,1]]),
sp.Matrix([[-1,1,0],[0,1,0],[0,0,1]]),
sp.Matrix([[-1,1,1],[0,1,0],[0,0,1]]),
sp.Matrix([[0,-1,0],[-1,0,0],[-1,-1,1]]),
sp.Matrix([[0,0,-1],[-1,1,-1],[-1,0,0]]),
sp.Matrix([[1,-1,-1],[0,0,-1],[0,-1,0]]),
sp.Matrix([[1,0,0],[0,-1,1],[0,0,1]]),
sp.Matrix([[1,0,0],[0,1,0],[0,1,-1]]),
sp.Matrix([[1,0,0],[0,1,0],[1,0,-1]]),
sp.Matrix([[1,0,0],[0,1,0],[1,1,-1]]),
sp.Matrix([[1,0,0],[1,-1,0],[0,0,1]]),
sp.Matrix([[1,0,0],[1,-1,1],[0,0,1]])]

Q0 = sp.diag(12,-1,-5)
P = sp.Matrix([[1,1,0],[3,4,0],[0,0,1]])
QSUB = sp.Matrix([[2,0,5],[0,1,0],[3,0,8]])
R = sp.Matrix([[-4,0,5],[0,-1,0],[-3,0,4]])
SMAC = sp.Matrix([[-2,1,0],[-3,2,0],[0,0,-1]])
T = P*R*P.inv()
U = QSUB*SMAC*QSUB.inv()

WALLS = ('p01','p02','p03','p12','p13','p23')
WALL_TO_MOVES = {
    'p01': (3,6),
    'p02': (5,12),
    'p03': (4,10),
    'p12': (1,7),
    'p13': (2,8),
    'p23': (9,11),
}


def _int_matrix(M):
    return [[int(sp.Integer(M[i,j])) for j in range(M.cols)] for i in range(M.rows)]

def gram_from_six(a,b,c,g,h,k):
    return sp.Matrix([[a,k,h],[k,b,g],[h,g,c]])

def six_from_gram(M):
    return (sp.expand(M[0,0]),sp.expand(M[1,1]),sp.expand(M[2,2]),
            sp.expand(M[1,2]),sp.expand(M[0,2]),sp.expand(M[0,1]))

def conorms(M):
    a,b,c,g,h,k=six_from_gram(M)
    vals=(a+h+k,b+k+g,c+h+g,-k,-h,-g)
    return dict(zip(WALLS,map(sp.expand,vals)))

def homogeneous(M):
    """Selling homogeneous six (g,h,k,l,m,n) for a positive ternary form."""
    a,b,c,g,h,k=six_from_gram(M)
    l=sp.expand(-a-h-k)
    m=sp.expand(-b-k-g)
    n=sp.expand(-c-g-h)
    return {'g':g,'h':h,'k':k,'l':l,'m':m,'n':n}

def is_reduced_positive(M):
    """Exact when entries are ordered rationals/integers; otherwise returns symbolic predicates."""
    hh=homogeneous(M)
    return all(bool(sp.Le(v,0)) for v in hh.values())

def positive_correspondent(Qind, xi, eta, zeta):
    """Selling relation Q+ = 2 vv^T - Q_ind in the chosen marked chart."""
    v=sp.Matrix([xi,eta,zeta])
    return sp.expand(2*(v*v.T)-Qind)

def active_walls(M):
    return [k for k,v in conorms(M).items() if sp.simplify(v)==0]

def wall_signs(M):
    out={}
    for k,v in conorms(M).items():
        vv=sp.simplify(v)
        out[k]=0 if vv==0 else (1 if bool(sp.Gt(vv,0)) else -1 if bool(sp.Lt(vv,0)) else None)
    return out

def adjacent_candidates(M, wall):
    if wall not in WALL_TO_MOVES: raise KeyError(wall)
    ans=[]
    for idx in WALL_TO_MOVES[wall]:
        A=S[idx-1]
        ans.append({'move':f'S{idx}','matrix':_int_matrix(A),'target':_int_matrix(A.T*M*A),
                    'target_conorms':{k:str(v) for k,v in conorms(A.T*M*A).items()}})
    return ans

# Generate S4 relabellings of a superbase (v0=-v1-v2-v3).
def relabelling_matrices():
    V=[sp.Matrix([-1,-1,-1]),sp.eye(3)[:,0],sp.eye(3)[:,1],sp.eye(3)[:,2]]
    mats=[]
    for pi in permutations(range(4)):
        A=sp.Matrix.hstack(V[pi[1]],V[pi[2]],V[pi[3]])
        if abs(int(A.det()))==1:
            key=tuple(map(int,list(A)))
            if key not in {tuple(map(int,list(B))) for B in mats}: mats.append(A)
    return mats
RELABEL=relabelling_matrices()

def numerical_form_key(M):
    """Canonical exact key under the 24 source superbase relabellings."""
    keys=[]
    for A in RELABEL:
        N=A.T*M*A
        keys.append(tuple(int(x) if x.q==1 else str(x) for x in list(N)))
    return min(keys,key=lambda x:tuple(map(str,x)))

@dataclass(frozen=True)
class FieldRecord:
    field_id:str
    gram:list
    numerical_key:list
    repetition_of:str|None

@dataclass(frozen=True)
class EdgeRecord:
    edge_id:str
    source:str
    target:str
    wall:str
    move:str
    matrix:list

@dataclass(frozen=True)
class VertexRecord:
    vertex_id:str
    kind:str
    vanishing:list
    incident_edges:list


def blank_historical_complex_schema():
    return {
      'C0_vertices':[],
      'C1_oriented_boundary_arcs':[],
      'C2_fields':[],
      'ordered_boundary_cycles':{},
      'repetition_identifications':{},
      'stabilizer_words':{},
      'Binet_sign_transports':{},
      'Ad_transports':{},
      'guards':[
        'A field/edge/vertex is historical only after source construction or source-pointed validation.',
        'Raster morphology alone never promotes a historical cell.',
        'H1/H2 are forbidden until every C2 boundary is ordered and d1*d0=0 is exact.'
      ]}


def source_known_fig1_seed():
    return {
      'figure':'Fig.1','invariant':60,
      'Q0':_int_matrix(Q0),
      'P':_int_matrix(P),'Q':_int_matrix(QSUB),'R':_int_matrix(R),'S':_int_matrix(SMAC),
      'T':_int_matrix(T),'U':_int_matrix(U),
      'checks':{
        'det_Q0':int(Q0.det()),
        'det_P':int(P.det()),'det_Q':int(QSUB.det()),'det_R':int(R.det()),'det_S':int(SMAC.det()),
        'T_stabilizes_Q0':T.T*Q0*T==Q0,
        'U_stabilizes_Q0':U.T*Q0*U==Q0,
        'TU_commute':T*U==U*T,
      }}



# ---------------------------------------------------------------------------
# v50: exact source-feasibility layer from Selling, 1877, pp. 160, 166--170.
# ---------------------------------------------------------------------------
def superbase_coefficients(M):
    """Return Selling's a,b,c,d,g,h,k,l,m,n for v0=-v1-v2-v3.

    The 3x3 Gram uses (a,b,c,g,h,k); the fourth diagonal and three
    missing pairings are exact consequences of the superbase relation.
    """
    a,b,c,g,h,k=six_from_gram(M)
    d=sp.expand(a+b+c+2*(g+h+k))
    l=sp.expand(-a-h-k)
    m=sp.expand(-b-k-g)
    n=sp.expand(-c-g-h)
    return {'a':a,'b':b,'c':c,'d':d,'g':g,'h':h,'k':k,'l':l,'m':m,'n':n}

def adjoint_coefficients(M):
    """Diagonal/coefficient data of the adjoint ternary form."""
    a,b,c,g,h,k=six_from_gram(M)
    return {
      'A':sp.expand(b*c-g*g),
      'B':sp.expand(c*a-h*h),
      'C':sp.expand(a*b-k*k),
      'G':sp.expand(h*k-a*g),
      'H':sp.expand(k*g-b*h),
      'K':sp.expand(g*h-c*k),
      'I':sp.expand(M.det()),
    }

def binary_reduced_feasible(u,v,w):
    """Selling's existence criterion for an indefinite binary reduced form.

    Part I, p.31: the two extreme coefficients have opposite signs.
    The optional first-positive convention is only a representative choice,
    not an existence condition, so it is deliberately not imposed here.
    """
    return bool((sp.Integer(u)>0 and sp.Integer(w)<0) or
                (sp.Integer(u)<0 and sp.Integer(w)>0))

def _same_sign_product_gate(x,y,I):
    """Selling p.166 sign/product gate used for candidate crossings."""
    x=sp.Integer(x); y=sp.Integer(y); I=sp.Integer(I)
    same=(x>0 and y>0) or (x<0 and y<0)
    if not same: return False
    return bool(x*y<=I) if x<0 and y<0 else True

def standard_hk_crossing_feasible(M):
    """Exact feasibility of a crossing h'=k'=0 for the positive correspondent.

    Selling (19): a and A=bc-g^2 must have the same sign; if both are
    negative, a*A <= I; and the binary form (b,g,c) must be reduced.
    """
    a,b,c,g,h,k=six_from_gram(M); A=sp.expand(b*c-g*g); I=sp.expand(M.det())
    return _same_sign_product_gate(a,A,I) and binary_reduced_feasible(b,g,c)

def hk_same_field_candidate_tests(M):
    """Three source-exact candidate crossings from h'=k'=0 (Selling p.166).

    Candidate labels are the positive-correspondent homogeneous lines:
      h'=l'=0, h'=g'=0, h'=n'=0.
    The arithmetic gates are respectively
      a ~ cd-n^2 with binary (c,n,d),
      c ~ ab-k^2 with binary (a,k,b),
      c ~ ad-l^2 with binary (a,l,d).
    Here '~' denotes the source same-sign/product gate.
    """
    q=superbase_coefficients(M); I=sp.expand(M.det())
    a,b,c,d,g,h,k,l,m,n=(q[x] for x in ('a','b','c','d','g','h','k','l','m','n'))
    data={
      'h_l': {'pair':['h','l'],'lhs':a,'rhs':sp.expand(c*d-n*n),'binary':[c,n,d]},
      'h_g': {'pair':['h','g'],'lhs':c,'rhs':sp.expand(a*b-k*k),'binary':[a,k,b]},
      'h_n': {'pair':['h','n'],'lhs':c,'rhs':sp.expand(a*d-l*l),'binary':[a,l,d]},
    }
    out={}
    for name,r in data.items():
        sg=_same_sign_product_gate(r['lhs'],r['rhs'],I)
        br=binary_reduced_feasible(*r['binary'])
        out[name]={
          'pair':r['pair'],'lhs':int(r['lhs']),'rhs':int(r['rhs']),
          'same_sign_product_gate':bool(sg),'binary_reduced':bool(br),
          'feasible':bool(sg and br),'binary':[int(x) for x in r['binary']]
        }
    return out

def q0_anisotropy_mod5_certificate():
    """Infinite-descent certificate: 12x^2-y^2-5z^2=0 has no Q-point."""
    # Mod 5: y^2 = 2 x^2. Since 2 is a nonsquare mod 5, 5|x,y;
    # substituting then forces 5|z, contradicting primitiveness.
    squares={0,1,4}
    return {
      'modulus':5,'squares_mod5':sorted(squares),'two_is_nonsquare':2 not in squares,
      'step1':'y^2 = 2 x^2 (mod 5) => 5|x and 5|y',
      'step2':'substitution in 12x^2-y^2-5z^2=0 => 5|z',
      'conclusion':'NO_NONZERO_RATIONAL_ISOTROPIC_VECTOR'
    }

def q0_hyperboloid_equation():
    xi,eta,zeta=sp.symbols('xi eta zeta')
    v=sp.Matrix([xi,eta,zeta])
    eq=sp.expand(60*((v.T*Q0.inv()*v)[0]-1))
    return sp.expand(eq) # 5 xi^2 - 60 eta^2 - 12 zeta^2 - 60

def v50_source_feasibility_report():
    # Deterministic relabelling placing Q0 in an hk-crossing-admissible marking.
    Qstart=sp.diag(-1,12,-5)
    cand=hk_same_field_candidate_tests(Qstart)
    cert=q0_anisotropy_mod5_certificate()
    checks={
      'q0_hyperboloid_exact': q0_hyperboloid_equation()==sp.Symbol('xi')**2*5-sp.Symbol('eta')**2*60-sp.Symbol('zeta')**2*12-60,
      'qstart_det60': Qstart.det()==60,
      'qstart_hk_crossing': standard_hk_crossing_feasible(Qstart),
      'qstart_three_candidate_records': set(cand)=={'h_l','h_g','h_n'},
      'qstart_candidates_source_gates': all(v['feasible'] for v in cand.values()),
      'q0_anisotropic_mod5': cert['two_is_nonsquare'] and cert['conclusion']=='NO_NONZERO_RATIONAL_ISOTROPIC_VECTOR',
    }
    return {
      'version':'v50','status':'PASS' if all(checks.values()) else 'FAIL',
      'checks':checks,'Q0':_int_matrix(Q0),'Qstart':_int_matrix(Qstart),
      'hyperboloid':'5 xi^2 - 60 eta^2 - 12 zeta^2 = 60',
      'anisotropy_certificate':cert,'hk_candidate_tests':cand,
      'source_branch_rules':[
        'At a crossing h=k=0, test h=l, h=g, h=n by the three arithmetic gates serialized above.',
        'If a candidate on the current branch is consecutive, follow its other boundary line first, then return to the current line.',
        'If no crossing occurs before the branch endpoint, continue through the fissure along the paired boundary line.',
        'A complete field boundary must reconnect to its initial boundary line.',
        'Development may stop only after numerical repetition is detected and a stabilizer word is serialized.'
      ],
      'total_fig1_compile':'NT: fissure-ordering/branch-order and total repetition closure are not yet serialized',
      'historical_H1_H2':'NT_GUARDED'
    }

def selftest():
    checks={}
    def ck(k,v): checks[k]=bool(v)
    ck('twelve_moves',len(S)==12 and len({tuple(map(int,list(x))) for x in S})==12)
    ck('all_moves_involutions',all(A*A==sp.eye(3) for A in S))
    ck('all_moves_det_minus1',all(A.det()==-1 for A in S))
    ck('relabel_count24',len(RELABEL)==24)
    # Verify wall-flip typing symbolically.
    pp=sp.symbols('p01 p02 p03 p12 p13 p23')
    a=pp[0]+pp[3]+pp[4]; b=pp[1]+pp[3]+pp[5]; c=pp[2]+pp[4]+pp[5]
    k=-pp[3]; h=-pp[4]; g=-pp[5]
    G=gram_from_six(a,b,c,g,h,k)
    for wall,(i,j) in WALL_TO_MOVES.items():
        pos=WALLS.index(wall)
        for idx in (i,j):
            cp=list(conorms(S[idx-1].T*G*S[idx-1]).values())
            ck(f'{wall}_S{idx}_flips',sp.expand(cp[pos]+pp[pos])==0)
    # Positive reduced synthetic chart from positive conorms.
    pvals=dict(zip(WALLS,[1,2,3,4,5,6]))
    Gred=gram_from_six(pvals['p01']+pvals['p12']+pvals['p13'],
                       pvals['p02']+pvals['p12']+pvals['p23'],
                       pvals['p03']+pvals['p13']+pvals['p23'],
                       -pvals['p23'],-pvals['p13'],-pvals['p12'])
    ck('synthetic_positive_definite',all(x>0 for x in Gred.cholesky().diagonal()))
    ck('synthetic_reduced',is_reduced_positive(Gred))
    ck('conorm_roundtrip',conorms(Gred)==pvals)
    # Wall fixture and candidate transitions.
    pwall=dict(pvals);pwall['p12']=0
    Gwall=gram_from_six(pwall['p01']+pwall['p12']+pwall['p13'],
                        pwall['p02']+pwall['p12']+pwall['p23'],
                        pwall['p03']+pwall['p13']+pwall['p23'],
                        -pwall['p23'],-pwall['p13'],-pwall['p12'])
    ck('wall_fixture_exact',active_walls(Gwall)==['p12'])
    ck('wall_candidates_two',len(adjacent_candidates(Gwall,'p12'))==2)
    # Historical seed identities.
    seed=source_known_fig1_seed()
    ck('fig1_invariant60',seed['checks']['det_Q0']==60)
    ck('fig1_T_stabilizer',seed['checks']['T_stabilizes_Q0'])
    ck('fig1_U_stabilizer',seed['checks']['U_stabilizes_Q0'])
    ck('fig1_TU_commute',seed['checks']['TU_commute'])
    # Positive correspondent identity fixture (symbolic).
    xi,eta,zeta=sp.symbols('xi eta zeta')
    qp=positive_correspondent(Q0,xi,eta,zeta)
    ck('positive_pullback_a',sp.expand(qp[0,0]-(2*xi**2-12))==0)
    ck('positive_pullback_b',sp.expand(qp[1,1]-(2*eta**2+1))==0)
    ck('positive_pullback_c',sp.expand(qp[2,2]-(2*zeta**2+5))==0)
    # m+ formula from homogeneous relations.
    hm=homogeneous(qp)['m']
    ck('m_pullback',sp.expand(hm-(-1-2*eta*(xi+eta+zeta)))==0)
    return {'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
            'check_count':len(checks),'fig1_seed':seed,
            'wall_to_moves':{k:[f'S{i}',f'S{j}'] for k,(i,j) in WALL_TO_MOVES.items()},
            'historical_complex_schema':blank_historical_complex_schema()}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--selftest',action='store_true')
    ap.add_argument('--emit-schema')
    ap.add_argument('--v50-source-feasibility',action='store_true')
    ns=ap.parse_args()
    if ns.v50_source_feasibility:
        r=v50_source_feasibility_report(); print(json.dumps(r,indent=2,ensure_ascii=False)); raise SystemExit(0 if r['status']=='PASS' else 1)
    if ns.selftest:
        r=selftest(); print(json.dumps(r,indent=2,ensure_ascii=False)); raise SystemExit(0 if r['status']=='PASS' else 1)
    if ns.emit_schema:
        Path(ns.emit_schema).write_text(json.dumps(blank_historical_complex_schema(),indent=2,ensure_ascii=False),encoding='utf-8')
        print(ns.emit_schema); return
    ap.print_help()

if __name__=='__main__': main()
