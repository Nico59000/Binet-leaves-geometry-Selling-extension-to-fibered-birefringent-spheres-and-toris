#!/usr/bin/env python3
"""Source-typed Selling field-construction primitives for v49.

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
    ns=ap.parse_args()
    if ns.selftest:
        r=selftest(); print(json.dumps(r,indent=2,ensure_ascii=False)); raise SystemExit(0 if r['status']=='PASS' else 1)
    if ns.emit_schema:
        Path(ns.emit_schema).write_text(json.dumps(blank_historical_complex_schema(),indent=2,ensure_ascii=False),encoding='utf-8')
        print(ns.emit_schema); return
    ap.print_help()

if __name__=='__main__': main()
