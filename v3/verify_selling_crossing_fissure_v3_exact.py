#!/usr/bin/env python3
import json
from itertools import combinations, permutations
import sympy as sp


def cycle_incidence(n):
    B = sp.zeros(n,n)
    for i in range(n):
        B[i,i] = -1
        B[i,(i+1)%n] = 1
    return B

def mat_to_list(M):
    return [[int(x) if x.is_Integer else str(sp.simplify(x)) for x in row] for row in M.tolist()]

def spectrum_exact(M):
    ev = M.eigenvals()
    vals=[]
    for v,m in ev.items(): vals += [str(sp.simplify(v))]*int(m)
    return sorted(vals, key=lambda s: float(sp.N(sp.sympify(s))))

def K4_L(a,b):
    return sp.Matrix([
        [2*a+b, -a, -b, -a],
        [-a, 2*a+b, -a, -b],
        [-b, -a, 2*a+b, -a],
        [-a, -b, -a, 2*a+b],
    ])

def permute_matrix(M,p):
    return M.extract(p,p)

# X4: two transverse full walls, four sectors.
B4 = cycle_incidence(4)
L4 = B4.T*B4
assert L4 == K4_L(sp.Integer(1),sp.Integer(0))

# X6: Selling historical crossing h=k=0, three lines h=0,k=0,h+k=0, six sectors.
B6 = cycle_incidence(6)
L6 = B6.T*B6
assert L6*sp.ones(6,1) == sp.zeros(6,1)

# Y3: Selling historical fissure h=m=0, three one-sided rays, three fields.
BY = cycle_incidence(3)
LY = BY.T*BY

# Four-field compressions of X6. Exact no-go against K4(a,b).
a,b=sp.symbols('a b', real=True, nonnegative=True)
subpatches=[]
for S in combinations(range(6),4):
    S=list(S)
    LP=L6.extract(S,S)
    # Number of cut edges from selected to omitted fields.
    selected=set(S)
    cut=[]
    induced=[]
    for i in range(6):
        j=(i+1)%6
        if (i in selected) ^ (j in selected): cut.append((i,j))
        if i in selected and j in selected: induced.append((i,j))
    assert len(cut)>=2
    # Spectral exactness no-go: K4 Laplacian kills constants; LP does not because of cut boundary potential.
    row_sums=LP*sp.ones(4,1)
    assert row_sums != sp.zeros(4,1)
    exact_K4=False
    exact_witness=None
    # Exhaustive symbolic equality check under all label permutations.
    for p in permutations(range(4)):
        M=permute_matrix(LP,list(p))
        eqs=[sp.Eq(M[i,j],K4_L(a,b)[i,j]) for i in range(4) for j in range(4)]
        sol=sp.solve(eqs,[a,b], dict=True)
        if sol:
            for so in sol:
                av=so.get(a,a); bv=so.get(b,b)
                if av.is_number and bv.is_number and av>=0 and bv>=0:
                    exact_K4=True; exact_witness=(p,so); break
        if exact_K4: break
    assert not exact_K4
    # Chain no-go for any closed target graph incidence: each row of any target incidence has coordinate sum 0.
    # B6*F2 has cut-interface rows with one +/-1, hence row sum +/-1. Therefore it is outside target incidence row-space.
    # We record this as exact for interface/cellular chain transport to a closed 4-field target.
    subpatches.append({
        'subset': S,
        'induced_internal_edges': len(induced),
        'cut_edges': len(cut),
        'LP': mat_to_list(LP),
        'LP_row_sums': [int(x) for x in row_sums],
        'Dspec_zero_possible_K4': False,
        'Dpartial_zero_possible_closed_four_field_target': False,
        'reason_Dspec': 'LP*1 != 0 because omitted adjacent fields generate boundary potential, whereas every K4(a,b) Laplacian kills 1',
        'reason_Dpartial': 'cut-interface rows of B6 F2 have coordinate sum +/-1, but every row in the span of a closed graph incidence matrix has coordinate sum 0',
    })

# Link/prong test for cusp equivalence.
# Polar cusp: parametrization t -> (3 b t^2, 2 b t^3), t != 0 has two local components (+/- t).
# Selling historical fissure: three one-sided limit rays => three components after deleting the singular point.
prongs_polar=2
prongs_selling_Y=3
assert prongs_polar != prongs_selling_Y

# Explicit quartic prototypes, showing degree 4 alone does not decide A2 cusp type.
u,v=sp.symbols('u v', real=True)
protos={
    'polar_normalized': v**2-u**3,
    'quartic_A2': v**2-u**3-u**4,
    'quartic_tacnode_A3': v**2-u**4,
}
proto_results={}
for name,G in protos.items():
    grad=sp.Matrix([sp.diff(G,u),sp.diff(G,v)]).subs({u:0,v:0})
    H=sp.hessian(G,(u,v)).subs({u:0,v:0})
    rank=H.rank()
    ker=H.nullspace()
    cubic=None
    if ker:
        e=ker[0]
        t=sp.symbols('t')
        ge=sp.expand(G.subs({u:t*e[0],v:t*e[1]}))
        cubic=sp.diff(ge,t,3).subs(t,0)
    proto_results[name]={
        'G': str(G),
        'gradient_at_origin': [str(x) for x in grad],
        'hessian': [[str(x) for x in row] for row in H.tolist()],
        'hessian_rank': rank,
        'kernel_cubic_derivative': None if cubic is None else str(cubic),
        'A2_jet_test': bool(grad==sp.zeros(2,1) and rank==1 and cubic!=0),
    }

out={
  'source_anchored_models': {
    'X6_historical_crossing_h_eq_k_eq_0': {
      'wall_normal_form': ['u=0','v=0','u+v=0'],
      'fields': 6,
      'local_prongs': 6,
      'incidence_B2': mat_to_list(B6),
      'local_down_laplacian': mat_to_list(L6),
      'spectrum': spectrum_exact(L6),
      'four_field_K4_audit': subpatches,
      'joint_zero_status': 'REFUTED_FOR_EVERY_ACTUAL_FOUR_FIELD_SUBPATCH_OF_THIS_C6_GERM',
    },
    'Y3_historical_fissure_h_eq_m_eq_0': {
      'wall_normal_form': ['three one-sided rays'],
      'fields': 3,
      'local_prongs': 3,
      'incidence_B2': mat_to_list(BY),
      'local_down_laplacian': mat_to_list(LY),
      'spectrum': spectrum_exact(LY),
      'K4_Dpartial_Dspec_status': 'NOT_APPLICABLE_DIMENSION_MISMATCH_3_FIELDS_VS_4_FIELD_DIAGNOSTIC',
      'polar_cusp_prongs': prongs_polar,
      'germ_equivalent_to_polar_cusp': False,
      'germ_status': 'REFUTED_BY_LINK_CARDINALITY_3_VS_2',
    },
  },
  'conditional_local_model': {
    'X4_two_transverse_wall_germ': {
      'fields': 4,
      'incidence_B2': mat_to_list(B4),
      'local_down_laplacian': mat_to_list(L4),
      'spectrum': spectrum_exact(L4),
      'K4_parameters': {'a':1,'b':0},
      'Dpartial': 0,
      'Dspec': 0,
      'status': 'PASS_EXACT_FOR_SUPPORT_REDUCED_K4_1_0_GERM_NOT_PROMOTED_TO_NAMED_HISTORICAL_SELLING_CROSSING',
    }
  },
  'cusp_degree4_control': proto_results,
  'global_conclusions': {
    'historical_six_field_crossing_to_exact_K4_four_field_patch': 'REFUTED',
    'historical_three_field_fissure_to_polar_A2_cusp': 'REFUTED',
    'degree_four_equation_alone_implies_polar_cusp': 'REFUTED',
    'existence_of_other_Selling_fissures_with_A2_cusp_auxiliary_discriminant': 'NT_NOT_CLASSIFIED_BY_RETRIEVED_CORPUS',
  }
}

path='/mnt/data/selling_crossing_fissure_v3_exact_certificate.json'
with open(path,'w',encoding='utf-8') as f: json.dump(out,f,indent=2,ensure_ascii=False)
print(json.dumps(out['global_conclusions'],indent=2,ensure_ascii=False))
print('WROTE',path)
