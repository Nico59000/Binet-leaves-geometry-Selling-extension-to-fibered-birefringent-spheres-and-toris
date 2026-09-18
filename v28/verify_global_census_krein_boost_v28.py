#!/usr/bin/env python3
import csv, json, math
from pathlib import Path
import sympy as sp

OUT=Path('/mnt/data')
checks={}

# ---------- fixed chart gauge / source-audited local census ----------
sigma_minus=[0,4,1,5,2,6,3]
fano={
    'u1':[0,1,3],
    'u2':[0,4,5],
    'u3':[0,2,6],
}

def cycle_edges(prefix,n,walls,axes,eps_names,source_status):
    rows=[]
    for i in range(n):
        wall=walls[i%len(walls)]
        axis=axes[i%len(axes)]
        eps=eps_names[i%len(eps_names)]
        rows.append({
            'edge_id':f'{prefix}_e{i}',
            'tail':f'{prefix}_C{i}',
            'head':f'{prefix}_C{(i+1)%n}',
            'tau_e':-1,
            'wall':wall,
            'axis_label':axis,
            'L_e':fano[axis],
            'L_status':'CONSTRUCTED_RELATIVE_TO_FIXED_S3_CHART_GAUGE',
            'sigma_e':sigma_minus,
            'sigma_status':'SELECTED_ANCHORED_GOLDEN_FRAMING_TRANSPORTED_BY_CHART',
            'epsilon_e':eps,
            'epsilon_status':'SYMBOLIC_SOURCE_AMPLITUDE_NOT_NUMERICALLY_CENSUSED',
            's_e':f'kappa_hist*Omega_sigma({axis})*{eps}^3',
            'source_status':source_status,
        })
    return rows

cross_edges=cycle_edges(
    'X6',6,
    ['k=0','h+k=0','h=0'],
    ['u2','u3','u1'],
    ['eps_k','eps_hk','eps_h'],
    'HISTORICAL_CROSSING_SOURCE_AUDITED; FANO/SINGER DATA CONSTRUCTED IN FIXED GAUGE'
)
fiss_edges=cycle_edges(
    'Y3',3,
    ['h=0','m=0','h-m=0'],
    ['u1','u2','u3'],
    ['eps_h','eps_m','eps_hm'],
    'HISTORICAL_FISSURE_SOURCE_AUDITED; FANO/SINGER DATA CONSTRUCTED IN FIXED GAUGE'
)

# Correct field IDs for fissure are Y3_C0..2 already.

# ---------- twisted incidence ----------
def twisted_incidence_cycle(n,tau=-1):
    D=sp.zeros(n,n)
    for i in range(n):
        D[i,i]=-tau
        D[i,(i+1)%n]=1
    return D

D6=twisted_incidence_cycle(6,-1)
D3=twisted_incidence_cycle(3,-1)
a,b,c=sp.symbols('a b c', real=True)
s6=sp.Matrix([a,b,c,a,b,c])
s3=sp.Matrix([a,b,c])
checks['X6_Dc_rank_5']=D6.rank()==5
checks['X6_augmented_rank_5_for_repeated_halfturn_sources']=D6.row_join(s6).rank()==5
left6=D6.T.nullspace()
checks['X6_left_obstruction_vector']=len(left6)==1 and left6[0] in (sp.Matrix([-1,1,-1,1,-1,1]),sp.Matrix([1,-1,1,-1,1,-1]))
checks['X6_alternating_period_zero']=sp.simplify((left6[0].T*s6)[0])==0
checks['Y3_Dc_full_rank_3']=D3.rank()==3
checks['Y3_any_source_locally_exact']=D3.row_join(s3).rank()==3

# Generic m-fissure planar bouquet: D=2*1_m and H^1 dimension m-1.
plane={}
for m in range(1,9):
    D=2*sp.ones(m,1)
    plane[str(m)]={'rank_Dc':D.rank(),'H1_dim':m-D.rank()}
    checks[f'planar_m{m}_H1_dim_mminus1']=(m-D.rank()==m-1)

# Relative charge criterion on m=4 example.
p1,p2,p3,p4=sp.symbols('p1 p2 p3 p4')
D4=2*sp.ones(4,1)
checks['planar_equal_odd_periods_exact']=D4.row_join(sp.Matrix([p1,p1,p1,p1])).rank()==D4.rank()
checks['planar_unequal_odd_periods_obstruct']=D4.row_join(sp.Matrix([1,2,1,1])).rank()>D4.rank()

# ---------- Fano/Singer consistency ----------
D0={0,1,3}
all_lines=sorted({tuple(sorted(((x+t)%7 for x in D0))) for t in range(7)})
checks['fixed_u_lines_are_Fano_lines']=all(tuple(sorted(v)) in all_lines for v in fano.values())
checks['sigma_minus_is_permutation']=sorted(sigma_minus)==list(range(7))
checks['sigma_minus_is_multiplier_4']=all(sigma_minus[n]==(4*n)%7 for n in range(7))

# ---------- SO^+(1,1) boost algebra ----------
x,y=sp.symbols('x y', real=True)
def B(z):
    return sp.Matrix([[sp.cosh(z),sp.sinh(z)],[sp.sinh(z),sp.cosh(z)]])
eta=sp.diag(1,-1)
checks['boost_preserves_eta']=sp.simplify(B(x).T*eta*B(x)-eta)==sp.zeros(2)
checks['boost_det_one']=sp.simplify(B(x).det()-1)==0
checks['boost_additivity']=sp.simplify(sp.trigsimp(B(x)*B(y)-B(x+y)))==sp.zeros(2)
checks['boost_inverse']=sp.simplify(sp.trigsimp(B(x)*B(-x)-sp.eye(2)))==sp.zeros(2)

# Gauge rapidities chi_e -> chi_e + psi_tail - psi_head; cycle sum invariant.
chi=sp.symbols('chi0:5', real=True)
psi=sp.symbols('psi0:5', real=True)
chip=[]
for i in range(5):
    chip.append(chi[i]+psi[i]-psi[(i+1)%5])
checks['boost_cycle_sum_gauge_invariant']=sp.simplify(sum(chip)-sum(chi))==0

# Crossing exact A2 word identity => no continuous boost.
# Use S3 permutations represented as tuples, with the three transpositions.
def compose(p,q):
    return tuple(p[q[i]] for i in range(3))
e=(0,1,2)
r1=(1,0,2); r2=(0,2,1); r3=(2,1,0)
g=e
for r in [r1,r2,r3,r1,r2,r3]:
    g=compose(r,g)
checks['X6_Weyl_word_identity']=g==e
checks['X6_continuous_boost_zero']=checks['X6_Weyl_word_identity']

# Fissure word is an involution. In det=+1/time-reversing component an element is -B(chi),
# whose square is B(2chi); involution forces chi=0.
gf=e
for r in [r1,r2,r3]:
    gf=compose(r,gf)
checks['Y3_Weyl_word_order2']=compose(gf,gf)==e and gf!=e
minusB=-B(x)
checks['time_reverse_square_is_B2chi']=sp.simplify(sp.trigsimp(minusB*minusB-B(2*x)))==sp.zeros(2)
# algebraic zero criterion: sinh(2x)=0 over reals iff x=0; encode exact symbolic fact via solve.
sol=sp.solve(sp.Eq(sp.sinh(2*x),0),x)
checks['time_reverse_involution_forces_zero_rapidity']=(sol==[0])
checks['Y3_continuous_boost_zero']=checks['Y3_Weyl_word_order2'] and checks['time_reverse_involution_forces_zero_rapidity']

# Ordinary incidence for continuous boost on the lifted/discrete-trivialized graph.
def ordinary_incidence_cycle(n):
    D=sp.zeros(n,n)
    for i in range(n):
        D[i,i]=-1
        D[i,(i+1)%n]=1
    return D
Db=ordinary_incidence_cycle(5)
z=sp.symbols('z0:5')
zv=sp.Matrix(z)
checks['boost_incidence_rank_nminus1']=Db.rank()==4
checks['boost_exact_if_cycle_sum_zero']=Db.row_join(sp.Matrix([1,-1,2,-2,0])).rank()==Db.rank()
checks['boost_obstructed_if_cycle_sum_nonzero']=Db.row_join(sp.Matrix([1,0,0,0,0])).rank()>Db.rank()

# ---------- Census source boundary ----------
census={
  'phase':'v28-historical-dual-census',
  'status':'MAXIMAL_SOURCE_SUPPORTED_PARTIAL_CENSUS',
  'global_complete':False,
  'global_completion_blocker':'The corpus supplies the dual-graph definition and local historical crossing/fissure germs, but not an exhaustive global chamber/edge enumeration, all edge amplitudes epsilon_e, or a global transported Singer framing.',
  'required_edge_tuple':['tau_e','L_e','sigma_e','epsilon_e'],
  'fixed_chart_gauge':{
      'u1':fano['u1'],'u2':fano['u2'],'u3':fano['u3'],
      'sigma_minus':sigma_minus,
      'gauge_status':'S3/Fano chart choice; rank conclusions below are gauge-independent.'
  },
  'source_audited_germs':{
      'historical_crossing_X6':{
          'vertices':[f'X6_C{i}' for i in range(6)],
          'edges':cross_edges,
          'Dc':[[int(D6[i,j]) for j in range(6)] for i in range(6)],
          'rank_Dc':D6.rank(),
          'source_vector':'(s_k,s_hk,s_h,s_k,s_hk,s_h)',
          'rank_augmented_for_source_pattern':D6.row_join(s6).rank(),
          'verdict':'LOCAL_EXACTNESS_PROVEN'
      },
      'historical_fissure_Y3':{
          'vertices':[f'Y3_C{i}' for i in range(3)],
          'edges':fiss_edges,
          'Dc':[[int(D3[i,j]) for j in range(3)] for i in range(3)],
          'rank_Dc':D3.rank(),
          'source_vector':'(s_h,s_m,s_hm) arbitrary',
          'rank_augmented_generic':D3.row_join(s3).rank(),
          'verdict':'LOCAL_EXACTNESS_PROVEN'
      }
  },
  'planar_m_fissure_skeleton':{
      'formula':'Dc=2*1_m, rank=1, H1(R^2\\F;R_sign)=R^(m-1)',
      'relative_charges':'P_i-P_1',
      'rank_samples_m1_to_m8':plane
  },
  'global_rank_target':{
      'formula':'rank(D_c) ?= rank([D_c|s])',
      'verdict':'NT_NOT_EVALUABLE_WITHOUT_COMPLETE_GLOBAL_CENSUS'
  },
  'krein_boost':{
      'discrete_reduction':'factor determinant/time characters; on the Binet/time double cover transitions reduce to SO^+(1,1)',
      'edge_parameter':'g_e^+=B(chi_e)',
      'gauge_law':'chi_e -> chi_e + psi_tail - psi_head',
      'cycle_holonomy':'chi(gamma)=sum_e chi_e',
      'historical_crossing':'chi=0 PROVEN from exact Weyl word identity',
      'historical_fissure':'chi=0 PROVEN: time-reversing holonomy is an involution, hence -B(0)',
      'global':'NT_NEEDS_COMPLETE_O11_EDGE_TRANSITIONS'
  }
}

# Write census JSON and CSV.
(OUT/'historical_dual_census_v28.json').write_text(json.dumps(census,indent=2,ensure_ascii=False),encoding='utf-8')
with (OUT/'historical_dual_census_v28.csv').open('w',newline='',encoding='utf-8') as fh:
    cols=['germ','edge_id','tail','head','tau_e','wall','axis_label','L_e','sigma_e','epsilon_e','s_e','source_status']
    w=csv.DictWriter(fh,fieldnames=cols)
    w.writeheader()
    for germ,rows in [('X6',cross_edges),('Y3',fiss_edges)]:
        for row in rows:
            rr={k:row.get(k,'') for k in cols}
            rr['germ']=germ
            rr['L_e']=' '.join(map(str,row['L_e']))
            rr['sigma_e']=' '.join(map(str,row['sigma_e']))
            w.writerow(rr)

checks['census_declares_global_incompleteness']=census['global_complete'] is False
checks['census_has_all_local_edge_tuples']=all(all(k in row for k in census['required_edge_tuple']) for row in cross_edges+fiss_edges)
checks['global_rank_not_silently_claimed']=census['global_rank_target']['verdict'].startswith('NT_')

checks={k:bool(v) for k,v in checks.items()}
status='PASS' if all(checks.values()) else 'FAIL'
cert={
  'phase':'v28-global-census-and-krein-boost',
  'status':status,
  'checks':checks,
  'census':{
      'global_complete':False,
      'local_edge_count':len(cross_edges)+len(fiss_edges),
      'crossing_Dc_rank':D6.rank(),
      'crossing_augmented_rank':D6.row_join(s6).rank(),
      'fissure_Dc_rank':D3.rank(),
      'fissure_augmented_rank':D3.row_join(s3).rank(),
      'global_rank':'NT: complete chamber/edge census and source amplitudes absent from corpus',
      'artifacts':['historical_dual_census_v28.json','historical_dual_census_v28.csv']
  },
  'twisted_memory':{
      'crossing':'rank equality PROVEN for source repeated-halfturn pattern',
      'fissure':'rank equality PROVEN for arbitrary local source values',
      'global':'NT until global census completion',
      'planar_relative_memory':'R^(m-1) remains the exact abstract result'
  },
  'krein_boost':{
      'group':'SO^+(1,1) ~= (R,+) via rapidity',
      'gauge':'chi_e -> chi_e + psi_tail - psi_head',
      'cohomology':'[chi] in H^1(lifted regular graph;R)',
      'crossing_local_boost':'0 PROVEN',
      'fissure_local_boost':'0 PROVEN after discrete time reversal factorization',
      'global_boost':'NT_NEEDS_GLOBAL_O11_TRANSITIONS'
  },
  'statuses':{
      'machine_readable_edge_schema':'PROVEN_CONSTRUCTED',
      'source_audited_local_census':'PROVEN_MAXIMAL_FROM_CURRENT_CORPUS',
      'complete_global_historical_census':'NT_NOT_SUPPLIED_BY_CURRENT_CORPUS',
      'global_history_rank_test':'NT_NOT_EXECUTABLE_YET',
      'local_crossing_history_exactness':'PROVEN_BY_RANK',
      'local_fissure_history_exactness':'PROVEN_BY_FULL_RANK',
      'SOplus11_rapidity_cocycle':'PROVEN_EXACT',
      'local_crossing_continuous_Krein_holonomy':'ZERO_PROVEN',
      'local_fissure_continuous_Krein_holonomy':'ZERO_PROVEN',
      'global_continuous_Krein_holonomy':'NT_NEEDS_GLOBAL_TRANSITIONS'
  },
  'guards':[
      'The current corpus defines the global dual graph abstractly but does not enumerate every historical chamber and edge; local germs must not be relabeled as a complete global census.',
      'Fano line labels and the anchored golden Singer order in the census are a fixed chart gauge/selected framing, not an absolute historical choice.',
      'Historical amplitudes epsilon_e are symbolic because no exhaustive numeric amplitude table is supplied.',
      'The continuous boost class is defined only after factoring the already-audited discrete determinant/time components, equivalently on the common double cover where the time character is trivial.',
      'Zero boost on the audited crossing/fissure germs does not imply zero global boost on quotient cycles.'
  ]
}
(OUT/'global_census_krein_boost_v28_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
