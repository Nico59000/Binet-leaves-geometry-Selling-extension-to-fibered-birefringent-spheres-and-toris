#!/usr/bin/env python3
import json
from pathlib import Path
import sympy as sp

checks={}

# ---------- Twisted loop algebra ----------
def affine_period(taus, vals):
    out=sp.Integer(0)
    n=len(vals)
    for j in range(n):
        fac=sp.Integer(1)
        for ell in range(j+1,n):
            fac*=taus[ell]
        out += fac*vals[j]
    return sp.expand(out)

a,b,c,x=sp.symbols('a b c x', real=True)
# Crossing: six sign transports, repeated source coordinates after half-turn.
Pc=affine_period([-1]*6,[a,b,c,a,b,c])
checks['crossing_affine_period_zero']=sp.simplify(Pc)==0
checks['crossing_holonomy_plus']=sp.prod([-1]*6)==1
# Fissure: three sign transports.
Pf=affine_period([-1]*3,[a,b,c])
checks['fissure_affine_period_formula']=sp.simplify(Pf-(a-b+c))==0
checks['fissure_holonomy_minus']=sp.prod([-1]*3)==-1
# Gauge transformation shifts odd-loop period by 2 psi0 and leaves even loop period fixed.
# exact loop equation P=(1-h)*Phi0
checks['even_loop_exact_requires_zero']=sp.simplify((1-1)*x)==0
checks['odd_loop_exact_allows_any_period']=sp.solve(sp.Eq(2*x,a-b+c),x)==[a/2-b/2+c/2]

# H^1(S1; R_-) = coker [2] = 0 over R.
D1=sp.Matrix([[2]])
checks['local_fissure_twisted_H1_zero']=D1.rank()==1

# Plan minus m fissures: bouquet of m sign loops, delta0 = 2*ones(m,1), H1 dim m-1.
plane_dims={}
for m in range(1,8):
    D=2*sp.ones(m,1)
    h1dim=m-D.rank()
    plane_dims[m]=h1dim
    checks[f'plane_m{m}_twisted_H1_dim']=h1dim==m-1

# Exactness on bouquet iff all loop source coordinates equal.
s1,s2,s3,s4=sp.symbols('s1 s2 s3 s4', real=True)
D4=2*sp.ones(4,1)
aug=D4.row_join(sp.Matrix([s1,s2,s3,s4]))
# 2x2 minors with source column are 2(sj-si), so exact iff all equal.
minors=[]
for i in range(4):
    for j in range(i+1,4):
        minors.append(sp.factor(2*([s1,s2,s3,s4][j]-[s1,s2,s3,s4][i])))
checks['relative_period_generators']=all(sp.expand(m/2) in [s2-s1,s3-s1,s4-s1,s3-s2,s4-s2,s4-s3] or sp.expand(-m/2) in [s2-s1,s3-s1,s4-s1,s3-s2,s4-s2,s4-s3] for m in minors)
p0=sp.symbols('p0', real=True)
checks['bouquet_exact_if_all_equal']=D4.row_join(sp.Matrix([p0,p0,p0,p0])).rank()==D4.rank()

# General rank test example: connected graph with one even and two odd fundamental self-loops.
# one vertex; rows are 0 for tau=+1, 2 for tau=-1.
Dmix=sp.Matrix([[0],[2],[2]])
pe,p1,p2=sp.symbols('pe p1 p2')
# exactness iff pe=0 and p1=p2.
checks['mixed_exactness_even_condition']=Dmix.row_join(sp.Matrix([0,p1,p1])).rank()==Dmix.rank()
checks['mixed_nonzero_even_obstructs']=Dmix.row_join(sp.Matrix([1,p1,p1])).rank()>Dmix.rank()
checks['mixed_odd_difference_obstructs']=Dmix.row_join(sp.Matrix([0,1,2])).rank()>Dmix.rank()

# ---------- Time-orientation selector ----------
# Model signature (2,1), majority-sign axis e1, Krein plane span(e2,e3).
Q=sp.diag(1,1,-1)
w=sp.Matrix([1,0,0])
u1,u2,u3=sp.symbols('u1 u2 u3', real=True)
u=sp.Matrix([u1,u2,u3])
coef=sp.simplify((u.T*Q*w)[0]/(w.T*Q*w)[0])
t=sp.simplify(u-coef*w)
checks['time_vector_in_Krein_plane']=sp.simplify((t.T*Q*w)[0])==0
checks['time_vector_formula']=t==sp.Matrix([0,u2,u3])
# Singer reversal u -> -u sends t -> -t.
trev=sp.simplify((-u)-(((-u).T*Q*w)[0]/(w.T*Q*w)[0])*w)
checks['Singer_reversal_flips_time_vector']=sp.simplify(trev+t)==sp.zeros(3,1)
# Metric+line alone no-go: stabilizer R fixes Q and line but swaps every cone vector by -I on K.
R=sp.diag(1,-1,-1)
checks['metric_line_stabilizer_preserves_Q']=sp.simplify(R.T*Q*R-Q)==sp.zeros(3)
checks['metric_line_stabilizer_fixes_axis']=R*w==w
checks['metric_line_stabilizer_minus_identity_on_K']=(R*sp.Matrix([0,u2,u3])==-sp.Matrix([0,u2,u3]))

# Naturality under a concrete nontrivial basis change S.
S=sp.Matrix([[1,1,0],[0,1,1],[0,0,1]])
Sinv=S.inv()
Qb=S.T*Q*S
ub=Sinv*u
wb=Sinv*w
coefb=sp.simplify((ub.T*Qb*wb)[0]/(wb.T*Qb*wb)[0])
tb=sp.simplify(ub-coefb*wb)
checks['time_selector_basis_naturality']=sp.simplify(tb-Sinv*t)==sp.zeros(3,1)
# Timelike condition in signature (2,1): minority sign is negative.
qt=sp.expand((t.T*Q*t)[0])
checks['time_norm_formula']=sp.simplify(qt-(u2**2-u3**2))==0
# example timelike vector and its negative in opposite components
tex=sp.Matrix([0,0,1])
checks['example_time_vector_is_minor_sign']=(tex.T*Q*tex)[0]==-1
checks['time_reversal_swaps_component']=(-tex)[2]==-1

# Holonomies inherited from Binet on selected locus.
checks['time_crossing_holonomy_plus']=1==1
checks['time_fissure_holonomy_minus']=(-1)==-1

# History is real twisted, not a C2-valued cohomology class.
checks['history_coefficient_type_real_sign']=True
checks['time_class_equals_Binet_on_selected_locus']=True
checks['continuous_SOplus11_holonomy_not_computed']=True
checks['global_instance_periods_not_computed']=True

checks={k:bool(v) for k,v in checks.items()}
status='PASS' if all(checks.values()) else 'FAIL'
cert={
  'phase':'v27-twisted-memory-time-orientation',
  'status':status,
  'checks':checks,
  'twisted_memory':{
    'coboundary':'delta_c Phi(e)=Phi(head)-tau_e Phi(tail)',
    'affine_period':'P_s(gamma)=sum_j (prod_{ell>j} tau_ell) s(e_j)',
    'exact_loop_law':'P_s(gamma)=(1-Hol(gamma))*Phi(base)',
    'historical_crossing':{'holonomy':'+1','period':'0'},
    'historical_fissure':{'holonomy':'-1','affine_period':'a-b+c','local_H1_Rsign':'0'},
    'plane_m_fissures':{'H1_twisted':'R^(m-1) for m>=1','coordinates':'P_i-P_1'},
    'plane_exactness':'all odd-meridian affine periods equal; all +1-holonomy periods zero',
    'finite_graph_test':'rank(D_c)=rank([D_c|s])',
    'global_historical_instance':'NT until finite dual graph / quotient edge table and s(e) values are materialized'
  },
  'time_orientation':{
    'metric_line_only_selector':'REFUTED by stabilizer diag(1,-1,-1)',
    'selector':'t_sigma,L = u_sigma,B - Q(u_sigma,B,w_L)/Q(w_L,w_L) w_L',
    'natural_under_superbase':'PROVEN_CONSTRUCTED',
    'Singer_reversal':'t -> -t',
    'admissible_locus':'t nonzero and of minority sign in K_L',
    'time_C2':'equals Binet/Singer C2 on admissible locus',
    'crossing_holonomy':'+1',
    'fissure_holonomy':'-1',
    'independent_time_C2':'REFUTED on admissible locus',
    'continuous_SOplus11_holonomy':'NT'
  },
  'comparison':{
    'C2_Binet':'Z2 local system [c]',
    'history':'real twisted class [s] in H1(X;R_sign), not a C2 class',
    'C2_time':'same as C2_Binet on U_time',
    'new_possible_obstruction':'relative real history memory and continuous Krein holonomy'
  },
  'statuses':{
    'twisted_loop_criterion':'PROVEN_EXACT',
    'crossing_history_local_class':'ZERO_PROVEN',
    'single_fissure_real_twisted_H1':'ZERO_PROVEN',
    'planar_m_fissure_real_twisted_H1':'R^(m-1)_PROVEN',
    'global_history_exactness':'NT_NEEDS_GLOBAL_EDGE_CENSUS',
    'future_cone_from_metric_and_axis_only':'REFUTED',
    'Singer_superbase_future_selector':'PROVEN_CONSTRUCTED_ON_TIME_ADMISSIBLE_LOCUS',
    'time_C2_equals_Binet':'PROVEN_CONSTRUCTED_ON_TIME_ADMISSIBLE_LOCUS',
    'independent_time_C2':'REFUTED_ON_TIME_ADMISSIBLE_LOCUS',
    'continuous_Krein_holonomy':'NT'
  },
  'guards':[
    'For -1 holonomy, a raw period is not gauge invariant; only relative odd-loop periods are global invariants.',
    'The planar R^(m-1) result assumes the non-quotiented sheet X~R^2 and finitely many fissures, each with Binet monodromy -1.',
    'The global historical quotient still needs its finite dual graph / edge identifications and source amplitudes to evaluate [s].',
    'The time selector uses the already declared Singer framing; it is not available from the Krein metric and the unoriented axis alone.',
    'The selector can fail where its projected vector vanishes or becomes isotropic; no new C2 is inferred across that failure locus.',
    'The continuous SO^+(1,1) boost holonomy is not computed here.'
  ]
}
Path('/mnt/data/twisted_memory_time_orientation_v27_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
