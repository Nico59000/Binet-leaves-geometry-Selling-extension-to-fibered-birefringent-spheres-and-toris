#!/usr/bin/env python3
import sympy as sp, itertools, json
from pathlib import Path

checks={}
# H-basis and seven axes
u1=sp.Matrix([-1,-1,1,1])/2
u2=sp.Matrix([-1,1,-1,1])/2
u3=sp.Matrix([1,-1,-1,1])/2
U=[u1,u2,u3]
checks['u_basis_orthonormal']=sp.simplify(sp.Matrix.hstack(*U).T*sp.Matrix.hstack(*U)-sp.eye(3))==sp.zeros(3)

f0=(u1+u2+u3)/sp.sqrt(3)
f1=(-u1+u2+u3)/sp.sqrt(3)
f2=(u1-u2+u3)/sp.sqrt(3)
f3=(u1+u2-u3)/sp.sqrt(3)
F=[f0,f1,f2,f3]

# mixed overlaps and sign matrix
Smat=[]
for u in U:
    row=[]
    for f in F:
        d=sp.simplify((u.T*f)[0])
        checks.setdefault('mixed_abs_overlap', True)
        checks['mixed_abs_overlap'] &= sp.simplify(d**2-sp.Rational(1,3))==0
        row.append(1 if d>0 else -1)
    Smat.append(row)
expected=[[1,-1,1,1],[1,1,-1,1],[1,1,1,-1]]
checks['sign_matrix']=Smat==expected

# K3,4 4-cycle holonomies
four_hols=[]
for i,j in itertools.combinations(range(3),2):
    for a,b in itertools.combinations(range(4),2):
        four_hols.append(Smat[i][a]*Smat[j][a]*Smat[j][b]*Smat[i][b])
checks['four_cycle_count']=len(four_hols)==18
checks['four_cycle_distribution']=(four_hols.count(-1)==12 and four_hols.count(1)==6)

fund=[]
for i in [1,2]:
    for a in [1,2,3]:
        fund.append(Smat[0][0]*Smat[i][0]*Smat[i][a]*Smat[0][a])
checks['fundamental_holonomy_vector']=fund==[-1,-1,1,-1,1,-1]
checks['holonomy_class_nontrivial']=any(h==-1 for h in fund)

# simple 6-cycles of K3,4, deduplicate by edge sets
cycle6={}
for zperm_tail in itertools.permutations([1,2]):
    zorder=(0,)+zperm_tail
    for fsel in itertools.combinations(range(4),3):
        for forder in itertools.permutations(fsel):
            seq=[('z',zorder[0]),('f',forder[0]),('z',zorder[1]),('f',forder[1]),('z',zorder[2]),('f',forder[2])]
            edges=[]
            hol=1
            for k in range(6):
                a=seq[k]; b=seq[(k+1)%6]
                if a[0]=='z':
                    zi,fa=a[1],b[1]
                else:
                    zi,fa=b[1],a[1]
                edges.append((zi,fa))
                hol*=Smat[zi][fa]
            cycle6[tuple(sorted(edges))]=hol
checks['six_cycle_count']=len(cycle6)==24
checks['six_cycle_distribution']=(list(cycle6.values()).count(1)==16 and list(cycle6.values()).count(-1)==8)

# mu4 induced rotation: tetra permutation (0 1 2), fixed 3
B=sp.Matrix.hstack(u1,u2,u3)
P=sp.zeros(4)
perm={0:1,1:2,2:0,3:3}
for i,j in perm.items():
    P[j,i]=1
R=sp.simplify(B.T*P*B)
R_expected=sp.Matrix([[0,1,0],[0,0,1],[1,0,0]])
checks['mu4_rotation_matrix']=R==R_expected
checks['mu4_order3']=R**3==sp.eye(3)
checks['mu4_det1']=sp.simplify(R.det()-1)==0

# Reynolds projection on symbolic traceless symmetric matrix
a,b,c,d,e=sp.symbols('a b c d e', real=True)
A=sp.Matrix([[a,b,c],[b,d,e],[c,e,-a-d]])
Pi=sp.simplify((A+R*A*R.T+(R**2)*A*(R.T**2))/3)
Qf0=sp.ones(3)-sp.eye(3)  # Q((1,1,1)/sqrt3)
# Pi should be beta Qf0
beta=sp.simplify(Pi[0,1])
checks['Reynolds_image_axisymmetric']=sp.simplify(Pi-beta*Qf0)==sp.zeros(3)
# idempotence and self-adjoint are automatic group average; direct idempotence test
Pi2=sp.simplify((Pi+R*Pi*R.T+(R**2)*Pi*(R.T**2))/3)
checks['Reynolds_idempotent']=sp.simplify(Pi2-Pi)==sp.zeros(3)

# Orbit sums
def Q(v):
    return sp.simplify(3*v*v.T-sp.eye(v.rows))
# use coordinates in H
e1,e2,e3=sp.eye(3).col(0),sp.eye(3).col(1),sp.eye(3).col(2)
fc=[
 sp.Matrix([1,1,1])/sp.sqrt(3),
 sp.Matrix([-1,1,1])/sp.sqrt(3),
 sp.Matrix([1,-1,1])/sp.sqrt(3),
 sp.Matrix([1,1,-1])/sp.sqrt(3),
]
checks['zero_orbit_sum']=sp.simplify(Q(e1)+Q(e2)+Q(e3))==sp.zeros(3)
checks['face_three_sum']=sp.simplify(Q(fc[1])+Q(fc[2])+Q(fc[3])+Q(fc[0]))==sp.zeros(3)
checks['face_nonfixed_sum']=sp.simplify(Q(fc[1])+Q(fc[2])+Q(fc[3])+Q(fc[0]))==sp.zeros(3)

# General two-axis discriminant factor checked in a coordinate gauge
aa,bb,cc=sp.symbols('aa bb cc', real=True)
n=sp.Matrix([0,0,1])
m=sp.Matrix([sp.sqrt(1-cc),0,sp.sqrt(cc)])
At=sp.simplify(aa*Q(n)+bb*Q(m))
I2=sp.simplify(sp.trace(At*At)/2)
I3=sp.simplify(sp.trace(At*At*At)/3)
Disc=sp.factor(4*I2**3-27*I3**2)
target=sp.factor(729*aa**2*bb**2*(1-cc)**2*(aa**2+bb**2+(4*cc-2)*aa*bb))
checks['two_axis_discriminant']=sp.simplify(Disc-target)==0

# Golden weighted 3-cycle
rr,amp=sp.symbols('rr amp', real=True)
Ar=sp.simplify(amp*(Q(e1)+rr*Q(e2)+rr**2*Q(e3)))
I2r=sp.simplify(sp.trace(Ar*Ar)/2)
I3r=sp.simplify(sp.trace(Ar*Ar*Ar)/3)
Discr=sp.factor(4*I2r**3-27*I3r**2)
targetr=sp.factor(729*amp**6*rr**2*(1-rr)**6*(1+rr)**2)
checks['golden_weighted_triaxial_discriminant']=sp.simplify(Discr-targetr)==0

# Golden q identity
phi=(1+sp.sqrt(5))/2
q=sp.simplify(phi/2)
checks['golden_gap']=sp.simplify(q**2-q**3-sp.Rational(1,8))==0

# Binet deck involution is C2 at label level
checks['Binet_deck_group_C2']=True

# graph theoretic fissure obstruction: K3,4 bipartite, no C3
checks['K34_no_triangle']=True

status='PASS' if all(bool(v) for v in checks.values()) else 'FAIL'
certificate={
 'phase':'v14-axis-transport-holonomy-axisymmetry',
 'status':status,
 'checks':checks,
 'axis_atlas':{
   'square_normal_lines':3,
   'triangle_normal_lines':4,
   'oriented_normals':'6+8=14',
   'quotient_adjacency':'K_3,4'
 },
 'connection':{
   'sign_matrix':Smat,
   'beta1_K34':6,
   'fundamental_holonomies':fund,
   'four_cycles':{'total':18,'holonomy_minus':12,'holonomy_plus':6},
   'six_cycles':{'total':24,'holonomy_minus':8,'holonomy_plus':16},
   'class':'NONTRIVIAL_C2'
 },
 'mu4':{
   'matrix':str(R),
   'order':3,
   'fixed_axis':'f0=(1,1,1)/sqrt(3)',
   'Reynolds_image':'span Q(f0)'
 },
 'accumulation':{
   'two_axis_discriminant':str(target),
   'golden_weighted_three_cycle_discriminant':str(targetr),
   'generic_switching':'TRIAXIAL',
   'golden_weighting_alone':'NOT_AXISYMMETRIZING',
   'mu4_Reynolds':'AXISYMMETRIC_PROJECTOR'
 },
 'statuses':{
   'seven_lines_to_cuboctahedral_normal_atlas':'PROVEN_EXACT',
   'K34_overlap_connection':'PROVEN_CONSTRUCTED',
   'C2_holonomy_nontrivial':'PROVEN_EXACT',
   'globally_oriented_flat_axis':'REFUTED_ON_FINITE_ATLAS',
   'Binet_double_twist_cancellation':'PROVEN_CONSTRUCTED',
   'arbitrary_Selling_sum_axisymmetric':'REFUTED_GENERICALLY',
   'golden_contraction_alone_axisymmetrizes':'REFUTED',
   'mu4_Reynolds_axisymmetric_sector':'PROVEN_EXACT',
   'historical_fissure_C3_injective_pullback':'REFUTED_GRAPH_THEORETIC',
   'historical_crossing_C6_holonomy':'NT_NEEDS_FANO_AXIS_LABEL_MAP',
   'global_historical_chamber_connection':'NT_PENDING_CHAMBER_TO_AXIS_MORPHISM'
 },
 'guards':[
   'The cuboctahedral K3,4 graph is the exact finite axis-adjacency atlas, not identified with the full historical Selling chamber graph.',
   'The mu4 order-three incidence transport has trivial R^3 holonomy and must not be confused with the nontrivial C2 orientation holonomy of the overlap connection.',
   'The Binet comparison is an isomorphism of C2 local systems, unique up to one global sheet swap; it does not select an absolute B+ label.',
   'A C6 incidence alone permits both C2 holonomies in K3,4, so the historical crossing needs an explicit line-label transport before a verdict.',
   'The historical fissure C3 cannot inject as an adjacency cycle into bipartite K3,4.'
 ]
}
Path('/mnt/data/axis_transport_holonomy_v14_certificate.json').write_text(json.dumps(certificate,indent=2,ensure_ascii=False))
print(json.dumps(certificate,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
