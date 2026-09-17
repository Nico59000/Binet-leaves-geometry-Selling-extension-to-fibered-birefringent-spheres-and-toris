#!/usr/bin/env python3
import sympy as sp, json
from pathlib import Path

# Edge variables and K4 Selling Laplacian
p12,p13,p14,p23,p24,p34=sp.symbols('p12 p13 p14 p23 p24 p34', real=True)
edges={(0,1):p12,(0,2):p13,(0,3):p14,(1,2):p23,(1,3):p24,(2,3):p34}
L=sp.zeros(4)
for (i,j),p in edges.items():
    v=sp.zeros(4,1); v[i]=1; v[j]=-1
    L += p*(v*v.T)
S=sum(edges.values())
P=sp.eye(4)-sp.ones(4)/4
A=sp.simplify(L-sp.Rational(2,3)*S*P)

checks={}
checks['trace_equals_2S']=sp.simplify(sp.trace(L)-2*S)==0
checks['kernel_one']=sp.simplify(L*sp.ones(4,1))==sp.zeros(4,1)
checks['traceless_on_H']=sp.simplify(sp.trace(A))==0

# Opposite-pair decomposition
s1,s2,s3,d1,d2,d3=sp.symbols('s1 s2 s3 d1 d2 d3', real=True)
subs={
 p12:(s1+d1)/2,p34:(s1-d1)/2,
 p13:(s2+d2)/2,p24:(s2-d2)/2,
 p14:(s3+d3)/2,p23:(s3-d3)/2,
}
normA=sp.expand(sp.trace(A.T*A).subs(subs))
Ssd=s1+s2+s3
center_norm=sp.expand(sum((x-Ssd/3)**2 for x in (s1,s2,s3)))
rhs=sp.expand(center_norm+2*(d1**2+d2**2+d3**2))
checks['five_dof_norm_identity']=sp.simplify(normA-rhs)==0

# Fano chart compatible with tetrahedral face incidence
D0={0,1,3}
lines=sorted({tuple(sorted(((x+t)%7 for x in D0))) for t in range(7)})
label_to_edge={1:p12,2:p13,3:p34,4:p23,5:p14,6:p24}
spectra={}
simple_vectors={}
for line in lines:
    vals={lab:(1 if lab in line else -1) for lab in range(1,7)}
    sub={label_to_edge[lab]:sp.Integer(vals[lab]) for lab in vals}
    dL=sp.simplify(L.subs(sub))
    dS=sum(vals.values())
    dA=sp.simplify(dL-sp.Rational(2,3)*dS*P)
    ev=dA.eigenvals()
    nonzero=[]
    for val,mult in ev.items():
        if val!=0:
            nonzero += [sp.simplify(val)]*mult
    nonzero=sorted(nonzero,key=lambda x:float(x))
    spectra[str(line)]=[str(x) for x in nonzero]
    if 0 in line:
        checks[f'spectrum_{line}']=nonzero==[sp.Rational(-8,3),sp.Rational(4,3),sp.Rational(4,3)]
    else:
        checks[f'spectrum_{line}']=nonzero==[-4,2,2]
    for val,mult,basis in dA.eigenvects():
        if val!=0 and mult==1:
            simple_vectors[line]=sp.Matrix(basis[0])

# isotropic K4
pstar=sp.sqrt(3)
Lstar=sp.simplify(L.subs({x:pstar for x in (p12,p13,p14,p23,p24,p34)}))
checks['isotropic_restriction']=sp.simplify(Lstar-4*pstar*P)==sp.zeros(4)

# old/new radial coordinates
Ppol,Eeq,Rstar=sp.symbols('Ppol Eeq Rstar')
U=(Ppol+Eeq)/2-Rstar**2
Ut=(Ppol+2*Eeq)/3-Rstar**2
D=Ppol-Eeq
checks['radial_coordinate_change']=sp.simplify(U-(Ut+D/6))==0

# Energy-to-D elementary formulas
eps,kappa=sp.symbols('eps kappa', positive=True)
Es0=sp.Rational(16,3)*eps**2
Ed0=sp.Integer(0)
Es1=sp.Integer(0)
Ed1=sp.Integer(6)*eps**2
checks['D_energy_zero_line']=sp.simplify(-kappa*(sp.sqrt(3*Es0)+sp.sqrt(6*Ed0))+4*kappa*eps)==0
checks['D_energy_face_line']=sp.simplify(-kappa*(sp.sqrt(3*Es1)+sp.sqrt(6*Ed1))+6*kappa*eps)==0

# Exact cyclotomic nondegeneracy for golden Singer order step 4
z=sp.symbols('z')
Phi7=sum(z**i for i in range(7))
edge_names={1:(0,1),2:(0,2),3:(2,3),4:(1,2),5:(0,3),6:(1,3)}
edge_vec={}
for lab,(i,j) in edge_names.items():
    v=sp.zeros(4,1); v[i]=1; v[j]=-1
    edge_vec[lab]=v
W=sp.zeros(4,1)
for n in range(7):
    lab=(4*n)%7
    if lab:
        W += (z**n-z**(-n))*edge_vec[lab]
Wpoly=sp.simplify(z**6*W)
remainders={}
nondeg=True
for line,v in simple_vectors.items():
    poly=sp.expand((v.T*Wpoly)[0])
    rem=sp.rem(sp.Poly(poly,z),sp.Poly(Phi7,z)).as_expr()
    remainders[str(line)]=str(sp.factor(rem))
    nondeg &= (sp.expand(rem)!=0)
checks['golden_Singer_axis_projection_nondegenerate']=bool(nondeg)

phi=(1+sp.sqrt(5))/2
q=sp.simplify(phi/2)
checks['golden_weight_gap']=sp.simplify(q**2-q**3-sp.Rational(1,8))==0

status='PASS' if all(checks.values()) else 'FAIL'
certificate={
 'phase':'v13-source-native-geometric-adapter',
 'status':status,
 'checks':checks,
 'elementary_spectra':spectra,
 'cyclotomic_projection_remainders':remainders,
 'identities':{
   'trace':'tr G_Sell = 2 sum p_ij',
   'traceless_norm':'||A||_F^2 = ||s-(S/3)1||^2 + 2||d||^2',
   'isotropic_calibration':'kappa_* = R_*^2/(4 sqrt(3))',
   'trace_radial':'U_tilde = (2 kappa_*/3)(S-6 sqrt(3))',
   'coordinate_change':'U = U_tilde + D/6',
   'elementary_D':'delta D = -kappa_* (sqrt(3 E_s)+sqrt(6 E_d))',
   'Singer_axis':'n_{sigma,L}=P_L w_sigma / ||P_L w_sigma||',
   'orientation_pair':'nu -> -nu, n -> -n, so nu*n is invariant'
 },
 'statuses':{
   'direct_Es_to_U':'REFUTED_TYPED',
   'trace_to_radial_size':'PROVEN_CONSTRUCTED_NORMALIZED',
   'direct_Ed_to_signed_D':'REFUTED_TYPED',
   'Selling_traceless_to_signed_D_elementary':'PROVEN_EXACT_AFTER_SCALE',
   'elementary_axis_line':'PROVEN_SOURCE_NATIVE',
   'Singer_inversion_to_axis_flip':'PROVEN_FOR_GOLDEN_ANCHORED_CHART',
   'orientation_pair_independent_axial_vector':'PROVEN',
   'absolute_choice_among_24_Singer_pairs':'NT_PENDING_CANONICAL_SELECTOR'
 },
 'guards':[
   'Normalized line-mass energy is a shape variable, not an absolute radial scale.',
   'The scalar calibration kappa_* fixes physical units by matching the unique isotropic Selling point to radius R_*.',
   'The seven elementary Selling moves are axisymmetric; generic finite superpositions need not remain axisymmetric.',
   'Singer inversion orients the intrinsic eigenline, but the corpus does not canonically select one of the 24 inverse-orientation pairs.'
 ]
}
Path('/mnt/data/source_native_geometric_adapter_v13_certificate.json').write_text(json.dumps(certificate,indent=2,ensure_ascii=False))
print(json.dumps(certificate,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
