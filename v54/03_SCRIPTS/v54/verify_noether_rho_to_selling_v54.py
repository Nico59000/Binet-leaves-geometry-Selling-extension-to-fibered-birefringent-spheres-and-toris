#!/usr/bin/env python3
import sympy as sp, math, itertools, json, hashlib
from pathlib import Path
ROOT=Path('/mnt/data/v54_gate2/noether_rho')
u,v,w=sp.symbols('u v w'); x,y,z=sp.symbols('x y z')
x1,y1,z1,x2,y2,z2,x3,y3,z3=sp.symbols('x1 y1 z1 x2 y2 z2 x3 y3 z3')
A12=sp.Matrix([[u,x1,x2],[v,y1,y2],[w,z1,z2]])
A23=sp.Matrix([[u,x2,x3],[v,y2,y3],[w,z2,z3]])
A31=sp.Matrix([[u,x3,x1],[v,y3,y1],[w,z3,z1]])
expr=sp.expand(A12.det()**2*A23.det()**2*A31.det()**2/sp.Integer(6))
P=sp.Poly(expr,u,v,w,x1,y1,z1,x2,y2,z2,x3,y3,z3,domain=sp.QQ)
basis_exps=[(4,0,0),(3,1,0),(2,2,0),(1,3,0),(0,4,0),(3,0,1),(2,1,1),(1,2,1),(0,3,1),(2,0,2),(1,1,2),(0,2,2),(2,0,2),(1,1,2),(0,2,2)]
# Correct canonical degree-4 basis in three variables
basis_exps=[(4,0,0),(3,1,0),(2,2,0),(1,3,0),(0,4,0),(3,0,1),(2,1,1),(1,2,1),(0,3,1),(2,0,2),(1,1,2),(0,2,2),(1,0,3),(0,1,3),(0,0,4)]
multinom=[math.factorial(4)//(math.factorial(a)*math.factorial(b)*math.factorial(c)) for a,b,c in basis_exps]

def rho_from_actual(actual):
    coeff={e:sp.Rational(actual.get(e,0),m) for e,m in zip(basis_exps,multinom)}
    psi=0
    for mon,k in P.terms():
        eu,ev,ew=mon[:3]; e1=mon[3:6]; e2=mon[6:9]; e3=mon[9:12]
        c1,c2,c3=coeff.get(e1,0),coeff.get(e2,0),coeff.get(e3,0)
        if c1 and c2 and c3:
            psi += k*c1*c2*c3*u**eu*v**ev*w**ew
    psi=sp.expand(psi)
    D=0
    for e,a in actual.items():
        if a: D += sp.Rational(a)*sp.diff(psi,u,e[0],v,e[1],w,e[2])
    return sp.factor(psi), sp.factor(sp.expand(D/sp.Integer(144)))

def actual_dict(poly):
    return {mon:coef for mon,coef in sp.Poly(sp.expand(poly),x,y,z,domain=sp.QQ).terms() if sum(mon)==4}

def qmat(poly):
    p=sp.Poly(sp.expand(poly),u,v,w,domain=sp.QQ); vs=[u,v,w]; C=sp.zeros(3)
    for i,t in enumerate(vs):C[i,i]=p.coeff_monomial(t**2)
    for i in range(3):
        for j in range(i+1,3):C[i,j]=C[j,i]=p.coeff_monomial(vs[i]*vs[j])/2
    return C

def quartic_tensor(actual):
    T=sp.MutableDenseNDimArray.zeros(3,3,3,3)
    for exp,a in actual.items():
        m=math.factorial(4)//math.prod(math.factorial(e) for e in exp)
        val=sp.Rational(a,m); inds=[0]*exp[0]+[1]*exp[1]+[2]*exp[2]
        for perm in set(itertools.permutations(inds,4)):T[perm]=val
    return T

def contract(actual,rho):
    T=quartic_tensor(actual);C=qmat(rho);Q=sp.zeros(3)
    for i in range(3):
        for j in range(3):Q[i,j]=sp.factor(sum(T[i,j,k,l]*C[k,l] for k in range(3) for l in range(3)))
    return Q,C

def js(expr): return str(sp.factor(expr))
def jmat(A): return [[js(A[i,j]) for j in range(A.cols)] for i in range(A.rows)]
def primitive_integer_matrix(Q):
    den=sp.ilcm(*[sp.denom(q) for q in Q])
    A=Q*den; nums=[abs(int(q)) for q in A if q]
    g=0
    for n in nums:g=math.gcd(g,n)
    if g:A=A/g
    # global sign chosen so signature has one positive/two negative when possible
    return sp.Matrix(A),int(den),int(g or 1)

def inertia_symmetric(A):
    # exact charpoly roots with multiplicities, sign isolated numerically only after algebraic root construction
    cp=sp.Poly(A.charpoly().as_expr(),sp.Symbol('lambda'))
    roots=sp.polys.polytools.all_roots(cp)
    pos=neg=zero=0
    for r in roots:
        if r.is_positive:pos+=1
        elif r.is_negative:neg+=1
        elif r.is_zero:zero+=1
        else:
            rv=sp.N(r,80)
            if rv>0:pos+=1
            elif rv<0:neg+=1
            else:zero+=1
    return [pos,neg,zero]

# Historical first printed Gordan/Noether automorphic quartic.
fG=x**3*y+y**3*z+z**3*x
aG=actual_dict(fG); psiG,rhoG=rho_from_actual(aG);QG,CG=contract(aG,rhoG)
assert rhoG==0 and QG==sp.zeros(3)
# Calibration specialization inside the same general ternary quartic law, not claimed printed by Noether.
fS=x**2*y**2-x**2*z**2-y**2*z**2
aS=actual_dict(fS);psiS,rhoS=rho_from_actual(aS);QS,CS=contract(aS,rhoS)
assert sp.factor(rhoS + (u**2+v**2-w**2)/648)==0
Qprim,den,gcd=primitive_integer_matrix(QS)
assert Qprim==sp.diag(-1,-1,1)
assert Qprim.det()==1 and inertia_symmetric(Qprim)==[1,2,0]
Gs=[sp.Matrix([[1,1,0],[0,1,0],[0,0,1]]),sp.Matrix([[1,0,0],[1,1,0],[0,0,1]]),sp.Matrix([[0,-1,0],[1,0,0],[0,0,1]])]
nat=[]
for G in Gs:
    Gi=G.inv();vec=Gi*sp.Matrix([x,y,z])
    fT=sp.expand(fS.subs({x:vec[0],y:vec[1],z:vec[2]},simultaneous=True));aT=actual_dict(fT)
    _,rhoT=rho_from_actual(aT);CT=qmat(rhoT); QT,_=contract(aT,rhoT)
    c_ok=sp.simplify(CT-G*CS*G.T)==sp.zeros(3)
    q_ok=sp.simplify(QT-G.inv().T*QS*G.inv())==sp.zeros(3)
    assert c_ok and q_ok
    nat.append({'G':[[int(e) for e in row] for row in G.tolist()],'rho_contravariance':c_ok,'Q_covariance':q_ok})

cert={
 'version':'v54','status':'PASS-EXACT-NOETHER-RHO-SPECIALIZATION-AND-GL3-SQUARE',
 'algorithm':{
   'psi':'Salmon/Ohno determinant construction det(A12)^2 det(A23)^2 det(A31)^2 / 6 with divided-power quartic coefficients',
   'rho':'D_f(psi)/144',
   'bridge':'Phi(f,rho)_{ij}=F_{ijkl} rho^{kl}',
   'psi_expanded_template_terms':len(P.terms())
 },
 'historical_printed_instance':{
   'source':'Noether 1908 introduction, Gordan special automorphic quartic',
   'f':'x^3*y + y^3*z + z^3*x','psi':js(psiG),'rho':js(rhoG),'Q':jmat(QG),
   'selling_gate':'STOP/RHO-ZERO/DEGENERATE-Q',
   'status':'PROVEN-HISTORICAL-INSTANCE-NONLANDING',
   'guard':'This is an instance result, not a refutation of the general GL3-equivariant bridge.'
 },
 'calibration_specialization':{
   'historical_printed':False,'f':'x^2*y^2 - x^2*z^2 - y^2*z^2','psi':js(psiS),'rho':js(rhoS),
   'rho_matrix':jmat(CS),'Q_rational':jmat(QS),'Q_primitive_integer':[[int(x) for x in row] for row in Qprim.tolist()],
   'primitive_det':int(Qprim.det()),'primitive_signature':inertia_symmetric(Qprim),
   'status':'PROVEN-SOURCE-LAW-CALIBRATION-LANDS-IN-SELLING-LOCUS'
 },
 'GL3_exact_naturality':{'cases':nat,'count':len(nat),'status':'PASS'},
 'guards':[
   'The calibration quartic is not promoted as a printed historical Noether example.',
   'Primitive integral rescaling is an additional normalization after Phi; it preserves the real signature but is not GL3-linear as a family normalization.',
   'The printed Gordan example is fail-closed before Selling because rho=0.'
 ]
}
(ROOT/'noether_rho_historical_and_calibration_v54.json').write_text(json.dumps(cert,indent=2,sort_keys=True)+"\n")
# Selling seed for primitive calibration Q
seed={
 'schema':'selling.seed.v1','name':'v54_noether_rho_calibration_primitive_Q',
 'provenance':{'kind':'Noether-rho-source-law-calibration','historical_printed':False,'certificate':'noether_rho_historical_and_calibration_v54.json'},
 'indefinite_gram':[[int(e) for e in row] for row in Qprim.tolist()],
 'chart_transform':[[1,0,0],[0,1,0],[0,0,1]],
 'quotient':{'mode':'superbase','deck_generators':[]},
 'compile':{'classifier':'exact','max_fields':8,'deck_group_max_size':256,'use_superbase_relabelling':True,'projective_sign':True,'fallback_sampled_classifier':False},
 'render':{'q_range':[-2.0,2.0],'r_range':[-1.2,1.2],'resolution':650,'labels':True,'fill_fields':False,'show_domain_boundary':True,'line_width':0.7}
}
(ROOT/'noether_rho_calibration_v54.seed.json').write_text(json.dumps(seed,indent=2)+"\n")
print(json.dumps({'status':cert['status'],'rho_historical':str(rhoG),'rho_calibration':str(rhoS),'Qprim':[[int(e) for e in row] for row in Qprim.tolist()],'naturality':len(nat)},indent=2))
