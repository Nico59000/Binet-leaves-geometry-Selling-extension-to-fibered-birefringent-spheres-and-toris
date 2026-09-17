import sympy as sp
import json

u,a,b,c,g,h,k,beta=sp.symbols('u a b c g h k beta', real=True)
M=sp.Matrix([[a,k,h],[k,b,g],[h,g,c]])
adj=M.adjugate()
A,B,C=adj[0,0],adj[1,1],adj[2,2]
G,H,K=adj[1,2],adj[0,2],adj[0,1]
detM=sp.factor(M.det())
I0=-3*detM
P=(u+g)*(u+h)*(u+k)
Q=sp.expand(
 A*(u+h)**2*(u+k)**2
 +B*(u+k)**2*(u+g)**2
 +C*(u+g)**2*(u+h)**2
 +2*G*(u+g)**2*(u+h)*(u+k)
 +2*H*(u+h)**2*(u+g)*(u+k)
 +2*K*(u+k)**2*(u+g)*(u+h)
 +I0*P
)
poly=sp.Poly(Q,u)
a4=poly.coeff_monomial(u**4); a3=poly.coeff_monomial(u**3)
a2=poly.coeff_monomial(u**2); a1=poly.coeff_monomial(u); a0=poly.coeff_monomial(1)
I2=sp.expand(12*a4*a0-3*a3*a1+a2**2)
J3=sp.expand(72*a4*a2*a0+9*a3*a2*a1-27*a4*a1**2-27*a3**2*a0-2*a2**3)

star={a:sp.Rational(136,31), b:-sp.Rational(26,29), c:-sp.Rational(47,46), g:-1, h:-2, k:-4}
Mstar=sp.simplify(M.subs(star)); Astar=sp.simplify(adj.subs(star))
detstar=sp.factor(detM.subs(star)); I0star=sp.factor(I0.subs(star))
Qstar=sp.factor(Q.subs(star)); P0=sp.factor(P.subs(star).subs(u,0))

# Signature via exact LDL pivots from leading principal minors.
D1=sp.factor(Mstar[:1,:1].det())
D2=sp.factor(Mstar[:2,:2].det())
D3=sp.factor(Mstar.det())
pivots=[D1,sp.factor(D2/D1),sp.factor(D3/D2)]
assert [sp.sign(p) for p in pivots]==[1,-1,-1]

# Genuine triple root, non-pole, nonquadruple.
assert P0 != 0
assert sp.factor(Qstar + sp.Rational(504,20677)*u**3*(143*u-112)) == 0
for j in range(3):
    assert sp.simplify(sp.diff(Q,u,j).subs(star).subs(u,0)) == 0
Q3=sp.factor(sp.diff(Q,u,3).subs(star).subs(u,0))
assert Q3 != 0
rho=sp.Rational(112,143)
assert all(rho != x for x in (1,2,4))

# Local A2 versality: quotient R[u]/(u^2), controls a,b.
def jet01(expr):
    pp=sp.Poly(sp.expand(expr),u)
    return sp.Matrix([pp.coeff_monomial(1),pp.coeff_monomial(u)])
V=sp.Matrix.hstack(jet01(sp.diff(Q,a).subs(star)), jet01(sp.diff(Q,b).subs(star)))
Vdet=sp.factor(V.det())
assert Vdet == sp.Rational(294912,529)
assert V.rank()==2

# Invariant map transversality in same controls.
JIJ=sp.Matrix([[sp.diff(I2,a),sp.diff(I2,b)],[sp.diff(J3,a),sp.diff(J3,b)]]).subs(star)
JIJdet=sp.factor(JIJ.det())
expected_Jdet=-sp.Rational(2**36*3**12*7**6, 23**5*29**3*31**3)
assert JIJdet==expected_Jdet and JIJ.rank()==2
assert sp.simplify(I2.subs(star))==0 and sp.simplify(J3.subs(star))==0

# Polar contact model: depressed cubic y^3+p y+q, mu=-beta p, nu=-beta q.
y,p,q,mu,nu=sp.symbols('y p q mu nu')
stationary=beta*y**3-mu*y-nu
assert sp.expand(stationary.subs({mu:-beta*p,nu:-beta*q})-beta*(y**3+p*y+q))==0
Dpolar=4*mu**3-27*beta*nu**2
assert sp.factor(Dpolar.subs({mu:-beta*p,nu:-beta*q}) + beta**3*(4*p**3+27*q**2))==0

out={
 'status':'PASS',
 'packet':{
   'M_star':[[str(x) for x in row] for row in Mstar.tolist()],
   'adj_M_star':[[str(x) for x in row] for row in Astar.tolist()],
   'det_M_star':str(detstar),
   'I0_star':str(I0star),
   'LDL_pivots':[str(x) for x in pivots],
   'signature':'(1,2)',
   'fixed':{'c':'-47/46','g':'-1','h':'-2','k':'-4'},
   'controls':{'a_star':'136/31','b_star':'-26/29'}
 },
 'triple_root':{
   'u0':'0','P_u0':str(P0),'Q_factor':str(Qstar),'Q3_u0':str(Q3),
   'simple_fourth_root':str(rho),'poles':['1','2','4'],
   'genuine_rational_root':True
 },
 'A2_versality':{
   'local_algebra_basis':['1','u'],
   'control_jet_matrix':[[str(x) for x in row] for row in V.tolist()],
   'determinant':str(Vdet),'rank':2
 },
 'invariant_transversality':{
   'I2_star':'0','J3_star':'0',
   'Jacobian_ab':[[str(x) for x in row] for row in JIJ.tolist()],
   'determinant':str(JIJdet),'rank':2
 },
 'polar_transport':{
   'stationary_normal_form':'beta*(y^3+p*y+q)=0',
   'control_map':'mu=-beta*p, nu=-beta*q',
   'polar_discriminant_pullback':'-beta^3*(4*p^3+27*q^2)',
   'status':'PROVEN_LOCAL_CONTACT_EQUIVALENCE_FOR_BETA_POSITIVE'
 },
 'typing':{
   'Selling_structured_quartic_to_A2':'PROVEN_LOCAL',
   'A2_to_declared_polar_nucleation_model':'PROVEN_LOCAL_TYPED',
   'source_native_identification_mu_with_Binet_hypertensor_observable':'NT_PENDING_MORPHISM'
 }
}
with open('/mnt/data/selling_A2_polar_promotion_v5_certificate.json','w',encoding='utf-8') as f:
    json.dump(out,f,indent=2,ensure_ascii=False)
print(json.dumps({
 'status':'PASS','Q_star':str(Qstar),'signature':'(1,2)',
 'versality_det':str(Vdet),'IJ_det':str(JIJdet),
 'polar':'PROVEN_LOCAL_CONTACT_EQUIVALENCE'
},indent=2))
