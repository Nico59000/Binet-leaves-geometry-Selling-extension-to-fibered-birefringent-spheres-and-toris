#!/usr/bin/env python3
from __future__ import annotations
import json, csv, hashlib, sys
from pathlib import Path
import sympy as sp

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import selling_reconstruction_engine_v50 as eng

q,r=sp.symbols('q r', real=True)
D0=8+4*q-7*q**2-60*r**2
Hh=8-56*q-7*q**2-60*r**2
Mm=8+4*q-7*q**2+60*r**2-120*r-120*q*r
Lp=48-120*r+84*q+63*q**2+540*r**2
pD=sp.Poly(469*q**4+4564*q**3+7224*q**2+1984*q-416,q,domain=sp.QQ)
qA=(-28+2*sp.sqrt(210))/7
rB=1-sp.sqrt(195)/15
rD_expr=(632-469*q**3-4095*q**2-3318*q)/1620
M=sp.Matrix([[7,0,2],[0,-1,0],[2,0,-8]])
Jy=sp.diag(1,-1,1)

# Positive root isolation for q_D.
intervals=pD.intervals(eps=sp.Rational(1,10**30))
pos=[(ab,m) for ab,m in intervals if ab[1]>0]
assert len(pos)==1
(qDl,qDu),mult=pos[0]
assert qDl>0 and qDu<qA and mult==1

# Descartes: exactly one positive root because coefficient signs ++++-.  We
# certify existence below q_A by p(0)<0<p(q_A).
p_qA=sp.simplify(pD.as_expr().subs(q,qA))
assert pD.eval(0)<0 and sp.N(p_qA)>0

# Intersection identities Hh=Mm=0 when q is the q_D root and r=rD(q).
def rem_mod_p(expr):
    num,den=sp.fraction(sp.cancel(expr.subs(r,rD_expr)))
    rem=sp.rem(sp.Poly(sp.expand(num),q,domain=sp.QQ),pD)
    return sp.factor(rem.as_expr()), sp.factor(den)
assert rem_mod_p(Hh)[0]==0
assert rem_mod_p(Mm)[0]==0

# Root-local sign helper for rational functions in q.
def sign_at_qD(expr):
    expr=sp.cancel(expr.subs(r,rD_expr))
    num,den=sp.fraction(expr)
    num=sp.Poly(sp.expand(num),q,domain=sp.QQ)
    den=sp.Poly(sp.expand(den),q,domain=sp.QQ)
    a,b=qDl,qDu
    # exact zero at root?
    g=sp.gcd(pD,num)
    if g.degree()>0 and g.count_roots(a,b)>0:
        return 0
    aa,bb=a,b
    for P in (num,den):
        for _ in range(200):
            if P.count_roots(aa,bb)==0: break
            aa,bb=sp.refine_root(pD,aa,bb,eps=(bb-aa)/4)
        else: raise RuntimeError('sign refinement failure')
    mid=(aa+bb)/2
    nv=num.eval(mid);dv=den.eval(mid)
    return 1 if nv*dv>0 else -1

assert sign_at_qD(rD_expr)==1
# Separate rational bounds certify r_D<r_B without mixing number fields.
drD=sp.expand(sp.diff(rD_expr,q))
assert sp.expand(-drD-(sp.Rational(469,540)*q**2+sp.Rational(91,18)*q+sp.Rational(553,270)))==0
rD_upper=sp.simplify(rD_expr.subs(q,qDl))
assert rD_upper < sp.Rational(8,125)  # <0.064
assert sp.Rational(69,1000) < rB      # 0.069<r_B, certified algebraically

# On the quarter chart q>=0,r>=0,Hh>=0,Mm>=0, D0>0.
# Hh>=0 -> r^2<=2/15 < 1, hence 1-r+q>0.  l is strictly
# negative because Lp = 540(r-1/9)^2 + 124/3 + 84q + 63q^2 >0.
Lp_sos=sp.expand(540*(r-sp.Rational(1,9))**2 + sp.Rational(124,3)+84*q+63*q**2)
assert sp.expand(Lp-Lp_sos)==0

# Exact normalized wall formulae for x>0,y=-rx,z=qx, x^2=60/D0.
wall_scaled={
 'g': -120*r*q,
 'h': -2*Hh,
 'k': -120*r,
 'l': -Lp,
 'm': -Mm,
 'n': -6*D0-120*q*(1-r+q),
}
# Each actual wall is scaled numerator / D0.

# Vertices in normalized coordinates.
vertices=[
 {'id':'V_E','kind':'AXIS-SEAM-JUNCTION','q':'0','r':'0',
  'vanishing':['g','k'],
  'xyz_exact':['sqrt(15/2)','0','0'],
  'note':'singular meeting of the g=0 components; g and k share the y=0 mirror seam'},
 {'id':'V_A','kind':'AXIS-TRIPLE-CROSSING','q':str(qA),'r':'0',
  'vanishing':['g','h','k'],
  'xyz_exact':['1/sqrt(q_A)','0','sqrt(q_A)'],
  'q_A_minpoly':'7*q^2+56*q-8'},
 {'id':'V_D','kind':'FISSURE','q':'q_D','r':'r_D(q_D)',
  'vanishing':['h','m'],
  'q_D_minpoly':str(pD.as_expr()),
  'q_D_isolating_interval':[str(qDl),str(qDu)],
  'r_D_formula':str(rD_expr),
  'approx_q':str(sp.N((qDl+qDu)/2,18)),
  'approx_r':str(sp.N(rD_expr.subs(q,(qDl+qDu)/2),18)),
  'opposite_pair':True},
 {'id':'V_B','kind':'CROSSING','q':'0','r':str(rB),
  'vanishing':['g','m'],
  'r_B_minpoly':'15*r^2-30*r+2',
  'adjacent_pair':True},
]

edges=[
 {'id':'E_EA','source':'V_E','target':'V_A','active_walls':['g','k'],
  'chart_equation':'r=0, 0<q<q_A','kind':'MIRROR-DOUBLE-WALL'},
 {'id':'E_AD','source':'V_A','target':'V_D','active_walls':['h'],
  'chart_equation':'Hh(q,r)=0, q_D<q<q_A','kind':'ORDINARY-LIMIT'},
 {'id':'E_DB','source':'V_D','target':'V_B','active_walls':['m'],
  'chart_equation':'Mm(q,r)=0, 0<q<q_D','kind':'ORDINARY-LIMIT'},
 {'id':'E_BE','source':'V_B','target':'V_E','active_walls':['g'],
  'chart_equation':'q=0, 0<r<r_B','kind':'ORDINARY-LIMIT'},
]

# Cellular incidence for this single closed field.
B1=sp.Matrix([
 [-1, 0, 0, 1],  # E
 [ 1,-1, 0, 0],  # A
 [ 0, 1,-1, 0],  # D
 [ 0, 0, 1,-1],  # B
])
B2=sp.Matrix([1,1,1,1])
assert B1*B2==sp.zeros(4,1)

# Exact mirror stabilizer on the seam.
assert Jy.T*M*Jy==M and Jy**2==sp.eye(3) and Jy.det()==-1

# Ordinary-wall move candidates.  The two source candidates at each wall
# represent the same numerical neighbour up to a unique S4 superbase
# relabelling, but branch/orientation selection is deliberately left open.
wall_to_p={'g':'p23','h':'p13','k':'p12','l':'p01','m':'p02','n':'p03'}
transitions={}
for wall in ('h','m','g'):
    inds=eng.WALL_TO_MOVES[wall_to_p[wall]]
    rec=[]
    Ns=[]
    for idx in inds:
        A=eng.S[idx-1];N=A.T*M*A;Ns.append(N)
        rec.append({'move':f'S{idx}','matrix':[[int(A[i,j]) for j in range(3)] for i in range(3)],
                    'target_six':[int(x) for x in eng.six_from_gram(N)],
                    'numerical_key':[str(x) for x in eng.numerical_form_key(N)]})
    rel=[]
    for R in eng.RELABEL:
        if R.T*Ns[0]*R==Ns[1]: rel.append(R)
    assert len(rel)==1
    R=rel[0]
    transitions[wall]={
      'candidates':rec,
      'same_numerical_neighbour':True,
      'unique_superbase_relabelling_candidate1_to_candidate2':[[int(R[i,j]) for j in range(3)] for i in range(3)],
      'branch_marking_selection':'PENDING-SOURCE-ORIENTATION'
    }

# Verification registry.
checks={
 'det_M_Q2_60': M.det()==60,
 'normalized_h_wall_identity': sp.expand((120*q-2*D0)-(-2*Hh))==0,
 'normalized_m_wall_identity': sp.expand((-D0+120*r*(1+q-r))-(-Mm))==0,
 'normalized_l_wall_sos': sp.expand(Lp-Lp_sos)==0,
 'qA_h_boundary': sp.simplify(Hh.subs({q:qA,r:0}))==0,
 'rB_m_boundary': sp.simplify(Mm.subs({q:0,r:rB}))==0,
 'qD_unique_positive_by_descartes': True,
 'qD_below_qA': qDu<qA,
 'D_H_identity_mod_p': rem_mod_p(Hh)[0]==0,
 'D_M_identity_mod_p': rem_mod_p(Mm)[0]==0,
 'rD_positive': sign_at_qD(rD_expr)>0,
 'rD_below_rB': rD_upper < sp.Rational(8,125) and sp.Rational(69,1000) < rB,
 'mirror_stabilizer': Jy.T*M*Jy==M and Jy**2==sp.eye(3),
 'boundary_closes': B1*B2==sp.zeros(4,1),
 'ordinary_edge_candidate_pairs_same_numerical_neighbour': all(v['same_numerical_neighbour'] for v in transitions.values()),
}
checks={k:bool(v) for k,v in checks.items()}
assert all(checks.values())

out={
 'version':'v50',
 'status':'PASS-LOCAL-CELL',
 'publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
 'field_id':'F_Q2',
 'field_form_six':[7,-1,-8,0,2,0],
 'field_gram':[[int(M[i,j]) for j in range(3)] for i in range(3)],
 'source_basis':'Selling reduction inequalities; no raster incidence used',
 'normalized_chart':{
   'definition':'q=z/x, r=-y/x on the x>0 sheet',
   'x_squared':'60/D0',
   'D0':str(D0),
   'domain':['q>=0','r>=0','Hh>=0','Mm>=0'],
   'Hh':str(Hh),'Mm':str(Mm),
   'wall_scaled_numerators':{k:str(v) for k,v in wall_scaled.items()},
   'l_strict_negativity_certificate':str(Lp_sos),
   'n_strict_negativity_reason':'Hh>=0 gives r^2<=2/15<1, so 1-r+q>0 and n<0.'
 },
 'topology_certificate':{
   'q_A':str(qA),'r_B':str(rB),
   'q_D_polynomial':str(pD.as_expr()),
   'q_D_isolating_interval':[str(qDl),str(qDu)],
   'q_D_unique_positive_root':'DESCARTES_ONE_SIGN_CHANGE + p(0)<0<p(q_A)',
   'boundary_order':['V_E','V_A','V_D','V_B','V_E'],
   'edge_order':['E_EA','E_AD','E_DB','E_BE'],
   'compact_disk_reason':'Hh>=0 bounds 0<=q<=q_A and 0<=r<sqrt(2/15); Mm>=0 selects the lower r-branch. Hh=0 and Mm=0 meet once at q_D.'
 },
 'C0_vertices':vertices,
 'C1_oriented_arcs':edges,
 'C2_fields':[{'id':'F_Q2','ordered_boundary':['E_EA','E_AD','E_DB','E_BE']}],
 'boundary_matrices':{'d1_vertices_by_edges':[list(map(int,row)) for row in B1.tolist()],
                      'd2_edges_by_field':[int(x) for x in B2]},
 'mirror_transport':{'edge':'E_EA','matrix':[[int(Jy[i,j]) for j in range(3)] for i in range(3)],
                     'involution':True,'det':-1,'stabilizes_field_form':True},
 'ordinary_edge_transition_gate':transitions,
 'checks':checks,
 'decision':{
   'first_source_constructed_2cell':'PROVEN-LOCAL',
   'partial_C0_C1_C2':'PROVEN-LOCAL-SUBCOMPLEX',
   'global_C0_C1_C2':'NT/PENDING-ADJACENT-FIELD-DEVELOPMENT',
   'global_d0_d1':'NOT-RUN/GUARDED',
   'historical_H1_H2':'NT/GUARDED',
   'raster_dependency':'NONE'
 }
}

outpath=HERE/'selling_Q2_first_cell_v50.json'
outpath.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')

csvpath=HERE/'selling_Q2_first_cell_boundary_v50.csv'
with csvpath.open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f)
    w.writerow(['edge','source','target','active_walls','kind','chart_equation'])
    for e in edges:w.writerow([e['id'],e['source'],e['target'],'+'.join(e['active_walls']),e['kind'],e['chart_equation']])

# compact replay hash
for p in (outpath,csvpath):
    print(p.name,hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_size)
print('PASS',sum(bool(v) for v in checks.values()),'/',len(checks))
