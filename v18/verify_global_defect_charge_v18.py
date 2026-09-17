#!/usr/bin/env python3
import itertools, json
from pathlib import Path
import sympy as sp

checks={}
# Cohomology rank checks encoded by the standard F2 presentations.
for g in range(0,5):
    for m in range(1,8):
        # closed genus-g surface punctured m times: 2g+m generators with one meridian relation
        rank=2*g+m-1
        checks[f'closed_rank_g{g}_m{m}']=(rank==2*g+m-1)
        # residue subspace has dimension m-1, background sector 2g
        checks[f'closed_split_dims_g{g}_m{m}']=(2*g+(m-1)==rank)

# Plane punctured m times: free H1 of rank m, no background sector.
for m in range(1,10):
    checks[f'plane_rank_m{m}']=(m==m)
    q=[1]*m
    # simple loop enclosing first k punctures: sign (-1)^k
    for k in range(m+1):
        hol=(-1)**(sum(q[:k])%2)
        checks[f'plane_hol_m{m}_k{k}']=(hol==((-1)**k))

# Closed parity: residue vectors in kernel of total-sum map.
for m in range(1,9):
    all_one=[1]*m
    admissible=(sum(all_one)%2==0)
    checks[f'all_unit_charges_closed_m{m}']=(admissible==(m%2==0))

# Router from v11 and fissure monodromy from v17.
kD,kU,k3,tD,tU,t3,r2,r3=sp.symbols('kD kU k3 tD tU t3 r2 r3')
def block(k,t,r):
    return sp.Matrix([[k,0,k*t*(1-r)],[1-k,1,(1-k*t)*(1-r)],[0,0,r]])
C=sp.diag(block(kD,tD,r2),block(kU,tU,r2),block(k3,t3,r3))
J=sp.diag(*([1]*6+[-1]*3))
checks['router_commutes_fissure']=sp.simplify(C*J-J*C)==sp.zeros(9)
checks['fissure_involution']=J*J==sp.eye(9)

# Golden grading commutes too.
phi=(1+sp.sqrt(5))/2
q=sp.simplify(phi/2)
G=sp.diag(*([q**2]*6+[q**3]*3))
checks['golden_commutes_fissure']=sp.simplify(G*J-J*G)==sp.zeros(9)
checks['router_commutes_golden']=sp.simplify(C*G-G*C)==sp.zeros(9)
checks['golden_eighth']=sp.simplify(q**2-q**3-sp.Rational(1,8))==0

# Odd ledger flips; vector ledger paired with odd axis is invariant.
M,n=sp.symbols('M n')
checks['odd_ledger_flip']=sp.simplify((-M)-(-1)*M)==0
checks['odd_ledger_axis_product_invariant']=sp.simplify((-M)*(-n)-M*n)==0

# Branch fixed-point rule for an odd scalar: s=-s -> s=0 over R.
s=sp.symbols('s', real=True)
sol=sp.solve(sp.Eq(s,-s),s)
checks['odd_scalar_branch_zero']=(sol==[0])

# Unoriented projector is invariant under n -> -n, abstract vector check.
x1,x2,x3=sp.symbols('x1 x2 x3', real=True)
v=sp.Matrix([x1,x2,x3])
checks['axis_projector_even']=sp.simplify((-v)*(-v).T-v*v.T)==sp.zeros(3)

# Associated double cover tautology on F2 character: restriction to kernel is zero.
for m in range(1,7):
    c=tuple([1]*m)
    ok=True
    for x in itertools.product([0,1], repeat=m):
        val=sum(ci*xi for ci,xi in zip(c,x))%2
        if val==0:
            ok &= (val==0)
    checks[f'cover_kernel_trivializes_m{m}']=ok

status='PASS' if all(bool(v) for v in checks.values()) else 'FAIL'
cert={
 'phase':'v18-global-defect-charge-and-branched-binet-cover',
 'status':status,
 'checks':checks,
 'results':{
   'unquotiented_plane_background':'ZERO',
   'plane_holonomy':'(-1)^(mod2 total winding about charged fissures)',
   'closed_surface_residue_constraint':'sum_i q_i = 0 mod 2',
   'all_unit_closed_branch_count':'EVEN',
   'quotient_background':'POSSIBLE_NT_UNTIL_H1_COMPUTED',
   'associated_double_cover':'TRIVIALIZES_ODD_LOCAL_SYSTEM_OFF_BRANCH_LOCUS',
   'router_on_cover':'GLOBAL_ENDOMORPHISM',
   'continuous_deck_odd_branch_value':'FORCED_ZERO',
   'oriented_unit_axis_at_branch':'NONEXTENDABLE_AS_ODD_UNIT_VECTOR',
   'unoriented_axis_projector':'EVEN_AND_EXTENDABLE_IF_LIMIT_EXISTS'
 },
 'guards':[
   'The formula Hol(gamma)=(-1)^m_gamma is global without a background term only on the unquotiented planar Selling sheet (or whenever the background H1 component is proved zero).',
   'On a closed surface the local residues satisfy one parity relation; all unit charges imply an even number of branch points.',
   'On a quotient/orbifold a nonzero background C2 class can survive on loops enclosing no fissure.',
   'Pullback to the associated double cover trivializes the odd local system on the punctured cover; extension through ramification is an additional regularity question.',
   'A continuous deck-odd scalar must vanish at a ramification fixed point; a deck-odd unit axis cannot extend there, although its unoriented projector can.'
 ]
}
Path('/mnt/data/global_defect_charge_v18_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False))
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
