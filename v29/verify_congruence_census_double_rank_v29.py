#!/usr/bin/env python3
import itertools, json, math, cmath
from pathlib import Path
import numpy as np
import sympy as sp

OUT=Path("/mnt/data")
checks={}

S_int=[
[[-1,0,1],[0,1,0],[0,0,1]],
[[-1,1,0],[0,1,0],[0,0,1]],
[[-1,1,1],[0,1,0],[0,0,1]],
[[0,-1,0],[-1,0,0],[-1,-1,1]],
[[0,0,-1],[-1,1,-1],[-1,0,0]],
[[1,-1,-1],[0,0,-1],[0,-1,0]],
[[1,0,0],[0,-1,1],[0,0,1]],
[[1,0,0],[0,1,0],[0,1,-1]],
[[1,0,0],[0,1,0],[1,0,-1]],
[[1,0,0],[0,1,0],[1,1,-1]],
[[1,0,0],[1,-1,0],[0,0,1]],
[[1,0,0],[1,-1,1],[0,0,1]],
]
S=[np.array(M,dtype=int) for M in S_int]
S2=[M%2 for M in S]
checks['all_integer_generators_involutions']=all(np.array_equal(M@M,np.eye(3,dtype=int)) for M in S)
checks['all_integer_generators_det_minus']=all(round(np.linalg.det(M))==-1 for M in S)

def rank_mod2(A):
    A=A.copy()%2; r=0
    for c in range(3):
        piv=next((i for i in range(r,3) if A[i,c]),None)
        if piv is not None:
            A[[r,piv]]=A[[piv,r]]
            for i in range(3):
                if i!=r and A[i,c]:
                    A[i]^=A[r]
            r+=1
    return r
def key(A): return tuple(int(x) for x in A.flatten())

G=[]
for bits in itertools.product([0,1], repeat=9):
    A=np.array(bits,dtype=np.uint8).reshape(3,3)
    if rank_mod2(A)==3: G.append(A)
G.sort(key=key); gid={key(A):i for i,A in enumerate(G)}
checks['GL32_order_168']=len(G)==168

I=np.eye(3,dtype=np.uint8)
seen={key(I)}; front=[I]
while front:
    A=front.pop()
    for B in S2:
        C=(A@B)%2; k=key(C)
        if k not in seen:
            seen.add(k); front.append(C)
checks['twelve_reductions_generate_GL32']=len(seen)==168
checks['triangle_S1S2S3']=np.array_equal((S2[0]@S2[1]@S2[2])%2,I)

meta={
1:("p12","u2"),7:("p12","u2"),
2:("p13","u1"),8:("p13","u1"),
3:("p01","u3"),6:("p01","u3"),
4:("p03","u2"),10:("p03","u2"),
5:("p02","u1"),12:("p02","u1"),
9:("p23","u3"),11:("p23","u3"),
}
checks['generator_line_partition_4_4_4']=sorted([sum(1 for j in meta if meta[j][1]==u) for u in ('u1','u2','u3')])==[4,4,4]

sigma=[0,4,1,5,2,6,3]
lines={"u1":[0,1,3],"u2":[0,4,5],"u3":[0,2,6]}
omega=cmath.exp(2j*math.pi/7)
Oms={}
for n,L0 in lines.items():
    L=set(L0)
    Z=sum((1 if sigma[k] in L else -1)*omega**k for k in range(7))
    Oms[n]=(Z**3).imag/(16*math.sqrt(2))
    checks[f'fourier_norm_{n}']=abs(abs(Z)**2-8)<1e-10

# exact cyclotomic distinctness
x=sp.symbols('x'); Phi=sum(x**k for k in range(7))
def rem7(expr):
    return sp.rem(sp.Poly(sp.expand(expr),x,domain=sp.QQ),
                  sp.Poly(Phi,x,domain=sp.QQ)).as_expr()
anti={}
for n,L0 in lines.items():
    L=set(L0)
    Z=sum((1 if sigma[k] in L else -1)*x**k for k in range(7))
    Z=rem7(Z); Zb=rem7(Z.subs(x,x**6))
    anti[n]=sp.expand(rem7(Z**3-Zb**3))
checks['Omega_u1_neq_u2_exact']=sp.expand(anti['u1']-anti['u2'])!=0
checks['Omega_u1_neq_u3_exact']=sp.expand(anti['u1']-anti['u3'])!=0
checks['Omega_u2_neq_u3_exact']=sp.expand(anti['u2']-anti['u3'])!=0

edges=[]
for u,A in enumerate(G):
    for j,B in enumerate(S2,start=1):
        v=gid[key((A@B)%2)]
        if u<v: edges.append((u,v,j,meta[j][1]))
checks['edge_count_1008']=len(edges)==1008

E=len(edges); V=len(G)
Dc=np.zeros((E,V)); D0=np.zeros((E,V)); s=np.zeros(E); chi=np.zeros(E)
for r,(u,v,j,L) in enumerate(edges):
    Dc[r,u]=Dc[r,v]=1
    D0[r,u]=-1; D0[r,v]=1
    s[r]=Oms[L]
rankDc=int(np.linalg.matrix_rank(Dc,tol=1e-10))
rankDca=int(np.linalg.matrix_rank(np.column_stack([Dc,s]),tol=1e-10))
rankD0=int(np.linalg.matrix_rank(D0,tol=1e-10))
rankD0a=int(np.linalg.matrix_rank(np.column_stack([D0,chi]),tol=1e-10))
checks['rank_Dc_168']=rankDc==168
checks['rank_Dc_aug_169']=rankDca==169
checks['rank_D0_167']=rankD0==167
checks['rank_D0_aug_167']=rankD0a==167

# Gauge reconciliation v15 -> all-minus v28/v29 lift gauge.
tau15=[-1,1,-1,-1,1,-1]; eps=[1,1,-1,-1,-1,1]
tau29=[eps[i]*tau15[i]*eps[(i+1)%6] for i in range(6)]
checks['crossing_tau_gauge_reconciliation']=tau29==[-1]*6

# Exact commuting square witnesses
pairs=[(1,2),(1,3),(2,3)]
charges=[]
for i,j in pairs:
    checks[f'commute_S{i}_S{j}_mod2']=np.array_equal((S2[i-1]@S2[j-1])%2,(S2[j-1]@S2[i-1])%2)
    Li,Lj=meta[i][1],meta[j][1]
    P=2*(Oms[Lj]-Oms[Li])
    charges.append({
        "word":[f"S{i}",f"S{j}",f"S{i}",f"S{j}"],
        "line_i":Li,"line_j":Lj,
        "period":P,
        "exact_formula":f"2*(Omega({Lj})-Omega({Li}))"
    })
checks['all_three_square_periods_nonzero']=all(abs(c['period'])>1e-8 for c in charges)

# Interpretation: normalized quotient history class is non-exact; flat continuous boost model is exact.
checks['history_nonexact_in_normalized_mod2_quotient']=rankDc!=rankDca
checks['krein_boost_zero_in_flat_normalized_model']=rankD0==rankD0a

status='PASS' if all(bool(v) for v in checks.values()) else 'FAIL'
cert={
 "phase":"v29-congruence-census-double-rank",
 "status":status,
 "checks":checks,
 "scope":{
   "primary_selling_global_census":"NT without a fixed-form chamber enumeration",
   "instantiated_quotient":"rho_2(Gamma_Sell^adj)=GL(3,2)",
   "normalization":"tau=-1 lift gauge, epsilon=1, kappa_hist=1, fixed Singer sigma_minus",
   "continuous_krein_component":"flat normalized B(0) model only; not a primary-Selling measured boost"
 },
 "census":{"vertices":168,"edges":1008,"generator_count":12},
 "ranks":{
   "Dc":rankDc,"Dc_aug_s":rankDca,
   "D0":rankD0,"D0_aug_chi":rankD0a
 },
 "decisions":{
   "history_class_in_instantiated_quotient":"NONZERO",
   "krein_continuous_class_in_flat_normalized_quotient":"ZERO",
   "historical_instance_history_class":"NT",
   "historical_instance_krein_boost":"NT"
 },
 "Omega":Oms,
 "first_real_memory_charges":charges,
 "guards":[
   "The 168/1008 census is an explicit congruence quotient of the source-anchored adjacent-superbase group, not the exhaustive chamber census of an arbitrary indefinite Selling form.",
   "The unit amplitude epsilon=1 and flat continuous Krein component B(0) are declared quotient normalizations.",
   "Non-exactness of s in this finite quotient is proved by gauge-invariant even commuting squares.",
   "No numerical claim about the actual historical amplitudes or continuous Krein boosts is inferred from the normalized quotient.",
   "The all-minus tau gauge is explicitly gauge-equivalent to the v15 crossing transport pattern."
 ]
}
Path('/mnt/data/congruence_census_double_rank_v29_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False))
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
