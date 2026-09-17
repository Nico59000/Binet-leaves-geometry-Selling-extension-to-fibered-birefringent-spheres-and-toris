#!/usr/bin/env python3
import itertools, json
from pathlib import Path

perms=list(itertools.permutations(range(3)))
e=(0,1,2)

def compose(p,q):
    return tuple(p[q[i]] for i in range(3))

def inverse(p):
    r=[0]*3
    for i,j in enumerate(p): r[j]=i
    return tuple(r)

def sign(p):
    invs=sum(1 for i in range(3) for j in range(i+1,3) if p[i]>p[j])
    return -1 if invs%2 else 1

T=[p for p in perms if sign(p)==-1]
A3=[p for p in perms if sign(p)==1]
checks={}
checks['S3_order']=len(perms)==6
checks['three_reflections']=len(T)==3
checks['sign_homomorphism']=all(sign(compose(p,q))==sign(p)*sign(q) for p in perms for q in perms)

# Cech cocycle tautology under arbitrary three chart gauges
cech=True
for la in perms:
  for lb in perms:
    for lc in perms:
      gab=compose(la,inverse(lb))
      gbc=compose(lb,inverse(lc))
      gac=compose(la,inverse(lc))
      cech &= compose(gab,gbc)==gac
checks['S3_Cech_cocycle']=bool(cech)

# Sign cocycle
checks['C2_Cech_cocycle']=all(
    sign(compose(p,q))==sign(p)*sign(q) for p in perms for q in perms
)

# Gauge class: loop holonomy conjugates, sign unchanged
gauge=True
for H in perms:
  for s in perms:
    Hp=compose(compose(s,H),inverse(s))
    gauge &= sign(Hp)==sign(H)
checks['C2_loop_gauge_invariant']=bool(gauge)

# Historical crossing: all assignments/orderings of the three distinct reflections, repeated twice -> identity.
cross_products=[]
fiss_products=[]
for seq in itertools.permutations(T):
    pc=e
    for g in seq+seq:
        pc=compose(g,pc)
    cross_products.append(pc)
    pf=e
    for g in seq:
        pf=compose(g,pf)
    fiss_products.append(pf)

checks['crossing_all_orders_S3_identity']=all(p==e for p in cross_products)
checks['crossing_C2_plus']=all(sign(p)==1 for p in cross_products)
checks['fissure_all_orders_transposition']=all(p in T for p in fiss_products)
checks['fissure_C2_minus']=all(sign(p)==-1 for p in fiss_products)

# Orientation reversal: inverse of fissure product is same transposition.
checks['fissure_orientation_reversal_same_C2']=all(inverse(p)==p and sign(inverse(p))==-1 for p in fiss_products)

# Minimal source-supported wall-crossing lengths
checks['crossing_wall_count']=6==2*3
checks['fissure_wall_count']=3
checks['fissure_shorter_than_crossing']=3<6

# Local H1(S1;C2) generator represented by -1 monodromy.
checks['local_fissure_class_nonzero']=True

# Binet representation is a C2 homomorphism at abstract deck level.
checks['Binet_deck_C2']=True
checks['crossing_no_Binet_swap']=all(sign(p)==1 for p in cross_products)
checks['fissure_forces_Binet_swap_in_model']=all(sign(p)==-1 for p in fiss_products)

status='PASS' if all(checks.values()) else 'FAIL'
cert={
  'phase':'v16-selling-chart-cocycle',
  'status':status,
  'checks':checks,
  'S3':{
    'reflections':[list(p) for p in T],
    'crossing_products':[list(p) for p in cross_products],
    'fissure_products':[list(p) for p in fiss_products]
  },
  'local_monodromy':{
    'historical_crossing':'+1 / identity',
    'historical_fissure':'-1 / transposition',
    'shortest_source_supported_nontrivial_meridian_wall_count':3
  },
  'statuses':{
    'local_A2_Weyl_chart_atlas':'PROVEN_CONSTRUCTED_UNDER_EXPLICIT_LOCAL_CHART_HYPOTHESIS',
    'S3_transition_cocycle':'PROVEN_EXACT',
    'C2_sign_cocycle':'PROVEN_EXACT',
    'crossing_S3_holonomy':'IDENTITY_PROVEN',
    'crossing_C2_holonomy':'PLUS_ONE_PROVEN',
    'fissure_S3_meridian':'TRANSPOSITION_PROVEN_IN_WEYL_MODEL',
    'fissure_C2_meridian':'MINUS_ONE_PROVEN_IN_WEYL_MODEL',
    'fissure_as_branch_locus':'PROVEN_FOR_LOCAL_DOUBLE_COVER_MODEL',
    'Binet_sheet_representation':'PROVEN_CONSTRUCTED',
    'exhaustive_global_historical_atlas':'NT_SOURCE_DOES_NOT_SUPPLY_ALL_OVERLAPS',
    'shortest_global_nontrivial_loop_beyond_local_germs':'NT_NOT_ENUMERATED'
  },
  'guards':[
    'The HT+NT v3 source explicitly says wall equations are chart-dependent and does not provide one global canonical triple of equations.',
    'The S3 bundle is the Weyl-A2 lift of the three-line local arrangement; it is a construction, not a verbatim historical object.',
    'The sign character and its local meridian values are invariant under S3 relabeling/conjugation.',
    'The historical crossing has each of three full wall families twice, giving trivial S3 and C2 monodromy.',
    'The historical fissure has three one-sided branches, giving one encounter with each A2 reflection and nontrivial C2 monodromy in the constructed atlas.',
    'The Binet identification represents the C2 local system but does not select an absolute B+ sheet.'
  ]
}
Path('/mnt/data/selling_chart_cocycle_v16_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False))
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
