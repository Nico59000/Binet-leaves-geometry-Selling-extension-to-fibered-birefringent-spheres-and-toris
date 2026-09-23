#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import json, hashlib
B=Path('/mnt/data/v51_work')
Z=json.loads((B/'selling_historical_cohomology_Z_2primary_audit_v50.json').read_text())
I=json.loads((B/'selling_integral_orbifold_cochains_v50.json').read_text())
L=json.loads((B/'selling_orbifold_links_v50.json').read_text())
G=json.loads((B/'selling_regular_TU_quotient_v50.json').read_text())
links={int(x['vertex']):x for x in L['links']}

def c1_map(sys):
    out=[]
    for rec in I['systems'][sys]['C1_saturated_bases']:
        e=int(rec['edge']); rank=int(rec['rank']); basis=rec['basis']
        for j in range(rank):
            vec=[row[j] for row in basis] if basis else []
            out.append({'edge':e,'local_index':j,'basis_vector':vec,'edge_kind':rec['kind']})
    assert len(out)==I['systems'][sys]['dimensions']['C1']
    return out

def sparse_vec(items,n):
    v=[0]*n
    for z in items:v[int(z['i'])]=int(z['value'])
    return v

def mat_from_sparse(entries,nrow,ncol):
    rows=[{} for _ in range(nrow)]
    for z in entries: rows[int(z['row'])][int(z['col'])]=int(z['value'])
    return rows

def dot_sparse(row,v):return sum(c*v[j] for j,c in row.items())

out={'version':'v51','status':'PASS-EXACT-TORSION-SUPPORT/BOCKSTEIN/LINK-PAIRING','source_version':'v50 CLOSED/APPEND-ONLY','systems':{},'cross_system':{},'binet_degree2_transfer_defect':{},'guard':[]}
for sys in ['Binet','Ad']:
    n0=I['systems'][sys]['dimensions']['C0']; n1=I['systems'][sys]['dimensions']['C1']; n2=I['systems'][sys]['dimensions']['C2']
    w=Z['systems'][sys]['explicit_H1_order2_witness']
    z=sparse_vec(w['z_in_C0_mod2_lift_support'],n0); tau=sparse_vec(w['tau_equals_d0_z_over_2_support'],n1)
    D0=mat_from_sparse(I['systems'][sys]['d0_sparse'],n1,n0)
    D1=mat_from_sparse(I['systems'][sys]['d1_sparse'],n2,n1)
    d0z=[dot_sparse(r,z) for r in D0]
    assert all(x%2==0 for x in d0z)
    assert [x//2 for x in d0z]==tau
    assert all(dot_sparse(r,tau)==0 for r in D1)
    # mod2 z is H0 cocycle and beta0 is tau
    cmap=c1_map(sys)
    edge_support=defaultdict(list)
    for i,x in enumerate(tau):
        if x:
            rec=dict(cmap[i]);rec.update({'coordinate':i,'coefficient':x});edge_support[rec['edge']].append(rec)
    # link V48 restriction
    lk=links[48]['chosen']['word']
    link_edges=[int(t['edge']) for t in lk]
    link_support=sorted(set(link_edges)&set(edge_support))
    # rows of C2 for vertex 48: derive offsets by ranks
    c2off=0; v48_rows=[]
    for rec in I['systems'][sys]['C2_saturated_left_fixed_bases']:
        rk=int(rec['rank']); v=int(rec['vertex'])
        if v==48:v48_rows=list(range(c2off,c2off+rk))
        c2off+=rk
    vals=[dot_sparse(D1[r],tau) for r in v48_rows]
    assert all(v==0 for v in vals)
    out['systems'][sys]={
      'H1_order2':'Z/2','H0_mod2_source_support':w['z_in_C0_mod2_lift_support'],
      'beta0_statement':'beta_0([z mod 2])=[d0 z / 2]=[tau] is the unique H1 order-2 class',
      'tau_C1_coordinate_support':w['tau_equals_d0_z_over_2_support'],
      'tau_edge_orbit_support_count':len(edge_support),
      'tau_edge_orbit_support':[{ 'edge':e,'coordinates':edge_support[e]} for e in sorted(edge_support)],
      'V48_link_edges':link_edges,'V48_link_supported_edges':link_support,
      'V48_link_support_count':len(link_support),'V48_d1_pairing_values':vals,
      'V48_boundary_pairing':'ZERO_BY_COCYCLE_CONDITION',
      'bockstein_checks':{'d0z_even':True,'tau_equals_half_d0z':True,'d1tau_zero':True}
    }
# support-only cross comparison, explicitly not a coefficient pairing
EB=set(x['edge'] for x in out['systems']['Binet']['tau_edge_orbit_support'])
EA=set(x['edge'] for x in out['systems']['Ad']['tau_edge_orbit_support'])
out['cross_system']={
 'Binet_edge_support_count':len(EB),'Ad_edge_support_count':len(EA),'common_edge_support_count':len(EB&EA),
 'common_edge_support':sorted(EB&EA),'union_edge_support_count':len(EB|EA),
 'status':'SUPPORT-INTERSECTION-ONLY/NOT-A-COHOMOLOGICAL-PAIRING',
 'reason':'No natural coefficient morphism or bilinear pairing Binet x Ad -> target has been constructed.'
}
# defect
lk=links[48]['chosen']
out['binet_degree2_transfer_defect']={
 'carrier_vertex':48,'monodromy':'TU','transfer_index':2,'fixed_lattice_rank':1,
 'downstairs_H2_torsion':False,'upstairs_invariant_H2_extra_torsion':'Z/2',
 'link_edge_word':lk['word'],'link_edge_ids':[int(t['edge']) for t in lk['word']],
 'pairing_with_Binet_H1_tau':{'typed_as':'d1 row evaluation on tau','value':out['systems']['Binet']['V48_d1_pairing_values'],'result':'ZERO'},
 'pairing_with_Ad_H1_tau':{'typed_as':'support incidence only; coefficient systems differ','common_link_edges':out['systems']['Ad']['V48_link_supported_edges'],'status':'NOT-A-COHOMOLOGICAL-PAIRING'},
 'separation':'DO-NOT-IDENTIFY-H1-TORSION-CARRIERS-WITH-H2-TRANSFER-DEFECT'
}
out['guard']=[
 'The two beta_0 classes live in different coefficient systems and are not identified.',
 'The V48 transfer defect is a degree-2 descent-index carrier, not the same object as either H1 torsion class.',
 'Cross-system edge support intersection is diagnostic only unless a natural coefficient pairing is constructed.'
]
p=B/'selling_torsion_support_bockstein_pairings_v51.json';p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({
 'status':out['status'],
 'Binet_tau_coords':len(out['systems']['Binet']['tau_C1_coordinate_support']),
 'Binet_tau_edges':out['systems']['Binet']['tau_edge_orbit_support_count'],
 'Ad_tau_coords':len(out['systems']['Ad']['tau_C1_coordinate_support']),
 'Ad_tau_edges':out['systems']['Ad']['tau_edge_orbit_support_count'],
 'common_edges':out['cross_system']['common_edge_support_count'],
 'V48_Binet_overlap':out['systems']['Binet']['V48_link_support_count'],
 'V48_Ad_overlap':out['systems']['Ad']['V48_link_support_count'],
 'sha256':hashlib.sha256(p.read_bytes()).hexdigest()
},indent=2))
