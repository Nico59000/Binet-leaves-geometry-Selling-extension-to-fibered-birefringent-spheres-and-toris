#!/usr/bin/env python3
import json,subprocess,sys,hashlib
from pathlib import Path
B=Path('/mnt/data/v50_work')
scripts=['build_field_orbigroupoid_v50.py','extract_orbifold_links_v50.py','build_orbifold_cochains_fast_v50.py','test_F29_field_functor_gate_v50.py','certify_F29_V4_stack_repair_v50.py']
files=['selling_field_orbigroupoid_v50.json','selling_orbifold_links_v50.json','selling_field_orbifold_cochains_v50.json','selling_F29_field_functor_gate_v50.json','selling_F29_V4_stack_repair_v50.json']
pre={f:hashlib.sha256((B/f).read_bytes()).hexdigest() for f in files}
run=[]
for s in scripts:
 p=subprocess.run([sys.executable,str(B/s)],cwd=B,capture_output=True,text=True)
 run.append({'script':s,'returncode':p.returncode,'stdout_tail':p.stdout[-2000:],'stderr_tail':p.stderr[-1000:]})
 if p.returncode!=0:
  print(json.dumps(run,indent=2));sys.exit(1)
post={f:hashlib.sha256((B/f).read_bytes()).hexdigest() for f in files}
checks={}
def ck(k,v):checks[k]=bool(v)
reg=json.load(open(B/'selling_regular_TU_quotient_v50.json'))
g=json.load(open(B/'selling_field_orbigroupoid_v50.json'))
l=json.load(open(B/'selling_orbifold_links_v50.json'))
c=json.load(open(B/'selling_field_orbifold_cochains_v50.json'))
f=json.load(open(B/'selling_F29_field_functor_gate_v50.json'))
r=json.load(open(B/'selling_F29_V4_stack_repair_v50.json'))
ck('primal_regular_62_123_62',reg['C0']==62 and reg['C1']==123 and reg['C2']==62)
ck('primal_boundary_squared_zero',reg['partial1_partial2_zero'] is True)
ck('dual_62_objects_123_arrows',g['objects']==62 and len(g['arrows'])==123)
ck('arrow_types_102_16_5',g['arrow_counts']=={'ORDINARY-CROSSING':102,'MULTIWALL-ISOTROPY':16,'DECK-MIRROR-ISOTROPY':5})
ck('links_62',l['vertex_count']==62 and len(l['links'])==62)
ck('link_cycle_counts_34_28',l['cycle_count_histogram']=={'2':34,'1':28} or l['cycle_count_histogram']=={2:34,1:28})
ck('cochain_status',c['status']=='PASS-EXACT-ORBIFOLD-COCHAIN-COMPLEX')
ck('binet_dims_62_107_37',c['dimensions']['Binet']=={'C0':62,'C1':107,'C2':37})
ck('ad_dims_186_348_174',c['dimensions']['Ad']=={'C0':186,'C1':348,'C2':174})
ck('binet_d1d0',c['identities']['Binet_d1d0_zero'] is True)
ck('ad_d1d0',c['identities']['Ad_d1d0_zero'] is True)
ck('cohomology_locked',c['cohomology_gate'].startswith('LOCKED'))
S=f['strict_field_functor_to_F29_can'];Q=f['canonical_repair']
ck('strict_F29_refuted_typed',S['status']=='REFUTED-TYPED')
ck('strict_arrow_count_111_123',S['arrow_endpoint_matches_strict']==111 and S['arrow_total']==123)
ck('strict_link_count_59_62',S['link_loops_strictly_closed']==59 and S['link_total']==62)
ck('witness_118_nontrivial_U',S['no_go_witness_edge_118']['edge']==118 and S['no_go_witness_edge_118']['V4_endpoint_witness']=='U' and not S['no_go_witness_edge_118']['strict_F29_endpoint'])
ck('right_action_free',S['free_right_action_exhaustive_check_168x3_nonidentity_V4'] is True)
ck('repair_42_objects',Q['object_orbits']==42 and Q['V4_order']==4)
ck('repair_all_arrows_close',Q['all_123_arrow_endpoints_close_in_quotient'] is True)
ck('repair_all_links_close',Q['all_62_orbifold_links_close_in_quotient'] is True)
ck('stack_repair_status',r['status']=='PASS-EXACT-V4-STACK-REPAIR-2CELL-COMPATIBILITY')
ck('stack_link_hist_59_1_1_1',r['V4_closure_histogram']=={'I':59,'TU':1,'T':1,'U':1})
ck('stack_all_62_fill',r['checks']['all_62_closed_lifts_fill_by_841_fundamental_cells'] is True)
ck('stack_841_basis',r['checks']['F29_fundamental_cell_count']==841)
ck('deterministic_replay_all_5',pre==post)
status='PASS' if all(checks.values()) else 'FAIL'
out={'version':'v50','status':status,'publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE','check_count':len(checks),'pass_count':sum(checks.values()),'checks':checks,'script_runs':run,'sha256':post,'deterministic_replay':pre==post,
'cohomology_gate':'LOCKED because strict F_Sell^field -> F29_can is REFUTED-TYPED; repaired [F29_can/V4] target is a distinct typed object.'}
(B/'v50_orbigroupoid_gate_verification.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({'status':status,'pass':sum(checks.values()),'total':len(checks),'deterministic_replay':pre==post},indent=2))
if status!='PASS':sys.exit(2)
