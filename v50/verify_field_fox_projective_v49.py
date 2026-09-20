#!/usr/bin/env python3
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path

B=Path(__file__).resolve().parent
checks=[]
def ck(name, cond):
    checks.append((name,bool(cond)))

# Replay sealed predecessor verifier.
p=subprocess.run([sys.executable,str(B/'verify_field_fox_projective_v48.py')],cwd=B,text=True,capture_output=True)
try: pred=json.loads(p.stdout)
except Exception: pred={}
ck('predecessor_v48_replay_pass',p.returncode==0 and pred.get('status')=='PASS' and pred.get('combined_checks')==8179)
# This replay gate is infrastructural and is not counted among the 67 new mathematical/schema checks.

# A. Exact constructor self-test: 29 checks.
e=subprocess.run([sys.executable,str(B/'selling_reconstruction_engine_v49.py'),'--selftest'],cwd=B,text=True,capture_output=True)
eng=json.loads(e.stdout)
engine_checks=eng['checks']
assert len(engine_checks)==29
for k,v in engine_checks.items(): ck('engine_'+k,v)

# B. Radius-four exact gate: 12 checks.
r4=json.loads((B/'gamma3_radius4_support_gate_v49.json').read_text())
ck('r4_status_pass',r4['status']=='PASS')
ck('r4_ball_size_10640',r4['exact_group_ball_radius4']['size']==10640)
ck('r4_depth_profile_exact',r4['exact_group_ball_radius4']['depth_profile']=={'0':1,'1':15,'2':147,'3':1213,'4':9264})
ck('r4_depth_profile_sums',sum(r4['exact_group_ball_radius4']['depth_profile'].values())==10640)
ck('r4_mod4_image_256',r4['radius4_image_in_mod4']['distinct_group_elements']==256)
ck('r4_mod4_min_profile_exact',r4['radius4_image_in_mod4']['minimal_depth_profile']=={'0':1,'1':9,'2':36,'3':84,'4':126})
ck('r4_mod4_profile_sums',sum(r4['radius4_image_in_mod4']['minimal_depth_profile'].values())==256)
ck('r4_old_dual_refuted_typed',r4['v48_dual_extension'].startswith('REFUTED-TYPED'))
cc=r4['first_exact_countercolumn']
ck('r4_countercolumn_relator',cc['relator_name']=='R1_F3^2' and cc['relator_index']==2)
ck('r4_countercolumn_support',cc['support_depth']==4 and cc['support_word']==['E13','E23','F1','F2'])
ck('r4_countercolumn_dual_target',cc['old_dual_value']==-2 and cc['target_value']==0)
ck('r4_lift_remains_nt',r4['radius4_lift_existence']=='NT_NOT_DECIDED' and 'not a positive lift' in r4['guard'])

# C. Thirteen radius-two cycles: 10 checks.
ra=json.loads((B/'gamma3_radius2_relation_audit_v49.json').read_text())
ck('r2_lift_count_13',ra['lift_count']==13)
ck('r2_term_count_466',ra['term_count']==466)
ck('r2_relator_coordinates_36',ra['relator_coordinates_used']==36)
ck('r2_critical_rank_13',ra['critical_character']['rank_of_13_lifts']==13)
ck('r2_critical_status',ra['critical_character']['status']=='PROVEN_INDEPENDENT_IN_CRITICAL_CHARACTER_SHADOW')
hist={str(k):v for k,v in ra['all_512_character_rank_histogram'].items()}
ck('r2_histogram_exact',hist=={'9':2,'10':29,'11':113,'12':209,'13':159})
ck('r2_histogram_total_512',sum(hist.values())==512)
ck('r2_rank_range_9_13',ra['minimum_character_rank']==9 and ra['maximum_character_rank']==13)
ck('r2_overlap_edges_71',ra['pairwise_overlap_edges']==71 and len(ra['pairwise_overlap_graph'])==71)
ck('r2_no_private_unit_certificate',ra['private_unit_coordinate_certificates']==[])

# D. H1/H2 compilation contract: 5 checks.
h12=json.loads((B/'selling_reconstruction_to_H12_contract_v49.json').read_text())
ck('h12_contract_version',h12['version']=='v49')
ck('h12_contract_C0_C1_C2',set(h12['required_total_outputs'])=={'C0','C1','C2'})
ck('h12_contract_d1d0',h12['cochains']['mandatory_identity']=='d1*d0 = 0 exactly before any cohomology is published')
ck('h12_contract_F29_three_gates',len(h12['F29_gate'])==3)
ck('h12_contract_current_state_nt',h12['current_v49_state']=='COMPILER-SCHEMA-PROVEN; INPUT_COMPLEX_NOT_YET_TOTAL; HISTORICAL_H1_H2_NT')

# E. Toolchain: 5 checks.
tc=json.loads((B/'selling_figures_1_6_reconstruction_toolchain_v49.json').read_text())
ck('toolchain_version',tc['version']=='v49')
ck('toolchain_seven_stages',len(tc['toolchain'])==7)
ck('toolchain_constructor_bound',tc['toolchain'][1]['artifact']=='selling_reconstruction_engine_v49.py')
ck('toolchain_feasibility_next',tc['toolchain'][2]['status']=='NEXT-IMPLEMENTATION')
ck('toolchain_raster_secondary',tc['toolchain'][5]['status']=='SECONDARY-ONLY')

# F. Figure recipes: 6 checks.
rec=json.loads((B/'selling_figure_recipe_registry_v49.json').read_text())
figs=rec['figures']; f1=figs[0]; f6=figs[5]
ck('recipes_version',rec['version']=='v49')
ck('recipes_six_figures',[x['figure'] for x in figs]==[1,2,3,4,5,6])
ck('recipes_fig1_invariant60',f1['known_exact']['invariant']==60)
ck('recipes_fig1_q0',f1['known_exact']['Q0']==[[12,0,0],[0,-1,0],[0,0,-5]])
ck('recipes_fig1_six_axes_centres',f1['known_exact']['external_symmetry_axes_count']==6 and f1['known_exact']['fourfold_centres_count']==6)
ck('recipes_fig6_counts',f6['known_exact']['external_symmetry_axes_count']==6 and f6['known_exact']['other_external_boundaries_count']==5)

# Count only the v49-specific block, excluding predecessor replay gate.
new_checks=checks[1:]
assert len(new_checks)==67, len(new_checks)
failed=[n for n,v in new_checks if not v]
report={
 'version':'v49',
 'predecessor_checks':8179,
 'new_checks':67,
 'combined_checks':8246,
 'new_passed':67-len(failed),
 'status':'PASS' if not failed and checks[0][1] else 'FAIL',
 'failed':failed,
 'predecessor_replay':pred,
 'new_check_names':[n for n,_ in new_checks]
}
print(json.dumps(report,indent=2,ensure_ascii=False))
raise SystemExit(0 if report['status']=='PASS' else 1)
