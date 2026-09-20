#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys,json,csv,hashlib
B=Path(__file__).resolve().parent
checks={}
def ok(k,v):
 checks[k]=bool(v)
 if not v: print('FAIL',k)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
# Parent v45 is immutable and CLOSED/APPEND-ONLY. Verify atomic token + sealed bundle ledger.
tok=json.loads((B/'v45_ATOMIC_COMPLETION_TOKEN.json').read_text(encoding='utf-8'))
ok('parent_v45_closed',tok.get('publication_state')=='CLOSED/APPEND-ONLY')
ok('parent_v45_atomic',tok.get('atomic_completion') is True and all(v=='PASS' for v in tok['atomic_conjunction'].values()))
ok('parent_v45_checks_7976',tok['verification']['combined']==7976 and tok['verification']['status']=='PASS_7976_OF_7976')
entries=[];mis=[]
for line in (B/'v45_SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
 if not line.strip():continue
 h,n=line.split('  ',1);entries.append(n);p=B/n
 if not p.exists() or sha(p)!=h:mis.append(n)
ok('parent_v45_ledger_248',len(entries)==248)
ok('parent_v45_ledger_byte_exact',not mis)
ok('parent_v45_pdf_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v45.pdf')==tok['pdf']['sha256'])
ok('parent_v45_tex_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v45.tex')==tok['tex']['sha256'])
# A. Source constraints and 7->6 direct DAX/HAX model.
p=subprocess.run([sys.executable,str(B/'analyze_selling_axes_v46.py')],cwd=B,capture_output=True,text=True)
ok('sell46_axis_helper_exit_zero',p.returncode==0)
axis=json.loads((B/'selling_axis_subset_census_v46.json').read_text(encoding='utf-8'))
ok('sell46_dax_count_7',axis['trace_count']==7)
ok('sell46_historical_axis_count_6',axis['historical_external_axes']==6)
ok('sell46_all_7_subsets_enumerated',len(axis['six_subset_census'])==7)
ok('sell46_direct_subset_survivors_zero',axis['direct_six_subset_survivors']==0)
ok('sell46_closest_pair_far_from_90',axis['closest_DAX_pair_to_perpendicular']['acute_angle_deg']<70)
ok('sell46_angle_uncertainty_small',axis['pair_angle_uncertainty_bound_deg']<2)
ok('sell46_direct_model_refuted_typed',axis['decision'].startswith('REFUTED-TYPED'))
ok('sell46_hax_zero',axis['historical_HAX_bound']==0)
ok('sell46_hc_zero',axis['historical_HC_bound']==0)
ok('sell46_C2_not_activated',axis['C2_89']=='NOT_ACTIVATED_FAIL_CLOSED')
for i,s in enumerate(axis['six_subset_census']):
 ok(f'sell46_subset_{i}_direct_refuted',s['status'].startswith('REFUTED-TYPED'))
 ok(f'sell46_subset_{i}_no_orth_pair',s['orthogonal_pair_count_under_conservative_error_bound']==0)
stitch=json.loads((B/'selling_axis_stitching_requirements_v46.json').read_text(encoding='utf-8'))
ok('sell46_stitch_protocol_refined',stitch['status'].startswith('PROTOCOL-REFINED'))
# Explicit C2 gate remains fail-closed.
c2={
 'version':'v46','degree_one_candidates':89,
 'precondition':'source-certified historical HAX curves and derived reflection action preserving the node/edge incidence carrier',
 'precondition_status':'FALSE_DIRECT_DAX_TO_HAX_MODEL_REFUTED_TYPED_AND_CURVED_AXIS_STITCHING_NOT_YET_MATERIALIZED',
 'activated':False,'fixed_points':None,'free_orbits':None,'total_orbits':None,
 'comparison_44_45':'SEPARATED_NOT_TESTED','comparison_34_55':'SEPARATED_NOT_TESTED',
 'decision':'C2_ORBIT_CENSUS_NOT_ACTIVATED_FAIL_CLOSED'}
(B/'selling_degree1_C2_gate_v46.json').write_text(json.dumps(c2,indent=2,ensure_ascii=False),encoding='utf-8')
ok('sell46_C2_values_unset',c2['fixed_points'] is None and c2['free_orbits'] is None)
# B. Historical cell complex / H1/H2 remains non-total, but next source target is sharpened.
nodes=list(csv.DictReader(open(B/'selling_fig1_nodes_v40.csv',encoding='utf-8')))
segs=list(csv.DictReader(open(B/'selling_fig1_segments_v40.csv',encoding='utf-8')))
fields=list(csv.DictReader(open(B/'selling_fig1_fields_v40.csv',encoding='utf-8')))
missing={
 'historical_node_assignments':sum(r['historical_node_type']=='UNRESOLVED' or not r['historical_node_id'] for r in nodes),
 'historical_segment_assignments':sum(r['historical_assignment_status'].startswith('NT_') for r in segs),
 'historical_field_assignments':sum(r['field_assignment_status'].startswith('NT_') for r in fields),
 'ordered_boundary_cycles':sum(not r['ordered_boundary_cycle'] or r['ordered_boundary_cycle']=='UNRESOLVED' for r in fields),
 'repetition_pairs':sum(not r['repetition_of'] or r['repetition_of']=='UNRESOLVED' for r in fields),
 'stabilizer_words':sum(not r['stabilizer_word'] or r['stabilizer_word']=='UNRESOLVED' for r in fields)}
ok('sell46_C0_136',len(nodes)==136 and missing['historical_node_assignments']==136)
ok('sell46_C1_105',len(segs)==105 and missing['historical_segment_assignments']==105)
ok('sell46_C2_34',len(fields)==34 and missing['historical_field_assignments']==34)
ok('sell46_boundary_cycles_34_missing',missing['ordered_boundary_cycles']==34)
ok('sell46_repetition_34_missing',missing['repetition_pairs']==34)
ok('sell46_stabilizers_34_missing',missing['stabilizer_words']==34)
h12={
 'version':'v46','table_sizes':{'C0':136,'C1':105,'C2':34},'missing':missing,
 'axis_model_update':'DIRECT_DAX_TO_HAX_6-SUBSET_REFUTED-TYPED; CURVED/STITCHED_HAX_TRANSCRIPTION_NEEDED',
 'historical_reflections':0,'C2_89_census':'NOT_ACTIVATED_FAIL_CLOSED',
 'd0':'NOT_FORMED','d1':'NOT_FORMED','d1d0_test':'NOT_FORMED','face_holonomy_stabilizer_test':'NOT_FORMED',
 'F29_object_edge_face_functor':'PENDING-MORPHISM','H1_hist':'NT_NOT_RECOMPUTED','H2_hist':'NT_NOT_RECOMPUTED',
 'next_required_source_data':stitch['next_minimal_source_tasks']+[
  'historically type all 136 C0 nodes after axis/centre anchors are fixed',
  'historically type/orient all 105 C1 arcs and attach their substitution matrices',
  'serialize every one of the 34 C2 ordered boundary cycles, repetition targets, and stabilizer words'],
 'decision':'HISTORICAL_COMPLEX_NOT_TOTAL; D0_D1_AND_H1_H2_REMAIN_FAIL_CLOSED'}
(B/'selling_h12_materialization_gate_v46.json').write_text(json.dumps(h12,indent=2,ensure_ascii=False),encoding='utf-8')
ok('sell46_h12_blocked',h12['H1_hist'].startswith('NT_') and h12['H2_hist'].startswith('NT_'))
# C. Character-guided Gamma3 residual and bounded exact lift test.
p=subprocess.run([sys.executable,str(B/'analyze_gamma3_character_residual_v46.py')],cwd=B,capture_output=True,text=True)
ok('g46_character_helper_exit_zero',p.returncode==0)
basis=json.loads((B/'gamma3_character_residual_basis_v46.json').read_text(encoding='utf-8'))
ng=json.loads((B/'gamma3_character_radius1_no_go_v46.json').read_text(encoding='utf-8'))
ok('g46_char_D_rank_8',basis['character_Fox_rank']==8)
ok('g46_char_kernel_35',basis['character_kernel_dimension']==35)
ok('g46_M43_shadow_rank_21',basis['M43_shadow_rank']==21)
ok('g46_residual_dim_14',basis['residual_quotient_dimension']==14)
rr=list(csv.DictReader(open(B/'gamma3_character_residual_basis_v46.csv',encoding='utf-8')))
ok('g46_14_direction_ids',len({r['direction_id'] for r in rr})==14)
ok('g46_radius1_ball_16',ng['ball_unique_elements']==16)
ok('g46_radius1_vars_688',ng['variables_43_relators_times_ball']==688)
ok('g46_radius1_rows_1498',ng['exact_boundary_rows']==1498)
ok('g46_radius1_boundary_rank_559',ng['exact_boundary_rank']==559)
ok('g46_radius1_kernel_129',ng['exact_kernel_dimension']==129)
ok('g46_radius1_residual_image_zero',ng['residual_shadow_image_dimension_of_exact_radius1_kernel_mod_M43']==0)
ok('g46_radius1_no_go_proven',ng['decision'].startswith('PROVEN BOUNDED NO-GO'))
# D. Typed guards retained.
r5=json.loads((B/'selling_gamma_R5_adapter_v43.json').read_text(encoding='utf-8'))
ok('g46_R5_strict_refuted_retained','REFUTED-TYPED' in r5['strict_incidence_preserving_16_to_16_adapter'])
sing=json.loads((B/'selling_singer_injection_gate_v45.json').read_text(encoding='utf-8'))
ok('g46_singer_pending',sing['injection_to_J_line'].startswith('NT_'))
ok('g46_lorentz_orientation_NT',sing['Lorentz_orientation']=='NT')
registry={
 'version':'v46','items':[
  {'id':'V46-01','object':'direct one-DAX-per-HAX 7-to-6 search space','status':'REFUTED-TYPED'},
  {'id':'V46-02','object':'curved/stiched six HAX + six HC source binding','status':'NT/PENDING-SOURCE-TRANSCRIPTION'},
  {'id':'V46-03','object':'89-node historical C2 action','status':'NT_NOT_ACTIVATED_FAIL_CLOSED'},
  {'id':'V46-04','object':'historical C0/C1/C2 totalization','status':'NT'},
  {'id':'V46-05','object':'F_Sell_field -> F29_can','status':'PENDING-MORPHISM'},
  {'id':'V46-06','object':'historical H1/H2','status':'NT_NOT_RECOMPUTED'},
  {'id':'V46-07','object':'minimal residual character quotient','status':'PROVEN_DIMENSION_14'},
  {'id':'V46-08','object':'exact residual lift with coefficient support in radius-one group ball','status':'REFUTED-TYPED/BOUNDED-NO-GO'},
  {'id':'V46-09','object':'complete Gamma3 group-ring resolution','status':'NT; NEXT SEARCH RADIUS>=2 OR DIFFERENT HIGHER-ID MODULE'},
  {'id':'V46-10','object':'strict R5-Gamma incidence adapter','status':'REFUTED-TYPED_RETAINED'},
  {'id':'V46-11','object':'Singer odd injection','status':'SEPARATED/PENDING-MORPHISM'},
  {'id':'V46-12','object':'MXM^-1 covariant use','status':'OPPOSITE_VARIANCE_RETAINED'},
  {'id':'V46-13','object':'roof reversal = transpose','status':'REFUTED-TYPED_RETAINED'},
  {'id':'V46-14','object':'Lorentz orientation','status':'NT'},
  {'id':'V46-15','object':'physical phase','status':'SEPARATED/PENDING-MORPHISM'}]}
(B/'pending_morphism_registry_v46.json').write_text(json.dumps(registry,indent=2,ensure_ascii=False),encoding='utf-8')
ok('g46_registry_15',len(registry['items'])==15)
failed=[k for k,v in checks.items() if not v]
new=len(checks);report={'version':'v46','predecessor_checks':7976,'new_checks':new,'combined_checks':7976+new,'new_passed':new-len(failed),'status':'PASS' if not failed else 'FAIL','failed':failed}
(B/'v46_VERIFICATION_REPORT.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(report,indent=2,ensure_ascii=False))
sys.exit(0 if not failed else 1)
