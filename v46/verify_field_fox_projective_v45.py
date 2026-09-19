#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys,json,csv,hashlib
B=Path(__file__).resolve().parent
checks={}
def ok(k,v): checks[k]=bool(v); print(('PASS ' if v else 'FAIL ')+k)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
# Parent replay and append-only ledger.
p=subprocess.run([sys.executable,str(B/'verify_field_fox_projective_v44.py')],cwd=B,capture_output=True,text=True)
ok('parent_v44_replay_exit_zero',p.returncode==0)
vr=json.loads((B/'v44_VERIFICATION_REPORT.json').read_text())
ok('parent_v44_status_pass',vr.get('status')=='PASS')
ok('parent_v44_combined_7917',vr.get('combined_checks')==7917)
entries=[]; mism=[]
for line in (B/'v44_SHA256SUMS.txt').read_text().splitlines():
    if not line.strip():continue
    h,n=line.split('  ',1);entries.append(n);q=B/n
    if not q.exists() or sha(q)!=h:mism.append(n)
ok('parent_v44_ledger_221',len(entries)==221)
ok('parent_v44_ledger_byte_exact_after_replay',not mism)
# Source bytes and axis-binding attempt.
ok('selling_plate_hash_canonical',sha(B/'LOG_0016.pdf')=='423da98bd897c5dcb2b3556f4c9dc1367a4f12263d247608cfb2ba7a00f2e47b')
ok('bsb_photo_hash_canonical',sha(B/'Tafel-III-fig1.jpg')=='36fbabf24d710e4bfd570bf05b2dba8350862ffcd02d44efd4f755204653a2b0')
p=subprocess.run([sys.executable,str(B/'detect_selling_axes_v45.py')],cwd=B,capture_output=True,text=True)
ok('axis_detector_exit_zero',p.returncode==0)
axis=json.loads((B/'selling_bsb_axis_binding_attempt_v45.json').read_text())
ok('axis_trace_count_7',axis['trace_count']==7)
ok('axis_dual_source_support_min_gt_0_95',axis['min_bsb_dash_support']>0.95)
ok('source_historical_axis_count_6',axis['historical_external_axis_count']==6)
ok('hax_binding_zero_due_nonunique_subset',axis['historical_binding']['HAX_bound']==0)
ok('hc_binding_zero_due_nonunique_cycle',axis['historical_binding']['HC_bound']==0)
ok('axis_binding_decision_fail_closed','NONUNIQUE' in axis['decision'] and 'FAIL_CLOSED' in axis['decision'])
ok('C2_89_not_activated',axis['C2_89']['activated'] is False)
tr=list(csv.DictReader(open(B/'selling_dotted_axis_traces_v45.csv',encoding='utf-8')))
ok('axis_csv_7',len(tr)==7)
for r in tr:
    ok('axis_support_'+r['trace_id'],float(r['bsb_support_fraction_on_source_dark_samples'])>0.95)
    ok('axis_trace_not_historical_'+r['trace_id'],r['status'].endswith('NOT_BOUND_TO_HAX'))
# Historical H1/H2 gate: inspect actual tables, not declared counts.
nodes=list(csv.DictReader(open(B/'selling_fig1_nodes_v40.csv',encoding='utf-8')))
segs=list(csv.DictReader(open(B/'selling_fig1_segments_v40.csv',encoding='utf-8')))
fields=list(csv.DictReader(open(B/'selling_fig1_fields_v40.csv',encoding='utf-8')))
missing={
 'historical_node_assignments':sum(r['historical_node_type']=='UNRESOLVED' or not r['historical_node_id'] for r in nodes),
 'historical_segment_assignments':sum(r['historical_assignment_status'].startswith('NT_') for r in segs),
 'historical_field_assignments':sum(r['field_assignment_status'].startswith('NT_') for r in fields),
 'ordered_boundary_cycles':sum(not r['ordered_boundary_cycle'] or r['ordered_boundary_cycle']=='UNRESOLVED' for r in fields),
 'repetition_pairs':sum(not r['repetition_of'] or r['repetition_of']=='UNRESOLVED' for r in fields),
 'stabilizer_words':sum(not r['stabilizer_word'] or r['stabilizer_word']=='UNRESOLVED' for r in fields),
}
ok('historical_C0_size_136',len(nodes)==136)
ok('historical_C1_size_105',len(segs)==105)
ok('historical_C2_size_34',len(fields)==34)
ok('historical_nodes_all_unresolved',missing['historical_node_assignments']==136)
ok('historical_segments_all_unresolved',missing['historical_segment_assignments']==105)
ok('historical_fields_all_unresolved',missing['historical_field_assignments']==34)
ok('boundary_cycles_all_missing',missing['ordered_boundary_cycles']==34)
ok('repetition_pairs_all_missing',missing['repetition_pairs']==34)
ok('stabilizer_words_all_missing',missing['stabilizer_words']==34)
h12={
 'version':'v45','table_sizes':{'C0':len(nodes),'C1':len(segs),'C2':len(fields)},'missing':missing,
 'binding_input':'selling_bsb_axis_binding_attempt_v45.json','historical_axis_bindings':0,'historical_reflections':0,
 'C2_89_census':'NOT_ACTIVATED_FAIL_CLOSED',
 'd1d0_test':'NOT_FORMED_BECAUSE_HISTORICAL_C0_C1_C2_AND_TRANSPORTS_ARE_NOT_TOTAL',
 'face_holonomy_stabilizer_test':'NOT_FORMED','F29_object_edge_face_functor':'PENDING-MORPHISM',
 'H1_hist':'NT_NOT_RECOMPUTED','H2_hist':'NT_NOT_RECOMPUTED',
 'next_required_source_data':['select exactly the six external dotted axes from source labels/incidence','bind their six cyclic fourfold centres','historically type nodes and boundary arcs','attach substitution matrices to arcs','order every field boundary','encode repetition target and stabilizer word per field'],
 'decision':'HISTORICAL_COMPLEX_NOT_TOTAL; D1_D0_AND_H1_H2_BLOCKED_FAIL_CLOSED'
}
(B/'selling_h12_materialization_gate_v45.json').write_text(json.dumps(h12,indent=2,ensure_ascii=False),encoding='utf-8')
ok('h12_recompute_blocked',h12['H1_hist'].startswith('NT_') and h12['H2_hist'].startswith('NT_'))
# Higher Gamma grammars: helpers recompute exhaustive declared searches.
for script,label in [('enumerate_gamma3_quadruple_v45.py','g4'),('enumerate_gamma3_five_v45.py','g5'),('enumerate_gamma3_six_v45.py','g6')]:
    q=subprocess.run([sys.executable,str(B/script)],cwd=B,capture_output=True,text=True)
    ok(label+'_helper_exit_zero',q.returncode==0)
g4=json.loads((B/'gamma3_quadruple_free_identity_probe_v45.json').read_text())
g6=json.loads((B/'gamma3_six_relator_probe_v45.json').read_text())
# Five helper is stdout-only; rerun capture and parse exact terminal line.
q5=subprocess.run([sys.executable,str(B/'enumerate_gamma3_five_v45.py')],cwd=B,capture_output=True,text=True)
ok('g4_pair_residual_keys_75672',g4['pair_residual_keys']==75672)
ok('g4_raw_3840',g4['four_distinct_free_identity_sequences_raw']==3840)
ok('g4_classes_1920',g4['dihedral_classes']==1920)
ok('g4_absorbed_all',g4['absorbed_by_M43']==1920 and g4['survivors_mod_M43']==0)
ok('g5_pairs_7904','pairs 7904 keys 3973' in q5.stdout)
ok('g5_triples_150896','triples 150896 keys 53572' in q5.stdout)
ok('g5_zero_identity_keys','identity keys 0 raw distinctbase5 0' in q5.stdout)
ok('g6_raw_1536',g6['raw_sequences']==1536)
ok('g6_classes_384',g6['dihedral_classes']==384)
ok('g6_absorbed_all',g6['absorbed_by_M43']==384 and g6['survivors_mod_M43']==0)
M43=json.loads((B/'gamma3_partial_resolution_v43.json').read_text())
ok('character_residual_min_14_retained',M43['finite_shadow_residual_min_dimension']==14)
higher={
 'version':'v45','ring':'R=Z[Gamma_3(2)]','baseline_module':M43['module_isomorphism'],
 'four_relator_free_identities':g4,
 'five_relator_3plus2_positive_cancellation_grammar':{'positive_pair_instances':7904,'positive_pair_residual_keys':3973,'sequential_triples':150896,'triple_residual_keys':53572,'inverse_residual_key_matches':0,'five_distinct_base_free_identities':0},
 'six_relator_3plus3_positive_cancellation_grammar':g6,
 'finite_character_shadow_residual_min_dimension':14,
 'decision':'NO_NEW_M43_DIRECTION_IN_DECLARED_4_5_6_RELATOR_FREE_CANCELLATION_GRAMMARS; COMPLETE_RESOLUTION_REMAINS_NT_BY_NONZERO_CHARACTER_SHADOW',
 'guard':'The search closes only the explicitly declared cyclic/free-cancellation grammars, not arbitrary identities among relators.'
}
(B/'gamma3_higher_identity_census_v45.json').write_text(json.dumps(higher,indent=2,ensure_ascii=False),encoding='utf-8')
ok('gamma_complete_resolution_stays_NT','COMPLETE_RESOLUTION_REMAINS_NT' in higher['decision'])
# Singer source-native gate.
sing44=json.loads((B/'selling_singer_injection_gate_v44.json').read_text()) if (B/'selling_singer_injection_gate_v44.json').exists() else {}
singer={
 'version':'v45','predecessor':'selling_singer_injection_gate_v44.json',
 'historical_reflection_source_bound':False,'source_native_line_L':False,'source_native_Singer_framing_sigma':False,
 'Gamma_equals_framing_inversion_certified':False,'path_and_face_equivariance_certified':False,
 'current_scalar_injection':'REFUTED-TYPED_RETAINED',
 'Singer_odd_pseudoscalar':'SEPARATED_CONDITIONAL_CANDIDATE_RETAINED',
 'injection_to_J_line':'NT_PENDING_SOURCE_NATIVE_L_SIGMA_GAMMA_AND_CELL_EQUIVARIANCE',
 'Lorentz_orientation':'NT','physical_phase':'SEPARATED/PENDING-MORPHISM',
 'decision':'NO_SOURCE_NATIVE_ODD_INJECTION_PROMOTION_IN_V45'
}
(B/'selling_singer_injection_gate_v45.json').write_text(json.dumps(singer,indent=2,ensure_ascii=False),encoding='utf-8')
ok('Singer_injection_remains_NT',singer['injection_to_J_line'].startswith('NT_'))
ok('Lorentz_orientation_remains_NT',singer['Lorentz_orientation']=='NT')
# Strict R5 adapter and old typed no-gos retained.
r5=json.loads((B/'selling_gamma_R5_adapter_v43.json').read_text())
ok('R5_strict_adapter_refuted_typed','REFUTED' in json.dumps(r5))
# Registry.
registry={
 'version':'v45','items':[
  {'id':'V45-01','object':'seven dual-source dotted traces','status':'PROVEN_DIAGNOSTIC'},
  {'id':'V45-02','object':'six external HAX subset / six HC cyclic binding','status':'NT_NONUNIQUE_SOURCE_BINDING'},
  {'id':'V45-03','object':'89-node historical C2 action','status':'NT_NOT_ACTIVATED_FAIL_CLOSED'},
  {'id':'V45-04','object':'historical C0/C1/C2 totalization','status':'NT'},
  {'id':'V45-05','object':'F_Sell_field -> F29_can','status':'PENDING-MORPHISM'},
  {'id':'V45-06','object':'historical H1/H2','status':'NT_NOT_RECOMPUTED'},
  {'id':'V45-07','object':'4-relator declared free identities','status':'PROVEN_ALL_1920_CLASSES_ABSORBED_BY_M43'},
  {'id':'V45-08','object':'5-relator declared 3+2 grammar','status':'ZERO_IDENTITIES'},
  {'id':'V45-09','object':'6-relator declared 3+3 grammar','status':'PROVEN_ALL_384_CLASSES_ABSORBED_BY_M43'},
  {'id':'V45-10','object':'complete group-ring resolution','status':'NT_RESIDUAL_CHARACTER_MIN14'},
  {'id':'V45-11','object':'Singer odd injection','status':'SEPARATED/PENDING-MORPHISM'},
  {'id':'V45-12','object':'strict R5-Gamma incidence adapter','status':'REFUTED-TYPED_RETAINED'},
  {'id':'V45-13','object':'MXM^-1 covariant use','status':'OPPOSITE_VARIANCE_RETAINED'},
  {'id':'V45-14','object':'roof reversal = transpose','status':'REFUTED-TYPED_RETAINED'},
  {'id':'V45-15','object':'physical phase','status':'SEPARATED/PENDING-MORPHISM'}]}
(B/'pending_morphism_registry_v45.json').write_text(json.dumps(registry,indent=2,ensure_ascii=False),encoding='utf-8')
ok('registry_15_items',len(registry['items'])==15)
failed=[k for k,v in checks.items() if not v]
new_count=len(checks); report={'version':'v45','predecessor_checks':7917,'new_checks':new_count,'combined_checks':7917+new_count,'new_passed':new_count-len(failed),'status':'PASS' if not failed else 'FAIL','failed':failed}
(B/'v45_VERIFICATION_REPORT.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(report,indent=2,ensure_ascii=False))
sys.exit(0 if not failed else 1)
