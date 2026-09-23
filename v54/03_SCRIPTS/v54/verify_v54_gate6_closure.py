#!/usr/bin/env python3
from pathlib import Path
import json, hashlib, subprocess, sys
ROOT=Path(__file__).resolve().parent
INP=ROOT/'v54_gate6_closure_inputs'
CP=json.loads((ROOT/'v54_GATE6_CLOSURE_CHECKPOINT.json').read_text(encoding='utf-8'))
checks=[]
def ok(name,cond):
    checks.append((name,bool(cond)))
    if not cond: raise AssertionError(name)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

# 1. Hash-lock all normalized closure inputs and inherited verifiers.
for name,h in CP['input_sha256'].items():
    p=INP/name
    ok('input_exists_'+name,p.exists())
    ok('input_hash_'+name,sha(p)==h)
for key,name in {'base_gate6':'verify_v54_gate6.py','hypertensor_bridge':'verify_v54_gate6_hypertensor_bridge.py','p46_global':'verify_v54_gate6_p46_global_repair.py'}.items():
    p=ROOT/name
    ok('verifier_exists_'+key,p.exists())
    ok('verifier_hash_'+key,sha(p)==CP['verifier_sha256'][key])

# 2. Destructive-style research replay of all controlling Gate-6 verifier branches.
for key,name in [('base','verify_v54_gate6.py'),('hyper','verify_v54_gate6_hypertensor_bridge.py'),('p46','verify_v54_gate6_p46_global_repair.py')]:
    p=subprocess.run([sys.executable,str(ROOT/name)],cwd=ROOT,capture_output=True,text=True,timeout=240)
    ok('replay_'+key+'_returncode',p.returncode==0)

base=json.loads((INP/'v54_GATE6_VERIFICATION_REPORT.json').read_text())
hyp=json.loads((INP/'v54_GATE6_HYPERTENSOR_STRICT_BRIDGE_REPORT.json').read_text())
p46=json.loads((INP/'v54_GATE6_P46_GLOBAL_REPAIR_REPORT.json').read_text())
p46cp=json.loads((INP/'v54_GATE6_P46_GLOBAL_REPAIR_CHECKPOINT.json').read_text())
state=json.loads((INP/'v54_GATE6_CURRENT_RESEARCH_STATE.json').read_text())
figobs=json.loads((INP/'selling_fig2_fig3_raster_binding_obstruction_v54.json').read_text())
r2inj=json.loads((INP/'r2pure13_to_sellingdeck17_object_injective_obstruction_v54.json').read_text())
r2q=json.loads((INP/'r2pure13_to_sellingdeck17_noninjective_path_quotient_v54.json').read_text())
v53=(INP/'v53_DELIVERY_INDEX.txt').read_text(encoding='utf-8')

# 3. Predecessor and inherited research verifier status.
ok('v53_closed_append_only','STATE: CLOSED / APPEND-ONLY' in v53)
ok('v53_atomic_contract','DOC ∧ MATH ∧ VERIFY ∧ PDF-QA ∧ BUNDLE-TREE ∧ SEAL ∧ DETERMINISTIC-RESEAL ∧ ROUNDTRIP ∧ POSTSEAL = PASS' in v53)
ok('base_31_31',base['status']=='PASS' and base['checks']==31 and base['passed']==31 and not base['failed'])
ok('hyper_all_pass',hyp['status']=='PASS' and hyp['checks']==hyp['passed'] and not hyp['failed'])
ok('p46_27_27',p46['status']=='PASS' and p46['checks']==27 and p46['passed']==27 and not p46['failed'])

# 4. Pascal/Noether and centered amplitude controlling promotions.
P=p46cp['Pascal_p46_global_replay']
ok('p46_full15_global',P['full_15_coefficient_extension']['status']=='PROVEN-GLOBAL-POLYNOMIAL-IDENTITY-IN-15-COEFFICIENTS')
ok('p46_identity',P['full_15_coefficient_extension']['identity']=='rho=(5/2)p_corrected')
ok('p46_gl3',P['GL3_covariance']['status']=='PASS-EXACT')
ok('p46_spinta',P['spinta_normalization'].startswith('PROVEN'))
ok('centered_amplitude',p46cp['centered_amplitude_adapter']['status']=='PROVEN-EXACT-CENTERED-JET-ADAPTER')
ok('no_time_identification',p46cp['centered_amplitude_adapter']['time_identification']=='SEPARATED/NOT-ASSERTED')

# 5. Fig2 exact no-go scope and retained source points.
F2=state['Fig2']
ok('fig2_P0',F2['retained_source_points']['P0']['status']=='PROVEN-SOURCE-RASTER')
ok('fig2_P1',F2['retained_source_points']['P1']['status']=='PROVEN-SOURCE-RASTER')
ok('fig2_P2',F2['retained_source_points']['P2']['status']=='PROVEN-SOURCE-DERIVED-BY-S')
ok('fig2_D5_order10',F2['current_typed_raster_carrier']['label_and_relation_type_preserving_automorphism_group_order']==10)
ok('fig2_current_unique_refuted',F2['point_level_binding_status'].startswith('REFUTED-TYPED'))
ok('fig2_S_NT',F2['S_specific_raster_edge'].startswith('NT'))
ok('fig2_Phi_NT',F2['Phi_raster_species'].startswith('NT'))
ok('fig2_O1_undecided',F2['O1_to_2_or_3']=='UNDECIDED')
FO=figobs['Fig2']
ok('fig2_two_distinct_maps',FO['maps_distinct'] is True)
ok('fig2_mapA_preserved',FO['map_A_audit']['preserved'] is True and FO['map_A_audit']['failures']==[])
ok('fig2_mapB_preserved',FO['map_B_audit']['preserved'] is True and FO['map_B_audit']['failures']==[])
ok('fig2_15_edges_each',FO['map_A_audit']['checked_edges']==15 and FO['map_B_audit']['checked_edges']==15)

# 6. Fig3 exact fail-closed frontier.
F3=state['Fig3']
ok('fig3_root19',F3['root_class']==19)
ok('fig3_exact_2_91',F3['direct_raster_bound']=='2/91')
ok('fig3_exact_89_91',F3['not_propagated']=='89/91')
ok('fig3_exact_two_edges',F3['retained_direct_raster_edges']==['19 -k-> 19','19 -n-> 19'])
ok('fig3_no_symmetry_propagation','NO-SYMMETRY-PROPAGATION' in F3['status'])

# 7. R2PURE13 exact obstruction + guarded noninjective witness.
ok('r2_injective_refuted',r2inj['status'].startswith('REFUTED-TYPED-OBJECT-INJECTIVE'))
ok('r2_source_13_connected',r2inj['source_graph']['vertices']==13 and r2inj['source_graph']['component_sizes']==[13])
ok('r2_source_54_edges',r2inj['source_graph']['edges']==54)
ok('r2_source_107_triangles',r2inj['source_graph']['triangle_count']==107)
ok('r2_target_17',r2inj['target_graph']['vertices']==17)
ok('r2_target_components_458',r2inj['target_graph']['component_sizes']==[4,5,8])
ok('r2_target_triangle_free',r2inj['target_graph']['triangle_count']==0)
ok('r2_strict_edge_refuted',r2inj['strict_edge_incidence_obstruction']['status'].startswith('REFUTED-TYPED'))
ok('r2_noninjective_exists',r2q['decision']['explicit_noninjective_edge_to_path_functor_on_1_skeleton'].startswith('PROVEN-EXISTS'))
ok('r2_54_congruences',r2q['quotient']['all_54_target_congruences']=='PASS' and len(r2q['quotient']['edge_images'])==54)
ok('r2_noninjective_true',r2q['quotient']['noninjective'] is True)
ok('r2_missing_source_2cells',r2q['source']['declared_2cell_presentation']=='NOT-MATERIALIZED-IN-FROZEN-R2PURE-INPUTS')
ok('r2_2cell_NT',r2q['decision']['source_2cell_compatibility'].startswith('NT/'))
ok('r2_full_separated',r2q['decision']['R2PURE13_to_SellingDeck17_adapter']=='SEPARATED/PENDING-MORPHISM')
ok('r2_whiskering_shadow_zero',r2q['whiskering_shadow']['status']=='ZERO-ALL-13-LIFTS' and r2q['whiskering_shadow']['odd_relator_augmentation_residuals']==0)

# 8. Frozen guards and closure semantics.
ok('r4_not_opened',state['retained_closed_guards']['r4']=='NOT-OPENED')
ok('B6_not_opened',state['retained_closed_guards']['B6']=='NOT-OPENED')
ok('3D5D_pending',state['Selling_3D_to_5D']=='CONTEXT/PENDING-MORPHISM')
ok('checkpoint_gate_closed',CP['state']=='CLOSED/APPEND-ONLY')
ok('version_not_closed',CP['scope'].startswith('RESEARCH-GATE-CLOSURE-ONLY'))
ok('preindustrial_not_publishable','NONPUBLISHABLE' in CP['publication_state'])
ok('closure_token',CP['closure_token']=='V54_GATE6_CLOSURE_PASS')

report={
 'version':'v54','gate':'Gate 6','status':'PASS','gate_closed':True,'v54_closed':False,
 'publication_state':CP['publication_state'],
 'inherited_check_count':base['checks']+hyp['checks']+p46['checks'],
 'closure_checks':len(checks),
 'closure_checks_passed':sum(v for _,v in checks),
 'combined_evidence_checks':base['checks']+hyp['checks']+p46['checks']+len(checks),
 'failed':[n for n,v in checks if not v],
 'closure_token':CP['closure_token'],
 'nonblocking_open_inventory':CP['nonblocking_open_inventory'],
 'next_frontier':CP['next_frontier']['label']
}
(ROOT/'v54_GATE6_CLOSURE_REPORT.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2))
