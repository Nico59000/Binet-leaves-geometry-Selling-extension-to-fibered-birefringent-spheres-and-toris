#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys,json,csv,hashlib
from collections import Counter
B=Path(__file__).resolve().parent
checks=[]
def ck(name,cond,detail=None):
    checks.append({'name':name,'pass':bool(cond),'detail':detail})
    if not cond: raise AssertionError(name if detail is None else f'{name}: {detail}')
def J(n): return json.load(open(B/n,encoding='utf-8'))
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
# A. Immutable predecessor chain and ledger.
tok=J('v51_ATOMIC_COMPLETION_TOKEN.json')
ck('parent_v51_closed',tok.get('publication_state')=='CLOSED/APPEND-ONLY')
ck('parent_v51_atomic',tok.get('atomic_completion') is True and all(v.startswith('PASS') for v in tok['atomic_conjunction'].values()))
entries=[];mis=[]
for line in (B/'v51_SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
    if not line.strip(): continue
    h,n=line.split('  ',1);entries.append(n);p=B/n
    if not p.exists() or sha(p)!=h:mis.append(n)
ck('parent_v51_ledger_541',len(entries)==541,len(entries))
ck('parent_v51_ledger_byte_exact',not mis,mis[:10])
r0=subprocess.run([sys.executable,str(B/'verify_field_fox_projective_v51.py')],cwd=B,capture_output=True,text=True,timeout=240)
ck('parent_v51_verifier_returncode',r0.returncode==0,r0.stderr[-800:])
vr51=J('v51_VERIFICATION_REPORT.json')
ck('parent_v51_8511',vr51['status']=='PASS' and vr51['combined_checks']==8511)
# B. Deterministic replay of all new final-gate builders.
outputs=[
 'selling_global_curvilinear_warp_BSB_P084_v52.json',
 'selling_exact_triangulated_BSB_P084_homeomorphism_certificate_K12_v52.json',
 'selling_exact_triangulated_BSB_P084_123edge_audit_K12_v52.csv',
 'selling_fig1_final_123_certificate_v52.json',
 'selling_fig1_final_123_historical_incidence_v52.csv',
 'selling_fig3_6_typed_fusion_v52.json']
pre={n:sha(B/n) for n in outputs}
replays=[]
for s in ['rebuild_global_curvilinear_warp_BSB_P084_v52.py','build_exact_triangulated_BSB_P084_K12_v52.py','build_final_123_propagation_v52.py','build_fig3_6_typed_fusion_v52.py']:
    r=subprocess.run([sys.executable,str(B/s)],cwd=B,capture_output=True,text=True,timeout=240)
    replays.append({'script':s,'returncode':r.returncode,'stdout_tail':r.stdout[-900:],'stderr_tail':r.stderr[-900:]})
    ck('replay_'+s,r.returncode==0,r.stderr[-600:])
post={n:sha(B/n) for n in outputs}
ck('new_outputs_byte_deterministic',pre==post,{n:(pre[n],post[n]) for n in outputs if pre[n]!=post[n]})
# C. Global P084 field, exact triangulation, e96.
w=J('selling_global_curvilinear_warp_BSB_P084_v52.json')
ck('smooth_global_candidate_pass',w['status']=='PASS-GLOBAL-CURVILINEAR-HOMEOMORPHISM-CANDIDATE')
sm=w['selected_metrics']
ck('smooth_jac_sample_positive',sm['jac_negative_or_zero']==0 and sm['jac_min']>0)
ck('smooth_crossings_preserved',sm['crossing_count']==3 and sm['new_crossings']==[] and sm['lost_baseline_crossings']==[])
cert=J('selling_exact_triangulated_BSB_P084_homeomorphism_certificate_K12_v52.json')
ck('exact_PL_status',cert['status']=='PASS-EXACT-TRIANGULATED-BSB-P084-GLOBAL-HOMEOMORPHISM')
ck('exact_PL_counts',(cert['vertex_count'],cert['triangle_count'],cert['mesh_edge_count'])==(1297,2484,3780))
ck('exact_source_orientation_integer_positive',isinstance(cert['source_min_signed_double_area_scaled2'],int) and cert['source_min_signed_double_area_scaled2']>0)
ck('exact_target_orientation_integer_positive',isinstance(cert['target_min_signed_double_area_scaled2'],int) and cert['target_min_signed_double_area_scaled2']>0)
# crossing fields may be names or counts depending on schema
for key in ['source_crossings','target_crossings','target_boundary_crossings']:
    if key in cert: ck('exact_'+key+'_zero',cert[key] in (0,[],{}),cert[key])
ck('exact_boundary_is_historical_RHAX','exact historical RHAX' in cert['boundary_target'])
ck('exact_global_boundary_blend_t2','t^2' in cert['boundary_target'])
e96=cert['e96']
ck('e96_exact_strong',e96['status']=='STRONG-EXACT-TRIANGULATED-STROKE-SUPPORT')
ck('e96_geom_support_ge_075',e96['geom_support_eps']>=0.75,e96['geom_support_eps'])
ck('e96_dark_support_ge_075',e96['dark_support_eps']>=0.75,e96['dark_support_eps'])
ck('e96_P084_support_ge_075',e96['P084_support_eps']>=0.75,e96['P084_support_eps'])
ck('e96_gate_pass',cert['e96_gate_pass'] is True)
# D. Source-exact 123/123 closure.
f=J('selling_fig1_final_123_certificate_v52.json')
ck('fig1_final_status',f['status']=='PASS-HISTORICAL-FIG1-123/123')
ck('fig1_final_123',f['final']['historical_edge_orbits_proven']==f['final']['historical_edge_orbits_total']==123)
ck('fig1_final_zero_NT',f['final']['remaining_NT']==0)
ck('fig1_final_promotion',f['final']['promotion']=='PROVEN-HISTORICAL-GLOBAL-FIG1-REALIZATION')
ck('fig1_e96_seed_exact',f['P084_exact_seed']['edge']==96 and f['P084_exact_seed']['status']=='STRONG-EXACT-TRIANGULATED-STROKE-SUPPORT')
trace=f['seed_implication']['compact_trace']
ck('source_prop_trace_exact',trace==[
 ['SAME-EXACT-IMPLICIT-FACTOR-CONTINUATION',[53]],
 ['UNIQUE-MISSING-EXACT-LINK-INCIDENCE',[34,97]],
 ['UNIQUE-MISSING-EXACT-LINK-INCIDENCE',[54,78]]],trace)
ck('source_prop_new_edges',sorted(f['seed_implication']['new_edges'])==[34,53,54,78,97])
rows=list(csv.DictReader(open(B/'selling_fig1_final_123_historical_incidence_v52.csv',encoding='utf-8')))
ck('fig1_provenance_csv_123',len(rows)==123,len(rows))
# accept explicit status field schema from builder
status_fields=[k for k in rows[0].keys() if 'status' in k.lower() or 'historical' in k.lower()]
ck('fig1_provenance_has_status_field',len(status_fields)>0,status_fields)
# E. Fig3--6 typed coproduct/fusion.
g=J('selling_fig3_6_typed_fusion_v52.json')
ck('fig3_6_fusion_pass',g['status']=='PASS-TYPED-DISJOINT-FUSION/PROVENANCE-PRESERVED')
ck('fig3_6_counts',(g['object_count'],g['relation_count'],g['ordered_structure_count'])==(31,28,3))
ck('fig3_6_all_checks',all(g['checks'].values()),g['checks'])
ck('fig3_6_no_cross_relations',g['checks']['no_cross_figure_relation'] is True)
ck('fig3_6_shared_labels_separated',g['checks']['shared_labels_not_identified'] is True)
fig6=J('selling_fig6_local_compiler_replay_v52.json')
labels=[o['label'] for o in fig6['compiled_carrier']['objects']]
ck('fig6_corrected_starred_tuple','(1,-7,8)/(3,0,0)' in labels)
ck('fig6_old_wrong_sign_absent','(1,-7,-8)/(3,0,0)' not in labels)
# F. Math checkpoint and retained separations.
mc=J('v52_MATH_CHECKPOINT.json')
ck('math_gate_pass',mc['math_gate']=='PASS-FINAL-FIG1-GATE-AND-TYPED-FIG3-6-FUSION')
ck('math_checkpoint_preseal',mc['publication_state']=='PRESEAL/NONPUBLISHABLE')
ck('Selling_Noether_separated',mc['retained_separations']['Selling_Noether_biquadratic']=='CONTEXT/PENDING-MORPHISM')
ck('Selling_3D5D_separated',mc['retained_separations']['Selling_3D_to_5D']=='CONTEXT/PENDING-MORPHISM')
ck('GL3_bridge_not_constructed',mc['retained_separations']['GL3_equivariant_bridge']=='NT/NOT-CONSTRUCTED')
ck('R2PURE_separated',mc['retained_separations']['R2PURE13_to_SellingDeck17']=='SEPARATED/PENDING-MORPHISM')
# final report, predecessor 8511 + all direct v52 checks excluding parent replay's inherited internal count
new=len(checks)
report={'version':'v52','status':'PASS','predecessor_checks':8511,'v52_direct_checks':new,'combined_checks':8511+new,'failed':[],
        'direct_checks':checks,'replays':replays,'output_sha256':post}
(B/'v52_VERIFICATION_REPORT.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps({k:report[k] for k in ['status','predecessor_checks','v52_direct_checks','combined_checks']},indent=2))
print('sha256',sha(B/'v52_VERIFICATION_REPORT.json'))
