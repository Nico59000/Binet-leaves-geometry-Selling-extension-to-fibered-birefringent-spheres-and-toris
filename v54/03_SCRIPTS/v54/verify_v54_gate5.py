#!/usr/bin/env python3
from pathlib import Path
import json, hashlib, sys
ROOT=Path(__file__).resolve().parent
G4=Path('/mnt/data/v54_gate4')
checks={}
def ok(name,cond):
    checks[name]=bool(cond)
def load(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

# Predecessor
p4=load(G4/'v54_GATE4_CHECKPOINT.json')
ok('predecessor_gate4_preseal',p4['publication_state']=='RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE')
ok('predecessor_fig2_closed','PASS-CLOSED-EXACT/8-CLASSES' in p4['results']['Selling_Fig2_source_quotient'])
ok('predecessor_fig3_closed','PASS-CLOSED-EXACT/27-CLASSES' in p4['results']['Selling_Fig3_fissure_source_quotient'])
ok('predecessor_B6_closed',p4['results']['B6']=='NOT-OPENED')
ok('predecessor_3D5D_pending','PENDING-MORPHISM' in p4['results']['Selling_3D_to_5D'])

# Numeric typed bridge
b=load(ROOT/'selling_fig2_fig3_numeric_typed_bridge_v54.json')
ok('bridge_status',b['status'].startswith('PASS-SOURCE-EXACT'))
ok('bridge_fig2_classes',b['counts']['Fig2_classes']==8)
ok('bridge_fig3_classes',b['counts']['Fig3_classes']==27)
ok('bridge_fig2_edges',b['counts']['Fig2_edges']==26)
ok('bridge_fig3_edges',b['counts']['Fig3_edges']==91)
ok('bridge_total_edges',b['counts']['total_quotient_edges']==117)
ok('bridge_source_117',b['counts']['all_source_certified_edges']==117)
ok('bridge_raster_zero',b['counts']['direct_raster_bound_edges']==0)
ok('bridge_fig2_root',b['Fig2']['root_class']==7)
ok('bridge_fig3_root',b['Fig3']['root_class']==19)
ok('bridge_fig2_root_anchor',any(x['class']==7 and x['object']=='SellCommon:Fig2:C3' and x['existing_v53_object'] for x in b['Fig2']['class_bindings']))
ok('bridge_fig3_root_new',any(x['class']==19 and not x['existing_v53_object'] for x in b['Fig3']['class_bindings']))
ok('bridge_all_edges_source_exact',all(e['source_incidence']=='PROVEN-EXACT-HISTORICAL-WALL-ADJACENCY' for F in ['Fig2','Fig3'] for e in b[F]['edge_certificates']))
ok('bridge_no_cross_figure',all(e['typed_source_object'].startswith('SellCommon:'+F+':') and e['typed_target_object'].startswith('SellCommon:'+F+':') for F in ['Fig2','Fig3'] for e in b[F]['edge_certificates']))

# Raster-binding obstruction
r=load(ROOT/'selling_fig2_fig3_raster_binding_obstruction_v54.json')
ok('raster_total_117','117/117' in r['source_edge_certification']['total'])
ok('fig2_mapA_preserves',r['Fig2']['map_A_audit']['preserved'] and r['Fig2']['map_A_audit']['checked_edges']==15)
ok('fig2_mapB_preserves',r['Fig2']['map_B_audit']['preserved'] and r['Fig2']['map_B_audit']['checked_edges']==15)
ok('fig2_maps_distinct',r['Fig2']['maps_distinct'])
ok('fig2_unique_binding_refuted','REFUTED-TYPED' in r['Fig2']['decision'])
ok('fig3_root_absent',not r['Fig3']['fissure_tuple_present_in_v53_raster_labels'])
ok('fig3_raster_binding_nt',r['Fig3']['decision'].startswith('NT/'))
ok('raster_overall_partial','PARTIAL/NT' in r['overall_status'])

# Fig4-6 exact no-extension gate
f=load(ROOT/'selling_fig4_6_gate5_no_extension_required_v54.json')
ok('fig456_gate_status',f['status'].startswith('PASS-'))
for F in ['Fig4','Fig5','Fig6']:
    x=f['figures'][F]
    ok(F+'_one_field',x['field_count']==1 and x['processed']==1)
    ok(F+'_closed',x['compile_status']=='PASS-CLOSED-EXACT' and x['exact_topology'])
    ok(F+'_no_fallback',x['fallback_count']==0)
    ok(F+'_no_extension',x['additional_repetition_presentation_triggered'] is False)

# Noether bounded source search / locus
n=load(ROOT/'noether_historical_quartic_search_and_selling_locus_v54.json')
ok('noether_source_hash',n['source']['sha256']=='ad9da4266cf035a25c3ca0a27cf87a9c06d9fa980edacae2994b725237532968')
ok('noether_two_controls',len(n['bounded_search']['source_audited_explicit_full_quartics'])==2)
ok('noether_controls_rho_zero',all(x['rho']=='0' for x in n['bounded_search']['source_audited_explicit_full_quartics']))
ok('noether_third_not_found',n['bounded_search']['third_distinct_explicit_full_quartic_in_current_source'].startswith('NOT-FOUND'))
ok('noether_scope_guard','not a no-go' in n['bounded_search']['scope_guard'].lower())
ok('noether_Zrho_proven',n['rho_carriers']['Z_rho'].startswith('PROVEN'))
ok('noether_Urho_proven',n['rho_carriers']['U_rho'].startswith('PROVEN'))
ok('noether_term_profile',n['rho_carriers']['rho_term_count_profile']==[66,66,66,68,68,68])
ok('noether_signature_locus','signature' in n['real_selling_admissible_locus']['equivalent_signature_statement'].lower())
ok('noether_historical_still_nt',n['real_selling_admissible_locus']['historical_point']=='OPEN/NT')

# Current Peiffer/higher presentation gate
p=load(ROOT/'r2pure13_existing_peiffer_shadow_gate_v54.json')
ok('peiffer_248',p['existing_families']['v44_triple_Peiffer']['comparison_cycles']==248)
ok('peiffer_248_absorbed',p['existing_families']['v44_triple_Peiffer']['absorbed_by_M43']==248)
ok('peiffer_no_new_shadow',p['existing_families']['v44_triple_Peiffer']['new_group_ring_direction_mod_M43']==0)
ok('peiffer_second_18',p['existing_families']['v43_second_syzygies']['count']==18)
ok('r2pure_shadow_13',p['R2PURE_shadow']['independent_residual_directions']==13 and len(p['R2PURE_shadow']['directions'])==13)
ok('peiffer_current_family_refuted','REFUTED-TYPED' in p['decision'])
ok('peiffer_new_family_open','OPEN/NT' in p['open_lane'])
ok('peiffer_no_radius4_open','Radius 4 is not opened' in p['open_lane'])

failed=[k for k,v in checks.items() if not v]
rep={'version':'v54','phase':'Gate5','checks':checks,'check_count':len(checks),'passed':len(checks)-len(failed),'failed':failed,'status':'PASS' if not failed else 'FAIL'}
(ROOT/'v54_GATE5_VERIFICATION_REPORT.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n")
print(json.dumps({'status':rep['status'],'passed':rep['passed'],'check_count':rep['check_count'],'failed':failed},indent=2))
sys.exit(0 if not failed else 1)
