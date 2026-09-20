#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,subprocess,sys,csv
import fitz
B=Path(__file__).resolve().parent
checks={}
def ok(k,v):
 checks[k]=bool(v)
 if not v: print('FAIL',k,file=sys.stderr)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
# Parent v47 closed and byte-exact.
tok=json.loads((B/'v47_ATOMIC_COMPLETION_TOKEN.json').read_text())
ok('parent_v47_closed',tok['publication_state']=='CLOSED/APPEND-ONLY')
ok('parent_v47_atomic',tok['atomic_completion'] is True and all(v=='PASS' for v in tok['atomic_conjunction'].values()))
ok('parent_v47_checks_8114',tok['verification']['combined']==8114)
entries=[];mis=[]
for line in (B/'v47_SHA256SUMS.txt').read_text().splitlines():
 if not line.strip():continue
 h,n=line.split('  ',1);entries.append(n);p=B/n
 if not p.exists() or sha(p)!=h:mis.append(n)
ok('parent_v47_ledger_305',len(entries)==305)
ok('parent_v47_ledger_byte_exact',not mis)
ok('parent_v47_pdf_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v47.pdf')==tok['pdf']['sha256'])
ok('parent_v47_tex_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v47.tex')==tok['tex']['sha256'])
# Editorial structural repair.
audit=json.loads((B/'v48_SECTION_ORDER_AUDIT.json').read_text())
ok('order_audit_pass',audit['status']=='PASS')
ok('order_v46_v47_v48_before_cert',audit['v46_v47_v48_before_certification'])
ok('order_appendix_after_conclusion',audit['appendix_after_conclusion'])
ok('order_bibliography_final_outline',audit['bibliography_physically_final_outline'])
ok('order_v43_name_repaired',audit['v43_naming_repaired'])
ok('order_references_numbered_removed',audit['references_numbered_section_removed'])
ok('order_references_unnumbered',audit['references_unnumbered_via_refname'])
ok('order_pages_254',audit['pdf_pages']==254)
tex=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v48.tex').read_text()
ix={k:tex.index(k) for k in ['\\section{v46 :','\\section{v47 :','\\section{v48 :','\\section{Certification HT+NT','\\section{Conclusion}','\n\\appendix\n','\\begin{thebibliography}']}
ok('tex_order_monotone',list(ix.values())==sorted(ix.values()))
ok('tex_v43_prefix_style','\\section{v43 :' in tex)
ok('tex_bibliography_after_appendix',ix['\\begin{thebibliography}']>ix['\n\\appendix\n'])
ok('tex_bibliography_before_end',ix['\\begin{thebibliography}']<tex.rindex('\\end{document}'))
# Selling centre transcription and fail-closed gates.
cent=json.loads((B/'selling_fourfold_centres_transcription_v48.json').read_text())
curve=json.loads((B/'selling_axis_curve_transcription_v48.json').read_text())
c2=json.loads((B/'selling_degree1_C2_gate_v48.json').read_text())
h12=json.loads((B/'selling_h12_materialization_gate_v48.json').read_text())
ok('sell48_source_exact_six',cent['source_exact_constraints'][0].startswith('exactly six'))
ok('sell48_diag_nodes_136',cent['diagnostic_node_table']['total_nodes']==136)
ok('sell48_solid_degree4_count2',cent['diagnostic_node_table']['solid_degree4_count']==2)
ok('sell48_solid_degree4_ids',{x['node_id'] for x in cent['diagnostic_node_table']['solid_degree4_nodes']}=={'N022','N109'})
ok('sell48_no_promoted_HC',cent['promoted_historical_centres']==[])
ok('sell48_HC_nt',cent['status'].startswith('NT/'))
ok('sell48_no_HAX_curves',curve['promoted_HAX_curves']==[])
ok('sell48_no_reflections',curve['promoted_reflections']==[])
ok('sell48_letters_unbound',curve['vanishing_variable_letters_bound_to_curves'] is False)
ok('sell48_PQRS_unbound',curve['PQRS_raster_routes_bound'] is False)
ok('sell48_C2_89',c2['degree_one_carrier_size']==89)
ok('sell48_C2_no_reflection',c2['historical_reflection_count_available']==0)
ok('sell48_C2_not_formable',c2['carrier_invariance_test']=='NOT-FORMABLE')
ok('sell48_C2_fixed_none',c2['fixed_points'] is None and c2['free_orbits'] is None and c2['total_orbits'] is None)
ok('sell48_44_45_not_tested',c2['comparison_44_45']=='SEPARATED_NOT_TESTED')
ok('sell48_34_55_not_tested',c2['comparison_34_55']=='SEPARATED_NOT_TESTED')
ok('sell48_h12_sizes',h12['diagnostic_cells']=={'C0_nodes':136,'C1_segments':105,'C2_fields':34})
ok('sell48_h12_d0_not_formed',h12['d0']=='NOT-FORMED')
ok('sell48_h12_d1_not_formed',h12['d1']=='NOT-FORMED')
ok('sell48_h12_not_recomputed',h12['historical_H1'].startswith('NT_') and h12['historical_H2'].startswith('NT_'))
# Gamma radius-3 exact bounded no-go.
p=subprocess.run([sys.executable,str(B/'verify_gamma3_radius3_mod4_v48.py')],cwd=B,capture_output=True,text=True)
ok('g48_dual_helper_exit_zero',p.returncode==0)
no=json.loads((B/'gamma3_character_radius3_no_go_v48.json').read_text())
cert=json.loads((B/'gamma3_radius3_mod4_dual_certificate_v48.json').read_text())
part=json.loads((B/'gamma3_partial_resolution_v48.json').read_text())
ok('g48_target_chires14',no['target']=='CHIRES_14 = R3b_123')
ok('g48_B3_size1376',no['exact_group_ball_radius3']['size']==1376)
ok('g48_B3_profile',no['exact_group_ball_radius3']['depth_profile']=={'0':1,'1':15,'2':147,'3':1213})
ok('g48_mod4_group512',no['finite_quotient']['group_size']==512)
ok('g48_B3_mod4_130',no['radius3_image_in_mod4']['distinct_group_elements']==130)
ok('g48_mod4_shape',no['bounded_mod4_Fox_matrix']=={'rows':3624,'columns':5590})
ok('g48_dual_terms127',cert['dual_certificate_nonzero_rows']==127 and len(cert['dual_terms'])==127)
ok('g48_dual_den1',cert['denominator']=='1')
ok('g48_dual_identity',cert['identity']=='y^T D_mod4 = l_CHIRES14 exactly over Q')
ok('g48_radius3_refuted','REFUTED-TYPED' in no['radius3_lift_status'])
ok('g48_min_radius_ge4',no['minimal_radius_if_lift_exists']=='>=4')
ok('g48_13_retained',len(no['retained_radius2_lifts'])==13)
ok('g48_partial_before14',part['character_residual_before_v47']==14)
ok('g48_partial_lifted13',part['radius2_lifted_directions']==13)
ok('g48_partial_remaining',part['remaining_character_direction']=='CHIRES_14 = R3b_123')
ok('g48_full_module_nt',part['complete_resolution']=='NT')
# Retain v47 exact cycles artifact and guards.
r2=list(csv.DictReader(open(B/'gamma3_radius2_pure_lifts_v47.csv',encoding='utf-8')))
ok('g48_r2_csv_nonempty',len(r2)>0)
ok('g48_r2_13_ids',len({r['lift_id'] for r in r2})==13)
r5=json.loads((B/'selling_gamma_R5_adapter_v43.json').read_text())
ok('g48_R5_refuted_retained','REFUTED-TYPED' in r5['strict_incidence_preserving_16_to_16_adapter'])
sing=json.loads((B/'selling_singer_injection_gate_v45.json').read_text())
ok('g48_singer_pending',sing['injection_to_J_line'].startswith('NT_'))
ok('g48_lorentz_nt',sing['Lorentz_orientation']=='NT')
reg=json.loads((B/'pending_morphism_registry_v48.json').read_text())
ok('g48_registry_12',len(reg['entries'])==12)
# PDF outline final check from canonical PDF.
pdf=B/'formalisation_feuillets_binet_birefringence_selling_polaire_v48.pdf';d=fitz.open(pdf);toc=d.get_toc(simple=True);lvl1=[(t,p) for l,t,p in toc if l==1]
ok('pdf48_pages254',d.page_count==254)
ok('pdf48_bibliography_last_outline',lvl1[-1][0]=='Références documentaires du corpus')
ok('pdf48_v48_before_cert',next(p for t,p in lvl1 if t.startswith('v48 :'))<next(p for t,p in lvl1 if t.startswith('Certification HT+NT')))
failed=[k for k,v in checks.items() if not v]
new=len(checks);report={'version':'v48','predecessor_checks':8114,'new_checks':new,'combined_checks':8114+new,'new_passed':new-len(failed),'status':'PASS' if not failed else 'FAIL','failed':failed}
(B/'v48_VERIFICATION_REPORT.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
print(json.dumps(report,indent=2,ensure_ascii=False))
sys.exit(0 if not failed else 1)
