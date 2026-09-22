#!/usr/bin/env python3
from pathlib import Path
import subprocess, tempfile, json, hashlib, sys, os
B=Path(__file__).resolve().parent
checks=[]
def ck(name,cond,detail=None):
    checks.append({'name':name,'pass':bool(cond),'detail':detail})
    if not cond: raise AssertionError(f'{name}: {detail}')
def J(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
# Parent sealed state and append-only audit.
tok=J(B/'v52_ATOMIC_COMPLETION_TOKEN.json')
ck('parent_v52_closed',tok['atomic_completion'] and tok['publication_state']=='CLOSED/APPEND-ONLY')
ck('parent_v52_bundle_hash',tok['bundle']['sha256']=='1936700dc187ae762a407a17f16019c0f871ff47962a5e01df0eed577c38787d')
ck('parent_v52_verify_8583',tok['verification']['combined']==8583)
ap=J(B/'v53_APPEND_ONLY_AUDIT.json')
ck('parent_594_byte_exact',ap['status']=='PASS' and ap['checked']==594 and ap['append_only_exact'])
# Exact CHIRES replay from raw mod-8 classes/Fox terms and raw dual solution.
with tempfile.TemporaryDirectory(prefix='v53chi_') as td:
    exe=Path(td)/'verify_chi'; out=Path(td)/'chi.json'
    c=subprocess.run(['g++','-O3','-std=c++17',str(B/'verify_chires14_b5_global_dual_v53.cpp'),'-o',str(exe)],capture_output=True,text=True)
    ck('chi_compile',c.returncode==0,c.stderr[-1000:])
    r=subprocess.run([str(exe),str(B/'06_CHIRES/v53/data/gamma3_b5_mod8_packed.bin'),str(B/'06_CHIRES/v53/data/gamma3_fox_mod8_terms.txt'),str(B/'06_CHIRES/v53/data/gamma3_radius5_mod8_raw_dual_solution_v53.bin'),str(out)],capture_output=True,text=True)
    ck('chi_replay_returncode',r.returncode==0,r.stderr[-1200:])
    chi=J(out)
ck('chi_raw_956320',chi['raw_equations']==956320 and chi['raw_variables']==747196 and chi['raw_incidence']==4403520)
ck('chi_raw_zero_fail',chi['raw_verification_failures']==0 and chi['target_verification_failures']==0 and chi['nontarget_verification_failures']==0)
ck('chi_core_counts',chi['degree_one_peeling']['core_equations']==856952 and chi['degree_one_peeling']['core_variables']==465651 and chi['degree_one_peeling']['core_incidence']==3528672)
ck('chi_core_fixed_point',chi['degree_one_peeling']['min_active_degree']==2 and chi['degree_one_peeling']['peeled_equations']==99368 and chi['degree_one_peeling']['core_rhs_one']==44)
ck('chi_solution_weight',chi['raw_solution_weight']==121219)
ck('chi_solution_hash',sha(B/'06_CHIRES/v53/data/gamma3_radius5_mod8_raw_dual_solution_v53.bin')=='1674a69d25c349e036701a9a01c057ca9bc3ea7b72e916dac03f62f773f48bf6')
ck('chi_archived_core_hash',sha(B/'06_CHIRES/v53/data/gamma3_radius5_mod8_core.bin')=='260ce4b5e0c09f0db2d338bc8d1173f8231b695d957a608eb63fe70c643c84c4')
cc=J(B/'gamma3_radius5_global_mod8_dual_obstruction_v53.json')
ck('chi_status_refuted_typed',cc['status'].startswith('PROVEN-GLOBAL-MOD8-DUAL-OBSTRUCTION') and cc['theorem']['chi_res_14_radius_le_5']=='REFUTED-TYPED')
ck('radius6_not_opened',cc['theorem']['radius6']=='NOT-OPENED')
# Selling common typed functor replay.
r=subprocess.run([sys.executable,str(B/'verify_selling_fig2_6_common_functor_v53.py')],cwd=B,capture_output=True,text=True)
ck('selling_functor_replay',r.returncode==0,r.stderr[-1000:])
sfr=J(B/'v53_SELLING_FUNCTOR_REPLAY_REPORT.json')
ck('selling_functor_24',sfr['status']=='PASS' and sfr['passed']==sfr['total']==24)
f=J(B/'selling_fig2_6_common_typed_functor_v53.json')
ck('selling_counts',f['object_count']==42 and f['morphism_count']==43 and f['ordered_structure_count']==3)
ck('selling_no_cross_merge',f['checks']['no_cross_figure_morphisms'] and f['checks']['all_provenance_present'])
# Append-only TeX insertion.
v52=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v52.tex').read_text(encoding='utf-8')
v53=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v53.tex').read_text(encoding='utf-8')
ins=(B/'v53_SECTION_INSERT.tex').read_text(encoding='utf-8').rstrip()+'\n\n'
ck('tex_has_v53','\\section{v53 :' in v53)
ck('tex_v53_before_certification',v53.index('\\section{v53 :')<v53.index('\\section{Certification HT+NT du présent formalisme}'))
ck('tex_append_only_exact',v53.replace(ins,'',1)==v52)
ck('tex_safe_v54','frontière sûre v54' in v53)
# Math/context/formal rigor.
mc=J(B/'v53_MATH_CHECKPOINT.json'); ck('math_gate',mc['math_gate']=='PASS-CHIRESON-B5-AND-FIG2-6-COMMON-TYPED-FUNCTOR')
ck('math_retained_separations',mc['retained_separations']['R2PURE13_to_SellingDeck17']=='SEPARATED/PENDING-MORPHISM' and mc['retained_separations']['Selling_Noether']=='CONTEXT/PENDING-MORPHISM' and mc['retained_separations']['Selling_3D_to_5D']=='CONTEXT/PENDING-MORPHISM')
fc=J(B/'v53_FORMAL_CONTEXT_ATTESTATION.json'); ck('formal_context',fc['status'].startswith('PASS-FORMAL-CONTEXT-LOCK') and fc['document_count']==9)
ra=J(B/'v53_formal_rigor_audit.json'); ck('formal_rigor',ra['status']=='PASS' and all(ra['checks'].values()))
# PDF QA/document order.
qa=J(B/'v53_PDF_QA_REPORT.json')
ck('pdf_qa_pass',qa['status']=='PASS' and qa['pdf_pages']==274 and qa['latex_overfull_count']==0)
ck('pdf_all_page_parity',qa['precanon_canonical_all_page_pixel_parity']['identical_pages']==274)
ck('pdf_scan_no_flags',qa['automated_render_scan']['flags']==[] and qa['pdfium_renderer']=='PASS_274_PAGES_AT_72_DPI')
ck('pdf_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v53.pdf')==qa['pdf_sha256'])
so=J(B/'v53_SECTION_ORDER_AUDIT.json')
ck('section_order',so['status']=='PASS' and so['v53_before_certification'] and so['bibliography_physically_final'])
ck('v53_printed_page_263',263 in so['pdf_pages_by_heading']['v53'])
# Bundle tree navigation.
bt=J(B/'00_INDEX/BUNDLE_TREE_AUDIT.json')
ck('tree_pass',bt['status']=='PASS' and bt['version_mirrors']>=35 and bt['script_mirrors']>=120 and bt['certificate_mirrors']>=180)
ck('tree_dirs',all((B/d).is_dir() for d in bt['directories']))
ck('tree_under_200MiB_preseal',bt['total_uncompressed_bytes_after_tree']<200*1024*1024)
# Summary.
new=len(checks); total=8583+new
rep={'version':'v53','status':'PASS','predecessor_checks':8583,'industrial_and_frontier_direct_checks':new,'combined_checks':total,'failed':[],'direct_checks':checks,'chi_replay':chi,'selling_functor_replay':sfr,'pdf_sha256':sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v53.pdf'),'tex_sha256':sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v53.tex')}
(B/'v53_VERIFICATION_REPORT.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps({k:rep[k] for k in ['status','predecessor_checks','industrial_and_frontier_direct_checks','combined_checks']},indent=2))
print('sha256',sha(B/'v53_VERIFICATION_REPORT.json'))
