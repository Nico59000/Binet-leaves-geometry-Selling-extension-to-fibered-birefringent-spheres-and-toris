#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys,json,fitz,hashlib
B=Path(__file__).resolve().parent
checks=[]
def ck(name,cond,detail=None):
 checks.append({'name':name,'pass':bool(cond),'detail':detail})
 if not cond: raise AssertionError(name)
def J(n):return json.load(open(B/n))
# predecessor replay in the current extracted/staged context
r0=subprocess.run([sys.executable,str(B/'verify_field_fox_projective_v50.py')],cwd=B,capture_output=True,text=True,timeout=180)
ck('v50_verifier_returncode',r0.returncode==0,r0.stderr[-500:])
v50=J('v50_VERIFICATION_REPORT.json');ck('v50_verification_pass',v50['status']=='PASS' and v50['combined_checks']==8438)
# v51 mathematical gate
r1=subprocess.run([sys.executable,str(B/'verify_v51_frontier_checkpoint.py')],cwd=B,capture_output=True,text=True,timeout=180)
ck('v51_frontier_verifier_returncode',r1.returncode==0,r1.stderr[-500:])
fv=J('v51_frontier_checkpoint_verification.json');ck('v51_math_55_of_55',fv['status']=='PASS' and fv['checks_passed']==fv['checks_total']==55)
# documentary / QA direct checks
tex=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v51.tex').read_text()
ck('source_has_v51_section','\\section{v51 :' in tex)
ck('source_retains_v50_section','\\section{v50 :' in tex)
ck('source_v51_before_certification',tex.index('\\section{v51 :')<tex.index('\\section{Certification HT+NT du présent formalisme}'))
ck('source_has_v52_frontier','frontière sûre v52' in tex)
ck('source_has_raster_NT','NT/PENDING-DIRECT' in tex)
ck('source_has_global_B5_NT','relèvement intégral global dans }B_5' in tex and r'\texttt{NT}' in tex)
pdf=B/'formalisation_feuillets_binet_birefringence_selling_polaire_v51.pdf';doc=fitz.open(pdf)
ck('pdf_267_pages',len(doc)==267)
qa=J('v51_PDF_QA_REPORT.json');ck('pdf_qa_pass',qa['status']=='PASS');ck('pdf_zero_overfull',qa['latex_overfull_count']==0);ck('pdf_all_page_parity',qa['precanon_canonical_all_page_pixel_parity']['identical_pages']==267);ck('pdf_render_scan_no_flags',qa['automated_render_scan']['flags']==[])
so=J('v51_SECTION_ORDER_AUDIT.json');ck('section_order_pass',so['status']=='PASS' and so['bibliography_physically_final'])
sc=J('v51_SOURCE_CONTEXT.json');ck('source_context_pass',sc['status']=='PASS' and sc['formal_context_document_count']==9 and sc['predecessor_artifacts']['formalisation_feuillets_binet_birefringence_selling_polaire_v50.pdf']['sha256']=='cfcf7d8692c30a39920f05592de1572c38f00c15f4c02589f5a20cef6cd3fa93')
fr=J('v51_formal_rigor_audit.json');ck('formal_rigor_pass',fr['status']=='PASS' and all(fr['checks'].values()))
mc=J('v51_MATH_CHECKPOINT.json');ck('math_gate_pass',mc['math_gate']=='PASS-WITH-EXPLICIT-NT-FRONTIERS')
ck('raster_binding_not_promoted',mc['selling_source']['raster_binding'].startswith('NT/'))
ck('gamma_global_lift_not_promoted',mc['gamma3_radius5']['global_integral_lift_in_B5']=='NT')
ck('R2PURE17_not_merged',mc['retained_separations']['R2PURE13_to_SellingDeck17']=='SEPARATED/PENDING-MORPHISM')
# count: inherited 8438 + frontier 55 + direct checks excluding the two meta calls already counted as direct too
new_direct=len(checks)-4 # first 4 are meta checks establishing inherited/frontier counts
new_total=55+new_direct
combined=8438+new_total
report={'version':'v51','status':'PASS','predecessor_checks':8438,'frontier_math_checks':55,'industrial_direct_checks':new_direct,'new_checks':new_total,'new_passed':new_total,'combined_checks':combined,'failed':[],'direct_checks':checks,'predecessor_replay':{'status':'PASS','stdout_tail':r0.stdout[-1200:]},'frontier_replay':{'status':'PASS','stdout_tail':r1.stdout[-1200:]}}
p=B/'v51_VERIFICATION_REPORT.json';p.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:report[k] for k in ['status','predecessor_checks','frontier_math_checks','industrial_direct_checks','new_checks','combined_checks']},indent=2))
print('sha256',hashlib.sha256(p.read_bytes()).hexdigest())
