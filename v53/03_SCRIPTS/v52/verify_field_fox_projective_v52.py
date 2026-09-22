#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys,json,fitz,hashlib
B=Path(__file__).resolve().parent
checks=[]
def ck(name,cond,detail=None):
    checks.append({'name':name,'pass':bool(cond),'detail':detail})
    if not cond: raise AssertionError(name if detail is None else f'{name}: {detail}')
def J(n):return json.load(open(B/n,encoding='utf-8'))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
# Full final mathematical/source gate, which itself replays v51 and all v52 builders.
r0=subprocess.run([sys.executable,str(B/'verify_v52_final_gate.py')],cwd=B,capture_output=True,text=True,timeout=600)
ck('v52_final_gate_returncode',r0.returncode==0,r0.stderr[-1000:])
fr=J('v52_VERIFICATION_REPORT.json')
ck('v52_final_gate_8561',fr['status']=='PASS' and fr['combined_checks']==8561,fr.get('combined_checks'))
# Append-only source insertion.
v51=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v51.tex').read_text(encoding='utf-8')
v52=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v52.tex').read_text(encoding='utf-8')
ins=(B/'v52_SECTION_INSERT.tex').read_text(encoding='utf-8')
ck('source_has_v52_section','\\section{v52 :' in v52)
ck('source_retains_v51_section','\\section{v51 :' in v52)
ck('source_v52_before_certification',v52.index('\\section{v52 :')<v52.index('\\section{Certification HT+NT du présent formalisme}'))
ck('append_only_insertion_byte_exact_text',v52.replace(ins+'\n','',1)==v51)
ck('source_has_final_123','123/123' in v52 and '70\\to0' in v52)
ck('source_has_v53_frontier','frontière sûre v53' in v52)
ck('source_retains_GL3_guard','GL(3)' in v52 and 'PENDING-MORPHISM' in v52)
# PDF and QA.
pdf=B/'formalisation_feuillets_binet_birefringence_selling_polaire_v52.pdf';doc=fitz.open(pdf)
ck('pdf_271_pages',len(doc)==271,len(doc))
qa=J('v52_PDF_QA_REPORT.json')
ck('pdf_qa_pass',qa['status']=='PASS')
ck('pdf_zero_overfull',qa['latex_overfull_count']==0)
ck('pdf_all_page_parity',qa['precanon_canonical_all_page_pixel_parity']['identical_pages']==271)
ck('pdf_render_scan_no_flags',qa['automated_render_scan']['flags']==[])
ck('pdf_hash_matches_qa',sha(pdf)==qa['pdf_sha256'])
so=J('v52_SECTION_ORDER_AUDIT.json')
ck('section_order_pass',so['status']=='PASS' and so['v52_before_certification'] and so['bibliography_physically_final'])
ck('v52_printed_page_260',260 in so['pdf_pages_by_heading']['v52'],so['pdf_pages_by_heading']['v52'])
# Context and formal rigor.
sc=J('v52_SOURCE_CONTEXT.json')
ck('source_context_pass',sc['status']=='PASS' and sc['formal_context_document_count']==9)
fa=J('v52_FORMAL_CONTEXT_ATTESTATION.json')
ck('formal_context_pass',fa['status'].startswith('PASS-FORMAL-CONTEXT-LOCK'))
ra=J('v52_formal_rigor_audit.json')
ck('formal_rigor_pass',ra['status']=='PASS' and all(ra['checks'].values()))
mc=J('v52_MATH_CHECKPOINT.json')
ck('math_gate_pass',mc['math_gate']=='PASS-FINAL-FIG1-GATE-AND-TYPED-FIG3-6-FUSION')
ck('cross_carriers_retained',mc['retained_separations']['Selling_Noether_biquadratic']=='CONTEXT/PENDING-MORPHISM' and mc['retained_separations']['Selling_3D_to_5D']=='CONTEXT/PENDING-MORPHISM')
# count cumulative
new=len(checks);combined=8561+new
report={'version':'v52','status':'PASS','predecessor_and_final_gate_checks':8561,'industrial_direct_checks':new,'new_checks_after_final_gate':new,'combined_checks':combined,'failed':[],
        'direct_checks':checks,'final_gate_replay':{'status':'PASS','stdout_tail':r0.stdout[-1600:]},
        'pdf_sha256':sha(pdf),'tex_sha256':sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v52.tex')}
(B/'v52_INDUSTRIAL_VERIFICATION_REPORT.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps({k:report[k] for k in ['status','predecessor_and_final_gate_checks','industrial_direct_checks','combined_checks']},indent=2))
print('sha256',sha(B/'v52_INDUSTRIAL_VERIFICATION_REPORT.json'))
