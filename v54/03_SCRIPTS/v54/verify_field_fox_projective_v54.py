#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,os,sys
B=Path(__file__).resolve().parent
checks=[]
def ck(n,c,d=None):
 checks.append({'name':n,'pass':bool(c),'detail':d})
 if not c: raise AssertionError(f'{n}: {d}')
def J(n):return json.loads((B/n).read_text(encoding='utf-8'))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for z in iter(lambda:f.read(1<<20),b''):h.update(z)
 return h.hexdigest()
# Parent exactness and sealed status.
pm=J('v54_PARENT_V53_1007_SHA256.json');bad=[]
for rel,m in pm['files'].items():
 p=B/rel
 if not p.exists() or p.stat().st_size!=m['bytes'] or sha(p)!=m['sha256']:bad.append(rel)
ck('parent_1007_byte_exact',not bad,bad[:10])
t=J('v53_ATOMIC_COMPLETION_TOKEN.json');ck('v53_closed',t['publication_state']=='CLOSED/APPEND-ONLY')
ck('v53_bundle_hash',t['bundle']['sha256']=='e120360f2bffad1a6eaf869822721c7fb3c515a54558b0ffa5e01a66143df74c')
# Freeze and gates.
f=J('v54_PREINDUSTRIAL_FREEZE_REPORT.json');ck('freeze_24',f['status']=='PASS' and f['passed']==f['direct_checks']==24)
for n,key,want in [('v54_RESEARCH_VERIFICATION_REPORT.json','passed',41),('v54_GATE2_VERIFICATION_REPORT.json','check_count',40),('v54_GATE4_VERIFICATION_REPORT.json','passed',56),('v54_GATE5_VERIFICATION_REPORT.json','passed',58)]:
 d=J(n);ck('gate_'+n,d['status']=='PASS' and d[key]==want)
g3=J('v54_GATE3_VERIFICATION_REPORT.json');ck('gate3_35',g3['status']=='PASS' and g3['new_passed']==35)
g6=J('v54_GATE6_CLOSURE_REPORT.json');ck('gate6_77',g6['status']=='PASS' and g6['closure_checks_passed']==77 and g6['gate_closed'])
# Controlling math statuses.
m=J('v54_MATH_CHECKPOINT.json');ck('math_gate',m['math_gate']=='PASS-GATES1-6-FROZEN')
ck('p46_global',m['pascal_noether']['status']=='PROVEN-GLOBAL-15-COEFFICIENT-IDENTITY' and m['pascal_noether']['identity']=='rho=(5/2)p_corrected')
ck('gl3_exact',m['pascal_noether']['GL3']=='PASS-EXACT')
ck('polar_projector',m['hypertensor']['polar_projector']=='PROVEN-EXACT')
ck('centered_jet',m['hypertensor']['centered_amplitude_jet']=='PROVEN-EXACT')
ck('no_time',m['hypertensor']['time_identification']=='SEPARATED/NOT-ASSERTED')
ck('fig2_nt',m['selling']['Fig2_absolute_binding'].startswith('NT/'))
ck('fig3_fail_closed',m['selling']['Fig3_direct_raster']=='2/91' and m['selling']['Fig3_not_propagated']=='89/91')
ck('r2_inj_refuted',m['r2pure13']['object_injective']=='REFUTED-TYPED')
ck('r2_strict_refuted',m['r2pure13']['strict_edge']=='REFUTED-TYPED')
ck('r2_1skel',m['r2pure13']['noninjective_1_skeleton'].startswith('PROVEN-EXISTS'))
ck('r2_2cell_nt',m['r2pure13']['source_2cells'].startswith('NT/'))
ck('r4_B6',m['retained_guards']['r4']=='NOT-OPENED' and m['retained_guards']['B6']=='NOT-OPENED')
ck('3D5D_pending',m['retained_guards']['Selling_3D_to_5D']=='CONTEXT/PENDING-MORPHISM')
# TeX append-only reconstruction.
v53=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v53.tex').read_text(encoding='utf-8')
v54=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v54.tex').read_text(encoding='utf-8')
ins=(B/'v54_SECTION_INSERT.tex').read_text(encoding='utf-8')
pascal='\\bibitem{Pascal1905}\nE. Pascal, \\emph{Contributo alla teoria della forma ternaria biquadratica}, Atti della Reale Accademia dei Lincei / memoria source utilisée dans le corpus, 1905.\n\n'
ck('tex_append_only',v54.replace(ins+'\n','',1).replace(pascal,'',1)==v53)
ck('tex_v54_before_cert',v54.index('\\section{v54 :')<v54.index('\\section{Certification HT+NT du présent formalisme}'))
# PDF QA.
q=J('v54_PDF_QA_REPORT.json');ck('pdf_qa',q['status']=='PASS' and q['pdf_pages']==277 and q['latex_overfull_count']==0)
ck('pdf_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v54.pdf')==q['pdf_sha256'])
ck('pdf_render_all',q['all_pages_rendered'] and q['automated_render_scan']['status']=='PASS' and q['automated_render_scan']['flags']==[])
o=J('v54_SECTION_ORDER_AUDIT.json');ck('section_order',o['status']=='PASS' and o['v54_before_certification'] and o['bibliography_physically_final'])
# Formal/source rigor.
fc=J('v54_FORMAL_CONTEXT_ATTESTATION.json');ck('formal_context',fc['status'].startswith('PASS-FORMAL-CONTEXT-LOCK') and fc['document_count']==9)
ra=J('v54_formal_rigor_audit.json');ck('formal_rigor',ra['status']=='PASS' and all(ra['checks'].values()))
sc=J('v54_SOURCE_CONTEXT.json');ck('source_context',sc['status']=='PASS')
# Tree overlay.
bt=J('v54_BUNDLE_TREE_AUDIT.json');ck('tree_status',bt['status']=='PASS')
ck('tree_dirs',all((B/d).is_dir() for d in bt['required_directories']))
ck('fragmentation_preseal',bt['total_uncompressed_bytes']<200*1024*1024)
# Summary.
industrial=len(checks); combined=8619+387+24+industrial
rep={'version':'v54','status':'PASS','predecessor_checks':8619,'v54_research_checks':387,'freeze_checks':24,'industrial_direct_checks':industrial,'combined_checks':combined,'failed':[],'direct_checks':checks,'pdf_sha256':sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v54.pdf'),'tex_sha256':sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v54.tex')}
if os.environ.get('V54_WRITE_REPORT')=='1':(B/'v54_VERIFICATION_REPORT.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps({k:rep[k] for k in ['status','predecessor_checks','v54_research_checks','freeze_checks','industrial_direct_checks','combined_checks']},indent=2))
