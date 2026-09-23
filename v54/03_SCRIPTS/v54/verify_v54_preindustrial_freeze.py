#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,subprocess,sys
B=Path(__file__).resolve().parent
checks=[]
def ck(n,c,d=None):
 checks.append({'name':n,'pass':bool(c),'detail':d})
 if not c: raise AssertionError(n)
def J(n):return json.loads((B/n).read_text(encoding='utf-8'))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
# Parent sealed anchors.
ck('v53_tex_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v53.tex')=='ba2a21715b00962d33e7ca7fe2a55c2928c7759b74d6aea659d04a7bb31b7f6a')
ck('v53_pdf_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v53.pdf')=='7db8718e541fab2d9e3ecbb0660c585f348b768ab7deb944d8d164b76bf75f0c')
tok=J('v53_ATOMIC_COMPLETION_TOKEN.json');ck('v53_closed',tok['publication_state']=='CLOSED/APPEND-ONLY')
# Research gate ledgers.
r1=J('v54_RESEARCH_VERIFICATION_REPORT.json'); ck('g1_41',r1['status']=='PASS' and r1['passed']==r1['total']==41)
r2=J('v54_GATE2_VERIFICATION_REPORT.json'); ck('g2_40',r2['status']=='PASS' and r2['check_count']==40)
r3=J('v54_GATE3_VERIFICATION_REPORT.json'); ck('g3_35',r3['status']=='PASS' and r3['new_passed']==r3['new_check_count']==35)
r4=J('v54_GATE4_VERIFICATION_REPORT.json'); ck('g4_56',r4['status']=='PASS' and r4['passed']==56)
r5=J('v54_GATE5_VERIFICATION_REPORT.json'); ck('g5_58',r5['status']=='PASS' and r5['passed']==58)
r6=J('v54_GATE6_CLOSURE_REPORT.json');ck('g6_77_closed',r6['status']=='PASS' and r6['closure_checks_passed']==77 and r6['gate_closed'])
# Replay gate6 closure directly.
p=subprocess.run([sys.executable,str(B/'verify_v54_gate6_closure.py')],cwd=B,capture_output=True,text=True)
ck('g6_closure_replay',p.returncode==0,p.stderr[-500:])
# Frozen statuses.
g=J('v54_GATE6_CLOSURE_CHECKPOINT.json')
ck('p46_closed',g['absorbed_results']['Pascal_Noether_lane_A']['status']=='CLOSED-ALGEBRAICALLY/PROVEN-GLOBAL')
ck('r2_injective_refuted',g['absorbed_results']['R2PURE13_to_SellingDeck17']['object_injective']=='REFUTED-TYPED')
ck('r2_1skel_proven',g['absorbed_results']['R2PURE13_to_SellingDeck17']['noninjective_1_skeleton'].startswith('PROVEN-EXISTS'))
ck('r2_2cell_nt',g['absorbed_results']['R2PURE13_to_SellingDeck17']['source_2cell_compatibility'].startswith('NT/'))
ck('fig2_nt',g['absorbed_results']['Fig2']['S_specific_raster_edge'].startswith('NT/'))
ck('fig3_2_91',g['absorbed_results']['Fig3']['direct_raster_bound']=='2/91' and g['absorbed_results']['Fig3']['not_propagated']=='89/91')
ck('r4_B6_closed',g['absorbed_results']['retained_guards']['r4']=='NOT-OPENED' and g['absorbed_results']['retained_guards']['B6']=='NOT-OPENED')
ck('old_pascal_blocker_superseded',g['supersession_ledger']['old_blocker_Pascal_dictionary_or_p46'].startswith('SUPERSEDED'))
ck('old_r2_blocker_superseded',g['supersession_ledger']['old_blocker_R2PURE_no_mixed_higher_cell'].startswith('SUPERSEDED'))
# Exact document append-only reconstruction.
v53=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v53.tex').read_text(encoding='utf-8')
v54=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v54.tex').read_text(encoding='utf-8')
ins=(B/'v54_SECTION_INSERT.tex').read_text(encoding='utf-8')
pascal='\\bibitem{Pascal1905}\nE. Pascal, \\emph{Contributo alla teoria della forma ternaria biquadratica}, Atti della Reale Accademia dei Lincei / memoria source utilisée dans le corpus, 1905.\n\n'
ck('v54_section_present','\\section{v54 :' in v54)
ck('v54_before_certification',v54.index('\\section{v54 :')<v54.index('\\section{Certification HT+NT du présent formalisme}'))
ck('v54_append_only_exact',v54.replace(ins+'\n','',1).replace(pascal,'',1)==v53)
# Parent 1007 bytes still exact.
pm=J('v54_PARENT_V53_1007_SHA256.json'); bad=[]
for rel,meta in pm['files'].items():
 p=B/rel
 if not p.exists() or sha(p)!=meta['sha256'] or p.stat().st_size!=meta['bytes']:bad.append(rel)
ck('parent_1007_byte_exact',not bad,bad[:10])
# Preindustrial freeze checkpoint.
f=J('v54_PREINDUSTRIAL_FREEZE_CHECKPOINT.json')
ck('freeze_status',f['status']=='PASS' and f['publication_state']=='PREINDUSTRIAL-FROZEN/NONPUBLISHABLE')
rep={'version':'v54','phase':'PREINDUSTRIAL FREEZE','status':'PASS','direct_checks':len(checks),'passed':sum(x['pass'] for x in checks),'failed':[x for x in checks if not x['pass']], 'parent_v53_combined_checks':8619,'v54_gate_research_checks':387,'next_action':'ATOMIC INDUSTRIAL CYCLE'}
(B/'v54_PREINDUSTRIAL_FREEZE_REPORT.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps(rep,indent=2))
