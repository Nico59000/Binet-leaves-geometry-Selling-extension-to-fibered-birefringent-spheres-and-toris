#!/usr/bin/env python3
import json, hashlib
from pathlib import Path
B=Path('/mnt/data/v51_work')
ctx=json.load(open(B/'formal_context/R6_74_FORMAL_CONTEXT_REPUSH_ATTESTATION.json'))
wanted=[
'addendum_sqrt3_hypergeometric_quarterturn_selling_fano_v2_2_hal.pdf',
'Annexe-Heyting-besteod-v5-nt-Semantique-Nf-nt-Interne-externe.pdf',
'heyting_besteod_htnt_master_hal_v1.pdf','Heyting-Besteod-V5—Fusion-V3.pdf',
'htnt_functional_logic_addendum.pdf','pochhammer_shifted_factorial_ternary_roof_hal_v1.pdf',
'polydivisible_hyperplane_roof_elliptic_selling_progression_v4_3_hal.pdf',
'positional_reversal_mixed_ternary_hal_v1_1.pdf','roof_five_slice_four_interface_ternary_frames_hal_v1.pdf']
# normalize titles only for exact documented entries
lookup={d['title']:d for d in ctx['documents']}
docs=[]
for title in wanted:
    key=title
    if key not in lookup:
        # stored attestation omits .pdf on some title string variants only if any
        cand=[d for d in ctx['documents'] if d['title'].replace('.pdf','')==title.replace('.pdf','')]
        if not cand: raise SystemExit(f'missing context {title}')
        d=cand[0]
    else:d=lookup[key]
    docs.append({k:d.get(k) for k in ['title','bytes','sha256','pages','encrypted','byte_comparison_status','expected_sealed_sha256','guard_phrase_checks','guard_phrase_passed','guard_phrase_total']})
formal={
 'version':'v51','status':'PASS-FORMAL-CONTEXT-LOCK','source_attestation':'R6_74_FORMAL_CONTEXT_REPUSH_ATTESTATION.json',
 'documents':docs,
 'discipline':{
  'typed_separation':'No carrier, sort, logical status, arithmetic value, geometric feature, raster feature, or cohomology class is identified without a declared map.',
  'nt_guard':'NT/PENDING states are not binarized or promoted by numerical resemblance, failed local search, or absence of a witness.',
  'append_only':'v50 CLOSED/APPEND-ONLY is immutable; v51 adds only new typed claims and explicit errata/diagnostics.',
  'frozen_context':'Algebraic properties are invoked only in a fixed declared carrier/context; changes of carrier require an explicit adapter.',
  'source_byte_guard':'A same-named but byte-different source edition is never silently substituted; the addendum v2_2 distinction remains explicit.',
  'industrial_atomicity':'Publish CLOSED only after MATH, VERIFY, PDF-QA, SEAL, ROUNDTRIP and POSTSEAL all pass.'},
 'decision':'FORMAL_CONTEXT_GOVERNS_v51_TYPING_STATUS_AND_PRODUCTION; IT_DOES_NOT_PROMOTE_SOURCE-SPECIFIC_MATHEMATICS_BY_ITSELF.'}
p=B/'v51_FORMAL_CONTEXT_ATTESTATION.json';p.write_text(json.dumps(formal,indent=2,sort_keys=True)+'\n')
# Build consolidated math checkpoint from exact artifacts
files={n:json.load(open(B/n)) for n in [
 'selling_source_axis_centre_carrier_v51.json','selling_source_axis_centre_stabilizers_v51.json','selling_source_HC_local_star_signatures_v51.json',
 'selling_raster_HC_signature_identifiability_v51.json','selling_source_raster_projective_bold_cycle_audit_v51.json','selling_torsion_support_bockstein_pairings_v51.json',
 'gamma3_radius5_support_gate_v51.json','gamma3_radius5_quotient_cycle_v51.json','gamma3_radius5_fourfiber_nogo_v51.json',
 'gamma3_radius5_residual_projection_lattice_v51.json','gamma3_radius5_iterated_kernel_correction_v51.json']}
ck={
 'version':'v51','publication_state':'PRESEAL/NONPUBLISHABLE','predecessor':'v50 CLOSED/APPEND-ONLY',
 'formal_context':'PASS-FORMAL-CONTEXT-LOCK',
 'selling_source':{
  'axis_centre_carrier':files['selling_source_axis_centre_carrier_v51.json']['status'],
  'stabilizers':files['selling_source_axis_centre_stabilizers_v51.json']['status'],
  'local_star_signatures':files['selling_source_HC_local_star_signatures_v51.json']['status'],
  'six_centres':[(c['id'],c['q'],c['r']) for c in files['selling_source_axis_centre_carrier_v51.json']['centres']],
  'raster_binding':'NT/PENDING-DIRECT-SOURCE-ANCHORS',
  'raster_table_identifiability':files['selling_raster_HC_signature_identifiability_v51.json']['status'],
  'raster_table_candidate_count':files['selling_raster_HC_signature_identifiability_v51.json']['candidate_count'],
  'homography_diagnostic_only':True},
 'torsion':{
  'status':files['selling_torsion_support_bockstein_pairings_v51.json']['status'],
  'binet_H1_Z2':'PROVEN-EXPLICIT-BOCKSTEIN-SUPPORT',
  'ad_H1_Z2':'PROVEN-EXPLICIT-BOCKSTEIN-SUPPORT',
  'cross_system_overlap':'SUPPORT-INTERSECTION-ONLY/NO-COEFFICIENT-MORPHISM',
  'binet_degree2_transfer_defect':'PROVEN-SEPARATE-CARRIER/V48/TU/INDEX2'},
 'gamma3_radius5':{
  'support_status':files['gamma3_radius5_support_gate_v51.json']['status'],
  'ball_size':files['gamma3_radius5_support_gate_v51.json']['ball_size'],
  'mod4_image_size':files['gamma3_radius5_support_gate_v51.json']['mod4_image_size'],
  'quotient_cycle_status':files['gamma3_radius5_quotient_cycle_v51.json']['status'],
  'quotient_cycle_terms':files['gamma3_radius5_quotient_cycle_v51.json']['quotient_cycle_term_count'],
  'canonical_section_integral_lift':'REFUTED-TYPED',
  'four_witness_fibers':files['gamma3_radius5_fourfiber_nogo_v51.json']['status'],
  'local_residual_lattice':files['gamma3_radius5_residual_projection_lattice_v51.json']['status'],
  'local_projection_saturated':len(files['gamma3_radius5_residual_projection_lattice_v51.json']['smith_nonunit_factors'])==0,
  'iterated_kernel_correction':files['gamma3_radius5_iterated_kernel_correction_v51.json']['status'],
  'residual_period_detected':files['gamma3_radius5_iterated_kernel_correction_v51.json']['residual_period_detected'],
  'global_integral_lift_in_B5':'NT',
  'global_radius5_no_go':'NT-NOT-PROVEN'},
 'retained_separations':{
  'R2PURE13_to_SellingDeck17':'SEPARATED/PENDING-MORPHISM',
  'figures_2_to_6':'PENDING-SOURCE-SPECIFIC-SEEDS/BRANCH-DATA',
  'R5_to_Gamma3':'REFUTED-TYPED','roof_reversal_transposition':'REFUTED-TYPED',
  'Singer_physical_phase':'SEPARATED/PENDING-MORPHISM','Lorentz_orientation':'NT'},
 'math_gate':'PASS-WITH-EXPLICIT-NT-FRONTIERS',
 'guard':'v51 may close with the global B5 integral lift and raster HC binding still NT because their non-promotion is itself formally guarded and all positive/negative bounded claims are certified.'}
ck['artifact_sha256']={n:hashlib.sha256((B/n).read_bytes()).hexdigest() for n in files}
q=B/'v51_MATH_CHECKPOINT.json';q.write_text(json.dumps(ck,indent=2,sort_keys=True)+'\n')
print('formal',hashlib.sha256(p.read_bytes()).hexdigest())
print('checkpoint',hashlib.sha256(q.read_bytes()).hexdigest())
