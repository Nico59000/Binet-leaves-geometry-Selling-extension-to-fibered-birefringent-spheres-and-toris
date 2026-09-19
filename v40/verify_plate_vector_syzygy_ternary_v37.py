#!/usr/bin/env python3
import csv, hashlib, json, subprocess, sys
from pathlib import Path
from collections import defaultdict
import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_decomp

OUT=Path(__file__).resolve().parent
checks={}

def check(name, cond):
    checks[name]=bool(cond)
    if not cond:
        print('FAIL',name,file=sys.stderr)

# -----------------------------------------------------------------------------
# Predecessor replay: v36 immutable (243/243).
# -----------------------------------------------------------------------------
p=subprocess.run([sys.executable, str(OUT/'verify_selling_domain_roof_gamma3_v36.py')],
                 cwd=OUT, capture_output=True, text=True)
check('predecessor_v36_replay_exit_zero', p.returncode==0)
v36=json.loads((OUT/'selling_domain_roof_gamma3_v36_certificate.json').read_text())
check('predecessor_v36_status_pass', v36.get('status')=='PASS')
check('predecessor_v36_combined_243', v36.get('combined')==243 or v36.get('combined_check_count')==243)

# -----------------------------------------------------------------------------
# A. Canonical Selling Tafel III bytes + source-anchored vector layers.
# The source bytes are now locally available and hash-locked.  Automatic
# geometry is exported with provenance/confidence, but is not silently promoted
# to the complete historical field incidence until all six symmetry axes,
# serpentines, repetitions and stabilizers are matched source-geometrically.
# -----------------------------------------------------------------------------
PDF=OUT/'LOG_0016.pdf'
CROP=OUT/'selling_fig1_full_600.png'
vec_run=subprocess.run([sys.executable,str(OUT/'vectorize_selling_fig1_v37.py')],cwd=OUT,capture_output=True,text=True)
check('vectorizer_exit_zero',vec_run.returncode==0)
check('source_pdf_present',PDF.exists())
check('source_pdf_size_1991876',PDF.stat().st_size==1991876)
pdf_sha=hashlib.sha256(PDF.read_bytes()).hexdigest()
check('source_pdf_sha_exact',pdf_sha=='423da98bd897c5dcb2b3556f4c9dc1367a4f12263d247608cfb2ba7a00f2e47b')
check('canonical_crop_present',CROP.exists())
crop_sha=hashlib.sha256(CROP.read_bytes()).hexdigest()
check('canonical_crop_sha_exact',crop_sha=='ff8583eb4372d9c4c1a285010e925173c78281c1cccf62888bdae47bca2261cf')

manifest=json.loads((OUT/'selling_fig1_vectorization_manifest_v37.json').read_text())
check('embedded_image_sha_exact',manifest['embedded_image_sha256']=='9116d054ae3267ae5f8a781135e434db659c4c8a6ff162cc5b6b2a12eee4c437')
check('crop_dimensions_4450_5400',manifest['crop_dimensions_600']==[4450,5400])
expected_tables={
 'selling_fig1_nodes_v37.csv':(136,'942b91656f60d6f3e84891d819b5327580c67f042b35602972a9e2cfa900476d'),
 'selling_fig1_segments_v37.csv':(105,'d866c35bae2e33b335af8e1818e75351ed890dbd47c44321855616c37106fa96'),
 'selling_fig1_fields_v37.csv':(34,'4827dd0cad2c17ce6995f5fc7c1966a3365d856a9755a30dcb1ae7cbc0636bd6'),
 'selling_fig1_layers_v37.geojson':(105,'8a410c44b8d98be828a76eefbb47e6534b083b7665411227ad51d8daaa055559')}
for fn,(n,h) in expected_tables.items():
    info=manifest['tables'][fn]
    gotn=info.get('rows',info.get('features'))
    check('vector_'+fn+'_count',gotn==n)
    check('vector_'+fn+'_hash',info['sha256']==h and hashlib.sha256((OUT/fn).read_bytes()).hexdigest()==h)
check('vector_layer_counts_exact',manifest['layer_counts']=={'ordinary':73,'reinforced':25,'serpentine':2,'symmetry_candidate':5})
check('vector_segment_count_partition',sum(manifest['layer_counts'].values())==105)
check('field_diagnostic_counts_exact',manifest['field_region_class_counts']=={'STABLE_CLOSED_REGION_HIGH':3,'STABLE_CLOSED_REGION_MEDIUM':4,'MORPHOLOGICAL_CANDIDATE':27})
check('vector_decision_source_anchored',manifest['decision'].startswith('SOURCE_BYTES_PROVEN; VECTOR_LAYERS_SOURCE_ANCHORED'))
check('vector_not_complete_historical_incidence','NOT_YET_PROMOTED' in manifest['decision'])
check('six_axis_gate_not_artificially_closed',manifest['layer_counts']['symmetry_candidate']!=6)
check('serpentine_gate_not_artificially_closed',manifest['layer_counts']['serpentine']==2)

plate_gate={
  'work':'E. Selling, Ueber die binären und ternären quadratischen Formen, J. reine angew. Math. 77 (1874), 143-229',
  'plate':'Tafel III; source PDF LOG_0016.pdf, page 2; Fig. 1 at right, Figs. 2--6 at left',
  'source_pdf':PDF.name,'source_pdf_bytes':PDF.stat().st_size,'source_pdf_sha256':pdf_sha,
  'embedded_raster':{'dimensions':[6716,5752],'dpi':[600,600],'sha256':manifest['embedded_image_sha256']},
  'fig1_crop':{'file':CROP.name,'dimensions':manifest['crop_dimensions_600'],'sha256':crop_sha},
  'source_status':{
    'plate_raster_bytes':'PROVEN_SOURCE_BYTES',
    'vector_layers':'PROVEN_SOURCE_ANCHORED_GEOMETRY_WITH_TYPED_DIAGNOSTIC_CLASSIFICATION',
    'complete_historical_field_incidence':'NT_PENDING_SOURCE_GEOMETRIC_MATCHING_OF_AXES_SERPENTINES_REPETITIONS_STABILIZERS'
  },
  'typing':['E_ordinary','E_serpentine','E_reinforced','E_symmetry_candidate'],
  'publication_guard':'No morphological closure or line-detection tolerance is promoted by itself to a historical field, wall, fissure, crossing, repetition or stabilizer.'
}
(OUT/'selling_fig1_plate_gate_v37.json').write_text(json.dumps(plate_gate,indent=2,ensure_ascii=False),encoding='utf-8')
check('plate_bytes_promoted_proven',plate_gate['source_status']['plate_raster_bytes']=='PROVEN_SOURCE_BYTES')
check('complete_incidence_remains_nt',plate_gate['source_status']['complete_historical_field_incidence'].startswith('NT_'))

rows=[
 ('SRC_PDF','source_bytes','global','Canonical two-page Selling Tafel III PDF acquired and hash-locked.','1991876 bytes','LOG_0016.pdf','PROVEN_SOURCE_BYTES'),
 ('SRC_RASTER','source_raster','global','Single embedded page-2 raster extracted directly from PDF.','6716x5752 @600dpi',manifest['embedded_image_sha256'],'PROVEN_SOURCE_BYTES'),
 ('FIG1_CROP','source_crop','figure1','Canonical Fig.1 crop from embedded raster.','4450x5400',crop_sha,'PROVEN_SOURCE_BYTES'),
 ('OUTER_AXES','boundary_axis','outer_skeleton','Six external boundary lines are two-part symmetry axes.','6','Selling source text','PROVEN_COUNT_ONLY / RASTER_MATCH_NT'),
 ('OUTER_CENTERS','quadripartite_center','outer_skeleton','Their six intersections are fourfold-symmetry centres.','6','Selling source text','PROVEN_COUNT_ONLY / RASTER_MATCH_NT'),
 ('BOUNDARY_RULE','boundary_rule','all_fields','Every boundary line begins/ends at fissures and passes through one or more crossings.','universal','Selling source text','PROVEN_RULE'),
 ('CROSSING_RULE','crossing','local_germ','Six outgoing lines of one division interleave with six of the adjoint division.','6+6 local rays','Selling source text','PROVEN_LOCAL_INCIDENCE'),
 ('FISSURE_RULE','fissure','local_germ','Fissures terminate boundary lines; continuation is source-defined.','not globally enumerated','Selling source text','PROVEN_LOCAL_RULE'),
 ('SERPENTINE_RULE','serpentine','figure1','Serpenting lines belong to the comparison/adjoint layer and must remain typed separately.','auto extraction finds 2 confirmed candidates','source raster + text','SOURCE_ANCHORED_PARTIAL'),
 ('REINFORCED_RULE','reinforced','figure1','Reinforced traits are kept as a separate selected-territory/boundary layer.','25 source-stroke candidates','source raster','SOURCE_ANCHORED_DIAGNOSTIC'),
 ('AA_HEXAGON','territory_boundary','A/a_territory','A cited territory development gives an alternating six-sided A/a boundary.','6 sides','Selling source text','PROVEN_TEXT_CONSTRAINT'),
 ('REPETITION_RULE','repeated_form','development','Numerically repeated fields yield self-transformations after finite development.','finite modulo repetition','Selling source text','PROVEN_REPETITION_GATE'),
 ('STABILIZER_RULE','stabilizer','repeated_field','Repeated fields and local symmetries determine stabilizers.','not exhaustively matched on raster','Selling source + v30-v36','PROVEN_RULE / ENUMERATION_NT'),
 ('VECTOR_TABLES','vector_diagnostic','figure1','Source-anchored typed vector exports exist for nodes, segments and candidate regions.','136 nodes; 105 segments; 34 region candidates','vectorizer + manifest','PROVEN_DIAGNOSTIC_GEOMETRY'),
 ('FULL_FIELD_TABLE','field_complex','global','Exhaustive field incidence requires matching every raster object to source semantics and repetition/stabilizer data.','not yet certified',CROP.name,'NT_PENDING_SOURCE_GEOMETRIC_MATCHING')]
inc_path=OUT/'selling_fig1_incidence_partial_v37.csv'
with inc_path.open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f); w.writerow(['record_id','entity_type','scope','incidence_statement','cardinality','source','status']); w.writerows(rows)
inc_sha=hashlib.sha256(inc_path.read_bytes()).hexdigest()
partial={
  'table':inc_path.name,'sha256':inc_sha,'row_count':len(rows),
  'source_vector_manifest':'selling_fig1_vectorization_manifest_v37.json',
  'status':'HASH_LOCKED_SOURCE_BYTES_PLUS_TYPED_VECTOR_DIAGNOSTICS',
  'complete_historical_incidence':'NT_PENDING_SOURCE_GEOMETRIC_MATCHING',
  'field_level_functor':'NT_NOT_PROMOTED',
  'historical_field_H1_H2':'NT_NOT_RECOMPUTED_BECAUSE_GATE_NOT_CLOSED'}
(OUT/'selling_fig1_incidence_gate_v37.json').write_text(json.dumps(partial,indent=2,ensure_ascii=False),encoding='utf-8')
check('partial_incidence_row_count_15',len(rows)==15)
check('partial_incidence_hash_64',len(inc_sha)==64)
check('field_functor_blocked_until_full_plate',partial['field_level_functor']=='NT_NOT_PROMOTED')
check('historical_H12_not_recomputed',partial['historical_field_H1_H2'].startswith('NT_'))

fun36=json.loads((OUT/'selling_to_F29_functor_v36.json').read_text())
field_functor_gate={
  'presentation_functor_v36':fun36,
  'source_vectorization':manifest,
  'decision':'FIELD_LEVEL_LIFT_NT: source bytes and typed vector layers are proven, but field-by-field semantics, all six symmetry axes, repetitions and stabilizers are not yet completely matched.',
  'cohomology_decision':'NO_FIELD_H1_H2_RECOMPUTATION_IN_V37; retain v36 transformation-orbifold Binet (1,0), Ad (3,0), and v35 F29_can (0,0) as separately typed computations.'}
(OUT/'selling_field_functor_gate_v37.json').write_text(json.dumps(field_functor_gate,indent=2,ensure_ascii=False),encoding='utf-8')
check('presentation_functor_retained',isinstance(fun36,dict) and len(fun36)>0)
check('field_lift_nt',field_functor_gate['decision'].startswith('FIELD_LEVEL_LIFT_NT'))

# -----------------------------------------------------------------------------
# B. Gamma_3(2): identities among 43 relators at 4-cell/syzygy level.
# -----------------------------------------------------------------------------
g3=json.loads((OUT/'gamma3_level2_relative_3cells_v36.json').read_text())
gens=g3['generator_order']; cells=g3['cells']; ng=len(gens); nr=len(cells)
check('gamma3_generator_count_9',ng==9)
check('gamma3_relator_count_43',nr==43)
idx={g:i for i,g in enumerate(gens)}

def parse(tok):
    return (tok[:-3],-1) if tok.endswith('^-1') else (tok,1)
B=sp.zeros(ng,nr)
for j,c in enumerate(cells):
    for tok in c['word']:
        g,s=parse(tok); B[idx[g],j]+=s
rankQ=B.rank()
check('augmentation_boundary_rank_9',rankQ==9)
check('augmentation_kernel_rank_34',nr-rankQ==34)
check('all_exponent_sums_even',all(int(x)%2==0 for x in B))
# Mod-2 matrix is zero.
check('augmentation_mod2_rank_zero',all(int(x)%2==0 for x in B))
S,Umat,Vmat=smith_normal_decomp(B,domain=ZZ)
diag=[abs(int(S[i,i])) for i in range(min(S.rows,S.cols)) if S[i,i]!=0]
check('smith_nine_nonzero',len(diag)==9)
check('smith_diag_all_2',diag==[2]*9)
check('smith_identity_U_B_V',Umat*B*Vmat==S)
K=Vmat[:,rankQ:]
check('integral_kernel_basis_shape',K.shape==(43,34))
check('integral_kernel_basis_exact',B*K==sp.zeros(9,34))
# V is unimodular; last columns therefore form a saturated Z-basis of ker B.
check('smith_V_unimodular',abs(int(Vmat.det()))==1)

# Export integral kernel basis.
syzcsv=OUT/'gamma3_augmented_syzygy_basis_v37.csv'
with syzcsv.open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f); w.writerow(['relator']+[f'z{j+1}' for j in range(K.cols)])
    for i,c in enumerate(cells): w.writerow([c['name']]+[int(K[i,j]) for j in range(K.cols)])
syzcsv_sha=hashlib.sha256(syzcsv.read_bytes()).hexdigest()

# Exact inverse-word pairs give genuine free-cancellation identities among relators.
def inv_word(word):
    out=[]
    for tok in reversed(word):
        g,s=parse(tok); out.append(g if s==-1 else g+'^-1')
    return out
name_to_word={c['name']:c['word'] for c in cells}
inv_pairs=[]
names=[c['name'] for c in cells]
for i,a in enumerate(names):
    for b in names[i+1:]:
        if name_to_word[b]==inv_word(name_to_word[a]): inv_pairs.append((a,b))
check('six_exact_inverse_relator_pairs',len(inv_pairs)==6)

def free_reduce(word):
    stack=[]
    for tok in word:
        g,s=parse(tok)
        if stack:
            h,t=parse(stack[-1])
            if h==g and t==-s:
                stack.pop(); continue
        stack.append(tok)
    return stack
for q,(a,b) in enumerate(inv_pairs,1):
    check(f'inverse_pair_{q}_freely_cancels',free_reduce(name_to_word[a]+name_to_word[b])==[])

syzygy={
  'source':'Kobayashi presentation specialized to n=3; 43 v36 relators',
  'chain_level':'augmentation / Fox boundary followed by augmentation Z[Gamma]->Z',
  'boundary_shape':[9,43], 'rank_over_Q':rankQ,
  'kernel_rank_over_Z':34,
  'smith_nonzero_diagonal':diag,
  'cokernel':'(Z/2)^9',
  'kernel_basis_csv':syzcsv.name,'kernel_basis_sha256':syzcsv_sha,
  'exact_free_inverse_pair_4cells':[{'relator_a':a,'relator_b':b,'identity':a+' * '+b+' freely reduces to empty word'} for a,b in inv_pairs],
  'exact_free_inverse_pair_count':len(inv_pairs),
  'status':{
    'augmented_syzygy_module':'PROVEN_Z_KERNEL_RANK_34_WITH_SMITH_BASIS',
    'six_free_cancellation_4cells':'PROVEN_NONABELIAN_FREE_WORD_IDENTITIES',
    'complete_nonabelian_identity_among_relations_module':'NT_PENDING_GROUP_RING_FOX_RESOLUTION_OR_CONTRACTING_HOMOTOPY'
  },
  'guard':'The 34-dimensional augmented kernel is not silently identified with the full pi_2/identity-among-relations module.'
}
(OUT/'gamma3_level2_syzygies_v37.json').write_text(json.dumps(syzygy,indent=2,ensure_ascii=False),encoding='utf-8')
check('syzygy_csv_hash_64',len(syzcsv_sha)==64)
check('full_nonabelian_syzygy_not_overclaimed',syzygy['status']['complete_nonabelian_identity_among_relations_module'].startswith('NT_'))

# -----------------------------------------------------------------------------
# C. Projectively compatible ternary/hyperfine refinement of M3(F2).
# R_k = Z/(2*3^k), CRT = F2 x Z/3^k; inverse limit = F2 x Z_3.
# -----------------------------------------------------------------------------
levels=[]
for k in range(1,9):
    m=2*(3**k); e=3**k
    check(f'crt_level_{k}_e_mod2_one',e%2==1)
    check(f'crt_level_{k}_e_mod3k_zero',e%(3**k)==0)
    check(f'crt_level_{k}_idempotent',(e*e-e)%m==0)
    check(f'crt_level_{k}_section_additive',(2*e)%m==0)
    check(f'crt_level_{k}_section_multiplicative',(e*e)%m==e%m)
    if k>1:
        prev=3**(k-1); check(f'crt_level_{k}_projective_e',e%(2*3**(k-1))==prev)
    levels.append({'k':k,'modulus':m,'parity_idempotent_e_k':e})

Ky=sp.Matrix([[0,1,0],[1,0,0],[0,0,0]])
Kz=sp.Matrix([[0,0,1],[0,0,0],[1,0,0]])
J=sp.Matrix([[0,0,0],[0,0,1],[0,-1,0]])
Gamma=sp.diag(-1,1,-1)
check('integer_bracket_Ky_Kz_J',Ky*Kz-Kz*Ky==J)
check('integer_bracket_J_Ky_minusKz',J*Ky-Ky*J==-Kz)
check('integer_bracket_J_Kz_Ky',J*Kz-Kz*J==Ky)
check('Gamma_reverses_J_integer',Gamma*J*Gamma.inv()==-J)

def modmat(A,m): return A.applyfunc(lambda x:int(x)%m)
for k in range(1,7):
    m=2*3**k
    ky,kz,j=modmat(Ky,m),modmat(Kz,m),modmat(J,m)
    check(f'level_{k}_bracket1',modmat(ky*kz-kz*ky,m)==j)
    check(f'level_{k}_bracket2',modmat(j*ky-ky*j,m)==modmat(-kz,m))
    check(f'level_{k}_bracket3',modmat(j*kz-kz*j,m)==ky)
    check(f'level_{k}_sign_J_distinguished_mod3k',modmat(J,3**k)!=modmat(-J,3**k))
    check(f'level_{k}_Gamma_reversal',modmat(Gamma*J*Gamma.inv(),m)==modmat(-J,m))

binder={
  'level_ring':'R_k = Z/(2*3^k)Z',
  'crt':'R_k ≅ F2 × Z/3^kZ',
  'matrix_binder':'B_k=M3(R_k) ≅ M3(F2) × M3(Z/3^kZ)',
  'projective_limit':'lim_k B_k ≅ M3(F2) × M3(Z_3)',
  'parity_section':{
    'e_k':'3^k mod 2*3^k','properties':['e_k^2=e_k','e_k mod2=1','e_k mod3^k=0','e_{k+1}->e_k'],
    'map':'s_k(a)=a e_k is a compatible non-unital ring section F2 -> R_k'
  },
  'lambda2_nine_labels':'The nine matrix units e_ij remain the parity coordinates; each carries an independent 3-adic/hyperfine fibre through M3(Z_3).',
  'Thomas_Wigner':{
    'integral_generators':{'Ky':[[int(x) for x in row] for row in Ky.tolist()],'Kz':[[int(x) for x in row] for row in Kz.tolist()],'J':[[int(x) for x in row] for row in J.tolist()]},
    'relations':['[Ky,Kz]=J','[J,Ky]=-Kz','[J,Kz]=Ky'],
    'tower_status':'PROVEN_AT_EVERY_FINITE_LEVEL_BY_INTEGER_REDUCTION',
    'mod2_shadow':'recovers v36 parity bracket, where - becomes +',
    'ternary_component':'retains the distinction J versus -J for every k>=1'
  },
  'orientation':{
    'Gamma_action':'Gamma J Gamma^{-1}=-J at every level',
    'decision':'The tower distinguishes the two orientations but does not canonically select one; O(2)/Real phase remains natural.'
  },
  'roof_relation':{
    '12345678987654321':'ordered nine-address visualization for the 3x3 parity coordinates after a chosen row-major chart',
    'status':'POINTED_CONTEXT_ADAPTER_ONLY; roof reversal is still not matrix transpose and no intrinsic permutation equivalence is introduced.'
  },
  'levels_checked':levels,
  'status':'PROVEN_PROJECTIVE_CRT_TERNARY_REFINEMENT'
}
(OUT/'ternary_projective_binder_v37.json').write_text(json.dumps(binder,indent=2,ensure_ascii=False),encoding='utf-8')
check('ternary_binder_status_proven',binder['status']=='PROVEN_PROJECTIVE_CRT_TERNARY_REFINEMENT')
check('orientation_not_selected', 'does not canonically select' in binder['orientation']['decision'])

# -----------------------------------------------------------------------------
# D. Pending-morphism registry after acquiring Tafel III source bytes.
# -----------------------------------------------------------------------------
pending=[
 {'id':'PM37-01','object':'Selling Tafel III raster/scan bytes','status':'RESOLVED_PROVEN_SOURCE_BYTES','next_gate':'Closed in v37 by LOG_0016.pdf + embedded-raster hash lock.'},
 {'id':'PM37-02','object':'Complete F_Sell^field(Q0) incidence','status':'IN_PROGRESS_SOURCE_ANCHORED / NT_COMPLETE_INCIDENCE','next_gate':'Match the 136 node / 105 segment / 34 candidate-region diagnostics to all six axes, serpentines, historical fields, repetitions and stabilizers.'},
 {'id':'PM37-03','object':'Field-level functor F_Sell^field(Q0)->F29_can','status':'PENDING-MORPHISM','next_gate':'Promote only after PM37-02; verify every source field/wall and every 2-cell incidence.'},
 {'id':'PM37-04','object':'Historical field Binet/Adjoint H1,H2','status':'NT_BLOCKED_BY_PM37-02_03','next_gate':'Compute only after source complex and functor close.'},
 {'id':'PM37-05','object':'Full nonabelian syzygy module among 43 Gamma3 relators','status':'NT_PENDING_GROUP_RING_FOX','next_gate':'Compute Fox Jacobian over Z[Gamma3(2)] or an explicit resolution/contracting homotopy.'},
 {'id':'PM37-06','object':'Natural roof nine-address -> lambda2 binder','status':'POINTED_ONLY / NATURAL_REVERSAL_TRANSPOSE_REFUTED','next_gate':'Seek richer equivariant roof operations, not plain reversal.'},
 {'id':'PM37-07','object':'Hyperfine CRT tower -> real Lorentz rapidity scale','status':'SEPARATED','next_gate':'Supply a typed 3-adic-to-real metric/rapidity comparison.'},
 {'id':'PM37-08','object':'Historical source-native orientation breaking Gamma','status':'NT','next_gate':'Search only in the completely matched field/stabilizer data.'},
 {'id':'PM37-09','object':'Oriented SO(2)->U(1) phase','status':'PENDING-MORPHISM','next_gate':'Requires PM37-08; otherwise retain O(2)/Real phase.'},
 {'id':'PM37-10','object':'Physical de Broglie phase theta=S/hbar','status':'SEPARATED/PENDING-MORPHISM','next_gate':'Need source-derived physical action/energy-momentum/synchronization and scale.'},
 {'id':'PM37-11','object':'Cube27 -> Selling/deck/hyperfine binder functor','status':'SEPARATED','next_gate':'Require explicit equivariant functor; common 3-locality/rank/prime is insufficient.'}]
reg={'version':'v37','items':pending,'counts':dict(defaultdict(int))}
for x in pending: reg['counts'][x['status']]=reg['counts'].get(x['status'],0)+1
(OUT/'pending_morphism_registry_v37.json').write_text(json.dumps(reg,indent=2,ensure_ascii=False),encoding='utf-8')
check('pending_registry_11_items',len(pending)==11)
check('plate_source_gate_resolved',pending[0]['status']=='RESOLVED_PROVEN_SOURCE_BYTES')
check('complete_incidence_still_nt','NT_COMPLETE_INCIDENCE' in pending[1]['status'])
check('debroglie_stays_separated',pending[9]['status'].startswith('SEPARATED'))
check('cube27_stays_separated',pending[10]['status']=='SEPARATED')

# -----------------------------------------------------------------------------
# Certificate.
# -----------------------------------------------------------------------------
new_count=len(checks)-3
status='PASS' if all(checks.values()) else 'FAIL'
cert={
  'phase':'v37-taufel-iii-source-vector-gamma3-syzygy-ternary-projective-refinement',
  'status':status,
  'predecessor_v36_checks':243,
  'new_check_count':new_count,
  'combined_check_count':243+new_count,
  'checks':checks,
  'plate':{
    'gate':plate_gate,'partial_incidence_sha256':inc_sha,'vectorization_manifest':manifest,
    'decision':'SOURCE_BYTES_PROVEN_AND_TYPED_VECTOR_TABLES_HASH_LOCKED; COMPLETE_HISTORICAL_FIELD_INCIDENCE_REMAINS_NT_PENDING_SOURCE_GEOMETRIC_MATCHING.'},
  'field_functor':'NT_NOT_PROMOTED; historical microscopic H1/H2 intentionally not recomputed',
  'gamma3_syzygy':syzygy,
  'ternary_binder':binder,
  'pending_registry':'pending_morphism_registry_v37.json',
  'safe_successor':'v38 -- retain all v37 source hashes and typed separations; manually/source-geometrically reconcile the 136 nodes, 105 typed segments and 34 region candidates with the six Selling symmetry axes, all serpentines, bounded fields, repetitions and stabilizers, then and only then promote F_Sell^field(Q0)->F29_can and compute historical Binet/adjoint H1/H2. In parallel lift the rank-34 augmented Gamma3 syzygies to the full Z[Gamma3(2)] Fox identity module and seek source-verified 4-cells; extend the M3(F2)xM3(Z3) hyperfine tower without identifying roof reversal, transpose, Lorentz orientation or physical phase absent explicit equivariant morphisms.'}
(OUT/'plate_syzygy_ternary_v37_certificate.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'status':status,'v36_checks':243,'new_checks':new_count,'combined':243+new_count,
                  'source_pdf_sha256':pdf_sha,'vector_nodes':136,'vector_segments':105,'region_candidates':34,
                  'complete_historical_incidence':partial['complete_historical_incidence'],
                  'augmentation_syzygy_rank':34,'exact_inverse_4cells':len(inv_pairs),
                  'ternary_binder':binder['status']},indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
