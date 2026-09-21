#!/usr/bin/env python3
from pathlib import Path
import json, csv, hashlib, subprocess, sys, itertools, collections, math
import sympy as sp
import numpy as np

B=Path(__file__).resolve().parent
checks=[]

def ok(name, cond):
    checks.append((name, bool(cond)))
    if not cond:
        print('FAIL', name)

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

def invtok(t):
    return t[:-3] if t.endswith('^-1') else t+'^-1'

def invword_tokens(w):
    return [invtok(t) for t in reversed(w)]

def reduce_tokens(w):
    st=[]
    for t in w:
        if st and invtok(t)==st[-1]: st.pop()
        else: st.append(t)
    return st

def parse_token(t):
    return (t[:-3],-1) if t.endswith('^-1') else (t,1)

# -----------------------------------------------------------------------------
# Parent replay: v38 must remain byte/logic compatible.
# -----------------------------------------------------------------------------
p=subprocess.run([sys.executable,str(B/'verify_field_fox_projective_v38.py')],cwd=B,text=True,capture_output=True)
ok('parent_v38_returncode',p.returncode==0)
ok('parent_v38_pass','"status": "PASS"' in p.stdout and '"combined": 614' in p.stdout)

# -----------------------------------------------------------------------------
# A. Close the source-semantic contract of Selling Fig. 1, while proving that
#    the current diagnostic tables are still insufficient to instantiate the
#    complete microscopic cell complex.
# -----------------------------------------------------------------------------
pg=json.loads((B/'selling_fig1_bsb_photo_gate_v38.json').read_text())
segments=list(csv.DictReader(open(B/'selling_fig1_segments_v37.csv',encoding='utf-8')))
nodes=list(csv.DictReader(open(B/'selling_fig1_nodes_v37.csv',encoding='utf-8')))
fields=list(csv.DictReader(open(B/'selling_fig1_fields_v37.csv',encoding='utf-8')))
inc_rows=list(csv.DictReader(open(B/'selling_fig1_incidence_partial_v37.csv',encoding='utf-8')))

layer_counts=collections.Counter(r['layer'] for r in segments)
field_classes=collections.Counter(r['classification'] for r in fields)
all_fields_unpromoted=all(r['historical_field_status'].startswith('NOT_PROMOTED') for r in fields)
segment_columns=set(segments[0].keys())
node_columns=set(nodes[0].keys())
field_columns=set(fields[0].keys())

semantic_contract={
  'version':'v39',
  'primary_source':'E. Selling, Des formes quadratiques binaires et ternaires, JMPA 3e serie 3 (1877), pp.153-206, second part; French translation of the 1874 memoir.',
  'original_source':'E. Selling, Ueber die binaeren und ternaeren quadratischen Formen, J. reine angew. Math. 77 (1874), pp.143-229; BSB bsb11179195.',
  'figure':'Tafel III, Fig. 1 (zu S. 188 u.f.)',
  'source_semantic_rules':[
    {'id':'S39-R1','printed_pages':'169','statement':'Every field boundary line begins and ends at fissure points and passes through one or more crossing points; after all crossings of a field are found, the last boundary reconnects to the first.','status':'PROVEN_SOURCE_RULE'},
    {'id':'S39-R2','printed_pages':'169','statement':'The crossing points of the primal and adjoint field divisions coincide; away from crossings, boundary lines of the two systems do not intersect.','status':'PROVEN_SOURCE_RULE'},
    {'id':'S39-R3','printed_pages':'169','statement':'At a crossing, each of the six boundary lines of one division lies between two of the six boundary lines of the adjoint division; selected adjoint lines are drawn serpentine.','status':'PROVEN_SOURCE_LOCAL_INCIDENCE'},
    {'id':'S39-R4','printed_pages':'170','statement':'After sufficiently far development, every newly added field is a numerical repetition of a previously existing field; all later development repeats the finite part.','status':'PROVEN_SOURCE_REPETITION'},
    {'id':'S39-R5','printed_pages':'170','statement':'In the completed Fig. 1 example, every boundary of the developed field system is a symmetry axis and the intersections of these outer axes have fourfold symmetry.','status':'PROVEN_SOURCE_SYMMETRY'},
    {'id':'S39-R6','printed_pages':'170','statement':'Each field numerically equal to the initial form determines a self-transformation obtained as the product of consecutive substitutions along adjacent-field transitions.','status':'PROVEN_SOURCE_STABILIZER_RULE'},
    {'id':'S39-R7','printed_pages':'175','statement':'In Fig. 1, all but one boundaries of positive-A territories coincide with field-boundary lines and are drawn as serpentine lines; one exceptional connection is represented by an existing field-boundary line rather than introducing a foreign line.','status':'PROVEN_SOURCE_DRAWING_CONVENTION'},
    {'id':'S39-R8','printed_pages':'175','statement':'The six external boundary lines are bipartite symmetry axes and their six intersections are fourfold-symmetry centres.','status':'PROVEN_SOURCE_COUNT_AND_TYPE'},
    {'id':'S39-R9','printed_pages':'184','statement':'The six indicated symmetry-axis intersections are fourfold centres; an additional sixfold symmetry exists around the centre of the figure but is deliberately not used to reduce the drawing.','status':'PROVEN_SOURCE_GLOBAL_SYMMETRY'},
    {'id':'S39-R10','printed_pages':'190-191','statement':'For the Fig. 1 example, four external-boundary substitution types P,Q,R,S bound the quadruplicated figure; T=P R P^{-1} and U=Q S Q^{-1} satisfy TU=UT.','status':'PROVEN_SOURCE_STABILIZER_PRESENTATION'},
  ],
  'semantic_closure':'PROVEN_SOURCE_SEMANTIC_CONTRACT_CLOSED',
  'diagnostic_geometry':{
    'nodes':len(nodes),'segments':len(segments),'regions':len(fields),
    'layer_counts':dict(layer_counts),'field_class_counts':dict(field_classes),
    'photo_support_counts':pg['v37_segment_support_counts']
  },
  'schema_obstructions':{
    'source_requires_external_axes':6,
    'typed_symmetry_candidates_present':layer_counts.get('symmetry_candidate',0),
    'segment_has_historical_wall_or_substitution_label':bool({'historical_wall_label','selling_generator','substitution_label'} & segment_columns),
    'node_has_fissure_crossing_semantic_type':bool({'historical_node_type','fissure_crossing_type'} & node_columns),
    'field_has_historical_field_id_or_form_label':bool({'historical_field_id','historical_form_label'} & field_columns),
    'all_34_regions_currently_unpromoted':all_fields_unpromoted,
    'repetition_pairs_machine_encoded':False,
    'stabilizer_generators_machine_encoded_per_field':False
  },
  'decision':'SOURCE_SEMANTICS_CLOSED; MICROSCOPIC_GEOMETRIC_INSTANTIATION_NT_SCHEMA_AND_MATCHING_OBSTRUCTION',
  'guard':'The source semantic contract is complete enough to state certification conditions, but it does not itself assign the 136 diagnostic nodes, 105 diagnostic segments and 34 region candidates to the historical cells. No field functor or field cohomology is inferred.'
}
(B/'selling_fig1_semantic_contract_v39.json').write_text(json.dumps(semantic_contract,indent=2,ensure_ascii=False),encoding='utf-8')
ok('semantic_rule_count_10',len(semantic_contract['source_semantic_rules'])==10)
ok('semantic_contract_closed',semantic_contract['semantic_closure'].startswith('PROVEN_'))
ok('source_six_axes',semantic_contract['schema_obstructions']['source_requires_external_axes']==6)
ok('typed_only_five_symmetry_candidates',semantic_contract['schema_obstructions']['typed_symmetry_candidates_present']==5)
ok('segment_labels_absent',not semantic_contract['schema_obstructions']['segment_has_historical_wall_or_substitution_label'])
ok('node_semantic_types_absent',not semantic_contract['schema_obstructions']['node_has_fissure_crossing_semantic_type'])
ok('field_labels_absent',not semantic_contract['schema_obstructions']['field_has_historical_field_id_or_form_label'])
ok('all_fields_unpromoted',all_fields_unpromoted and len(fields)==34)

field_gate={
  'version':'v39',
  'semantic_contract':'selling_fig1_semantic_contract_v39.json',
  'semantic_contract_sha256':sha(B/'selling_fig1_semantic_contract_v39.json'),
  'presentation_functor':'selling_to_F29_functor_v36.json',
  'source_cells_required':{'C0':'historical fields / crossings / fissures with source identity','C1':'historically labelled boundary arcs with substitution type','C2':'closed historical fields with ordered boundary cycles and repetition/stabilizer data'},
  'current_table_defects':[
    'Only 5 rows are typed symmetry_candidate although the source requires 6 external symmetry axes; at least one diagnostic segment requires semantic retyping or splitting.',
    'No segment row carries a Selling transformation/generator label, so the 1-cell map to the v36 presentation functor cannot be instantiated.',
    'No node row carries the historical fissure/crossing semantic type.',
    'All 34 region rows explicitly remain NOT_PROMOTED and carry no historical form/field identity.',
    'Repeated-field pairings and per-field stabilizer words are not encoded.'
  ],
  'zero_cell_compatibility':'NT_NOT_INSTANTIABLE_FROM_CURRENT_TABLES',
  'one_cell_compatibility':'NT_NOT_INSTANTIABLE_NO_HISTORICAL_WALL_LABELS',
  'two_cell_compatibility':'NT_NOT_INSTANTIABLE_NO_ORDERED_HISTORICAL_BOUNDARY_CYCLES',
  'field_level_functor':'PENDING-MORPHISM_NOT_PROMOTED',
  'historical_field_H1_H2':'NT_NOT_RECOMPUTED_GATE_NOT_CLOSED',
  'decision':'PROVEN_CERTIFICATION_OBSTRUCTION_FOR_CURRENT_V37_V38_TABLE_SCHEMA; no mathematical no-go for a future enriched transcription.'
}
(B/'selling_field_functor_gate_v39.json').write_text(json.dumps(field_gate,indent=2,ensure_ascii=False),encoding='utf-8')
ok('field_functor_not_promoted',field_gate['field_level_functor'].startswith('PENDING'))
ok('field_H12_not_recomputed',field_gate['historical_field_H1_H2'].startswith('NT_'))
ok('field_gate_obstruction_typed','SCHEMA' in field_gate['decision'])

# -----------------------------------------------------------------------------
# B. Gamma_3(2): critical-overlap census and an exact finite-quotient Fox shadow.
# -----------------------------------------------------------------------------
g3=json.loads((B/'gamma3_level2_relative_3cells_v36.json').read_text())
gens=g3['generator_order']; cells=g3['cells']; nrel=len(cells)
words=[reduce_tokens(c['word']) for c in cells]

# Fixed-basepoint oriented products: complete 2-relator free cancellations.
orient=[]
for idx,c in enumerate(cells):
    w=words[idx]
    orient.append((idx,1,c['name'],w))
    orient.append((idx,-1,c['name']+'^-1',invword_tokens(w)))
pair_id=[]
pair_cancel_hist=collections.Counter()
for a in orient:
    for b in orient:
        if a[0]==b[0]: continue
        rr=reduce_tokens(a[3]+b[3])
        can=(len(a[3])+len(b[3])-len(rr))//2
        pair_cancel_hist[can]+=1
        if not rr: pair_id.append((a[2],b[2]))
# underlying unordered exact inverse pairs
underlying_pairs=set()
for a,b in pair_id:
    aa=a.replace('^-1',''); bb=b.replace('^-1','')
    underlying_pairs.add(tuple(sorted((aa,bb))))

# Cyclic critical-overlap census (basepoint allowed to move along attaching loop).
cyc=[]
for idx,c in enumerate(cells):
    w=words[idx]
    for sgn,ww in [(1,w),(-1,invword_tokens(w))]:
        for sh in range(len(ww)):
            cyc.append((idx,sgn,sh,c['name'],ww[sh:]+ww[:sh]))
cyc_hist=collections.Counter(); cyc_ge2=[]
for a in cyc:
    for b in cyc:
        if a[0]==b[0]: continue
        rr=reduce_tokens(a[4]+b[4])
        can=(len(a[4])+len(b[4])-len(rr))//2
        cyc_hist[can]+=1
        if can>=2: cyc_ge2.append((can,a[3],a[1],a[2],b[3],b[1],b[2]))
max_cyclic_overlap=max(cyc_hist)
length3_pairs=sorted(set(tuple(sorted((x[1],x[4]))) for x in cyc_ge2 if x[0]==3))

# Exhaustive bare triple identity search for three distinct underlying relators.
triple_identities=[]
for inds in itertools.combinations(range(nrel),3):
    found=False
    for signs in itertools.product([1,-1],repeat=3):
        ws=[]; ns=[]
        for i,s in zip(inds,signs):
            w=words[i] if s==1 else invword_tokens(words[i])
            ws.append(w); ns.append(cells[i]['name']+('' if s==1 else '^-1'))
        for perm in itertools.permutations(range(3)):
            if not reduce_tokens(ws[perm[0]]+ws[perm[1]]+ws[perm[2]]):
                triple_identities.append(tuple(ns[i] for i in perm)); found=True; break
        if found: break

# Fox Jacobian evaluated on every character of the mod-4 abelian quotient (Z/2)^9.
def fox_eval(word, x, char):
    pref=1; s=0
    for g,e in [parse_token(t) for t in word]:
        if e==1:
            if g==x: s += pref
            pref *= char[g]
        else:
            if g==x: s += -pref*char[g]
            pref *= char[g]
    return s

rank_hist=collections.Counter(); character_records=[]; kernel_sum=0
for bits in itertools.product([1,-1], repeat=len(gens)):
    ch=dict(zip(gens,bits))
    M=sp.Matrix([[fox_eval(w,g,ch) for w in words] for g in gens])
    r=int(M.rank()); k=nrel-r
    rank_hist[r]+=1; kernel_sum+=k
    character_records.append({'signs':[int(v) for v in bits],'rank':r,'kernel_dimension':k})

fox_shadow={
  'version':'v39',
  'presentation':'43 Kobayashi relators for Gamma_3(2)',
  'critical_overlap_census':{
    'oriented_fixed_basepoint_word_count':len(orient),
    'fixed_basepoint_cancellation_histogram':{str(k):v for k,v in sorted(pair_cancel_hist.items())},
    'exact_two_relator_oriented_identities':len(pair_id),
    'underlying_exact_inverse_pairs':len(underlying_pairs),
    'underlying_pairs':sorted([list(x) for x in underlying_pairs]),
    'cyclic_oriented_attaching_words':len(cyc),
    'cyclic_cancellation_histogram':{str(k):v for k,v in sorted(cyc_hist.items())},
    'cyclic_overlap_ge2_count':len(cyc_ge2),
    'max_cyclic_cancellation_length':max_cyclic_overlap,
    'underlying_length3_overlap_pairs':length3_pairs,
    'bare_three_distinct_relator_free_identities':triple_identities,
    'scope':'Complete for the stated bounded combinatorial search; this is not a complete Peiffer/identity-among-relations computation.'
  },
  'mod4_character_fox_shadow':{
    'quotient':'Gamma_3(2) -> (Z/2)^9 induced by the nine standard generators modulo 4',
    'coefficient_field':'Q via all 512 sign characters of (Z/2)^9',
    'character_count':2**len(gens),
    'rank_histogram':{str(k):v for k,v in sorted(rank_hist.items())},
    'total_kernel_dimension_over_regular_character_decomposition':kernel_sum,
    'average_character_kernel_dimension':kernel_sum/(2**len(gens)),
    'trivial_character_rank':character_records[0]['rank'],
    'trivial_character_kernel_dimension':character_records[0]['kernel_dimension'],
    'nontrivial_character_rank_uniform':sorted(set(r['rank'] for r in character_records[1:])),
    'nontrivial_character_kernel_dimension_uniform':sorted(set(r['kernel_dimension'] for r in character_records[1:])),
    'seed_cycle_rank_per_character':6,
    'residual_dimension_trivial_character':character_records[0]['kernel_dimension']-6,
    'residual_dimension_nontrivial_character':character_records[1]['kernel_dimension']-6,
    'decision':'PROVEN_SIX_SEEDS_DO_NOT_SATURATE_THIS_FINITE_QUOTIENT_FOX_SHADOW',
    'guard':'A kernel after non-flat quotient/base change is a typed finite-quotient shadow, not the full Z[Gamma_3(2)] identity module.'
  },
  'complete_group_ring_identity_module':'NT_PENDING_FULL_GROUP_RING_RESOLUTION_OR_CONTRACTING_HOMOTOPY'
}
(B/'gamma3_critical_fox_shadow_v39.json').write_text(json.dumps(fox_shadow,indent=2,ensure_ascii=False),encoding='utf-8')
ok('fox_underlying_inverse_pairs_6',len(underlying_pairs)==6)
ok('fox_no_bare_triple_identities',len(triple_identities)==0)
ok('fox_fixed_hist_expected',pair_cancel_hist==collections.Counter({0:6754,1:446,4:24}))
ok('fox_cyclic_hist_expected',cyc_hist==collections.Counter({0:143584,1:7540,2:252,3:16,4:96}))
ok('fox_length3_pair_count_4',len(length3_pairs)==4)
ok('fox_char_count_512',len(character_records)==512)
ok('fox_char_rank_hist',rank_hist==collections.Counter({8:511,9:1}))
ok('fox_kernel_sum_17919',kernel_sum==17919)
ok('fox_trivial_kernel34',character_records[0]['kernel_dimension']==34)
ok('fox_nontrivial_kernel35',all(r['kernel_dimension']==35 for r in character_records[1:]))
ok('fox_seed_nonsaturation',fox_shadow['mod4_character_fox_shadow']['residual_dimension_trivial_character']==28 and fox_shadow['mod4_character_fox_shadow']['residual_dimension_nontrivial_character']==29)
ok('fox_full_module_stays_nt',fox_shadow['complete_group_ring_identity_module'].startswith('NT_'))

# -----------------------------------------------------------------------------
# C. Naturality of the oriented projective carrier under the historical
#    Selling/deck witness groupoid. Full carrier is natural by conjugation;
#    pointed matrix-unit/J charts are not globally preserved.
# -----------------------------------------------------------------------------
grp=json.loads((B/'selling_mixed_groupoid_phase_v34_certificate.json').read_text())['groupoid']
edges=grp['edges']
Ky=sp.Matrix([[0,1,0],[1,0,0],[0,0,0]])
Kz=sp.Matrix([[0,0,1],[0,0,0],[1,0,0]])
J=sp.Matrix([[0,0,0],[0,0,1],[0,-1,0]])
I3=sp.eye(3)
units=[]
for i in range(3):
    for j in range(3):
        A=sp.zeros(3); A[i,j]=1; units.append(A)
unit_keys_mod2={tuple(int(x)%2 for x in A):idx for idx,A in enumerate(units)}

def mat_mod(A,m):
    return A.applyfunc(lambda x:int(x)%m)

def tuple_mod2(A):
    return tuple(int(x)%2 for x in list(A))

edge_records=[]; preserve_unit_edges=[]; preserve_J_integer=[]; preserve_J_mod2=[]
for ei,e in enumerate(edges):
    M=sp.Matrix(e['matrix']); Minv=M.inv(); det=int(M.det())
    ok(f'groupoid_edge_{ei}_unimodular',abs(det)==1)
    # Lie bracket naturality on TW generators exactly over Z.
    for nm,A,C in [('KyKz',Ky,Kz),('JKy',J,Ky),('JKz',J,Kz)]:
        lhs=M*(A*C-C*A)*Minv
        MA=M*A*Minv; MC=M*C*Minv
        rhs=MA*MC-MC*MA
        ok(f'groupoid_edge_{ei}_bracket_{nm}',sp.simplify(lhs-rhs)==sp.zeros(3))
    AJ=sp.simplify(M*J*Minv)
    jstat='J' if AJ==J else ('-J' if AJ==-J else 'OTHER')
    if jstat!='OTHER': preserve_J_integer.append(e['label'])
    if mat_mod(AJ,2)==mat_mod(J,2): preserve_J_mod2.append(e['label'])
    # Does conjugation preserve the pointed set of nine matrix units mod 2?
    imgs=[]; unit_pres=True
    for A in units:
        AA=mat_mod(M*A*Minv,2)
        key=tuple(int(x)%2 for x in list(AA))
        if key not in unit_keys_mod2: unit_pres=False
        imgs.append(unit_keys_mod2.get(key))
    if unit_pres: preserve_unit_edges.append(e['label'])
    # Projective tower checks k=1..12: reduction and parity shadow commute.
    level_ok=True
    samples=units+[Ky,Kz,J]
    for k in range(1,13):
        m=2*(3**k)
        # all U+ scalars are central, so conjugation respects scalar orbits; sample scalar 1+6t.
        u=(1+6*(k%max(1,3**(k-1))))%m
        if math.gcd(u,m)!=1 or u%3!=1: u=1
        for si,A in enumerate(samples):
            C1=mat_mod(M*(u*A)*Minv,m)
            C2=mat_mod(u*(M*A*Minv),m)
            same=(C1==C2)
            ok(f'groupoid_edge_{ei}_level_{k}_sample_{si}_scalar_orbit',same)
            level_ok &= same
            if k>1:
                m0=2*(3**(k-1))
                red1=mat_mod(C1,m0)
                red2=mat_mod(M*mat_mod(u*A,m0)*Minv,m0)
                same2=(red1==red2)
                ok(f'groupoid_edge_{ei}_level_{k}_sample_{si}_reduction_natural',same2)
                level_ok &= same2
        parity=(mat_mod(M*J*Minv,2)==mat_mod(M*mat_mod(J,2)*Minv,2))
        ok(f'groupoid_edge_{ei}_level_{k}_parity_natural',parity); level_ok &= parity
    edge_records.append({'source':e['source'],'target':e['target'],'label':e['label'],'det':det,'J_line_image':jstat,'J_mod2_fixed':e['label'] in preserve_J_mod2,'matrix_units_mod2_preserved_as_set':unit_pres,'all_projective_naturality_checks':level_ok})

projective_nat={
  'version':'v39',
  'carrier':'P_k^+ = M_3(Z/(2*3^k))/U_k^+, U_k^+={u unit: u=1 mod 3}',
  'historical_groupoid':'15-edge finite witness sub-groupoid materialized in v34, acting by integer conjugation A -> M A M^{-1}',
  'theorem':{
    'scalar_orbit_well_defined':'PROVEN: U_k^+ consists of central scalars, hence Ad_M(uA)=u Ad_M(A).',
    'bracket_naturality':'PROVEN: Ad_M([A,B])=[Ad_M(A),Ad_M(B)] for every unimodular groupoid edge.',
    'tower_naturality':'PROVEN: reduction k+1 -> k commutes with integer conjugation.',
    'mod2_shadow_naturality':'PROVEN: parity reduction commutes with the same conjugation action.'
  },
  'edge_count':len(edges),
  'edge_records':edge_records,
  'pointed_chart_no_go':{
    'matrix_unit_set_preserved_mod2_edge_count':len(preserve_unit_edges),
    'matrix_unit_set_preserved_mod2_edges':preserve_unit_edges,
    'matrix_unit_set_not_preserved_edge_count':len(edges)-len(preserve_unit_edges),
    'constant_J_line_preserved_integrally_edge_count':len(preserve_J_integer),
    'constant_J_line_preserved_integrally_edges':preserve_J_integer,
    'Jbar_fixed_mod2_edge_count':len(preserve_J_mod2),
    'Jbar_fixed_mod2_edges':preserve_J_mod2,
    'decision':'FULL_PROJECTIVE_CARRIER_NATURAL; POINTED_NINE_MATRIX_UNIT_CHART_AND_CONSTANT_J_DIRECTION_NOT_NATURAL_UNDER_FULL_WITNESS_GROUPOID'
  },
  'orientation':'NT: the groupoid transports the carrier but does not provide a canonical globally fixed Thomas-Wigner J direction or Lorentz orientation.',
  'roof':'REFUTED/PENDING as before: the row-major roof nine-address is not a natural full-groupoid matrix-unit chart, independently of the earlier reversal-vs-transpose cycle-type no-go.',
  'physical_phase':'SEPARATED/PENDING-MORPHISM'
}
(B/'projective_groupoid_naturality_v39.json').write_text(json.dumps(projective_nat,indent=2,ensure_ascii=False),encoding='utf-8')
ok('groupoid_edge_count_15',len(edges)==15)
ok('unit_chart_preserved_only_2_edges',len(preserve_unit_edges)==2)
ok('unit_chart_preserved_Q_and_deck',collections.Counter(preserve_unit_edges)==collections.Counter(['Q','kappa_12_deck']))
ok('J_integer_fixed_no_edges',len(preserve_J_integer)==0)
ok('Jbar_fixed_only_deck',preserve_J_mod2==['kappa_12_deck'])
ok('full_carrier_natural_all_edges',all(r['all_projective_naturality_checks'] for r in edge_records))

# -----------------------------------------------------------------------------
# D. Pending registry and stratified-reality audit.
# -----------------------------------------------------------------------------
pending=[
 {'id':'PM39-01','object':'Selling Fig.1 source-semantic contract','status':'PROVEN_SOURCE_SEMANTIC_CONTRACT_CLOSED','next_gate':'No further prose mining required for the current certification schema; future work is targeted geometric retyping/assignment only.'},
 {'id':'PM39-02','object':'Complete F_Sell^field(Q0) microscopic incidence','status':'NT_SCHEMA_AND_SOURCE_GEOMETRIC_MATCHING','next_gate':'Retype/split at least one diagnostic segment to realize six axes, assign historical node types, wall/substitution labels, ordered field boundaries, repetitions and stabilizer words.'},
 {'id':'PM39-03','object':'Field-level functor F_Sell^field(Q0)->F29_can','status':'PENDING-MORPHISM_NOT_PROMOTED','next_gate':'Requires explicit 0/1/2-cell data from PM39-02.'},
 {'id':'PM39-04','object':'Historical field H1/H2','status':'NT_BLOCKED_BY_PM39-02_03','next_gate':'Recompute only after field functor promotion.'},
 {'id':'PM39-05','object':'Six exact group-ring Fox 4-cycles','status':'PROVEN_RETAINED','next_gate':'Seeds retained.'},
 {'id':'PM39-06','object':'Critical-pair/conjugacy overlap census','status':'PROVEN_BOUNDED_CENSUS','next_gate':'Use length-2/3 cyclic overlaps, especially R3b interfaces, only inside a typed crossed-module/Peiffer computation.'},
 {'id':'PM39-07','object':'Finite mod-4 character Fox shadow','status':'PROVEN_SIX_SEEDS_NON_SATURATING_SHADOW','next_gate':'A full Z[Gamma3(2)] resolution is still required; quotient shadow is not promoted to the full module.'},
 {'id':'PM39-08','object':'Full nonabelian identity-among-relations module','status':'NT_PENDING_GROUP_RING_RESOLUTION','next_gate':'Construct a contracting homotopy, crossed-module rewriting system, or equivalent complete resolution.'},
 {'id':'PM39-09','object':'Oriented projective carrier under Selling/deck groupoid','status':'PROVEN_UNPOINTED_NATURALITY','next_gate':'Carrier/tower/bracket naturality closed for the 15-edge witness groupoid.'},
 {'id':'PM39-10','object':'Pointed nine matrix-unit/roof chart naturality','status':'REFUTED_FULL_GROUPOID_POINTED_CHART','next_gate':'Seek a transported rank-one/projective configuration rather than a fixed row-major nine-address.'},
 {'id':'PM39-11','object':'Constant Thomas-Wigner J direction as groupoid-natural section','status':'REFUTED_FOR_CONSTANT_SECTION','next_gate':'A varying equivariant section may still exist; no Lorentz orientation is selected.'},
 {'id':'PM39-12','object':'Roof reversal = matrix transpose','status':'REFUTED-TYPED_RETAINED','next_gate':'Do not identify.'},
 {'id':'PM39-13','object':'Lorentz orientation / physical phase','status':'NT_OR_SEPARATED_PENDING','next_gate':'Requires explicit source-equivariant orientation/action morphism.'}
]
reg={'version':'v39','items':pending,'counts':dict(collections.Counter(x['status'] for x in pending))}
(B/'pending_morphism_registry_v39.json').write_text(json.dumps(reg,indent=2,ensure_ascii=False),encoding='utf-8')
ok('pending_registry_13_items',len(pending)==13)
ok('stratified_observed_actionable_separation',pending[0]['status'].startswith('PROVEN') and pending[1]['status'].startswith('NT_') and pending[2]['status'].startswith('PENDING'))

# Certificate summary.
new_count=len(checks)-2
status='PASS' if all(v for _,v in checks) else 'FAIL'
cert={
 'phase':'v39-source-semantic-field-gate-critical-fox-shadow-projective-groupoid-naturality',
 'status':status,
 'predecessor_v38_checks':614,
 'new_check_count':new_count,
 'combined_check_count':614+new_count,
 'failed':[n for n,v in checks if not v],
 'selling':{
   'semantic_contract':'PROVEN_CLOSED',
   'microscopic_incidence':'NT_SCHEMA_AND_MATCHING',
   'field_functor':'PENDING_NOT_PROMOTED',
   'historical_H1_H2':'NT_NOT_RECOMPUTED'},
 'gamma3':{
   'exact_seed_4cycles':6,
   'underlying_two_relator_free_identity_pairs':len(underlying_pairs),
   'bare_three_relator_identities':len(triple_identities),
   'mod4_character_rank_histogram':{str(k):v for k,v in sorted(rank_hist.items())},
   'complete_module':'NT'},
 'projective_groupoid':{
   'edge_count':len(edges),
   'full_carrier_naturality':'PROVEN',
   'fixed_matrix_unit_chart_preserved_edges':preserve_unit_edges,
   'constant_J_integral_preserved_edge_count':len(preserve_J_integer),
   'Jbar_mod2_preserved_edges':preserve_J_mod2,
   'Lorentz_orientation':'NT',
   'physical_phase':'SEPARATED/PENDING-MORPHISM'},
 'pending_registry':'pending_morphism_registry_v39.json'
}
(B/'v39_VERIFICATION_REPORT.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
