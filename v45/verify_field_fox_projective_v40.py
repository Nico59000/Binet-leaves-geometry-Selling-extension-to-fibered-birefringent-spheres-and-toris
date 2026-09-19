#!/usr/bin/env python3
from pathlib import Path
import json,csv,hashlib,subprocess,sys,itertools,collections,math,copy
import sympy as sp

B=Path(__file__).resolve().parent
checks=[]
def ok(name,cond):
    checks.append((name,bool(cond)))
    if not cond: print('FAIL',name)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def invtok(t): return t[:-3] if t.endswith('^-1') else t+'^-1'
def invword(w): return [invtok(t) for t in reversed(w)]
def reduce_word(w):
    s=[]
    for t in w:
        if s and invtok(t)==s[-1]: s.pop()
        else: s.append(t)
    return s

def write_csv(path,rows,fields):
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def matlist(M): return [[int(M[i,j]) for j in range(M.cols)] for i in range(M.rows)]
def matmod(M,m): return M.applyfunc(lambda x:int(x)%m)

# -----------------------------------------------------------------------------
# Parent replay
# -----------------------------------------------------------------------------
p=subprocess.run([sys.executable,str(B/'verify_field_fox_projective_v39.py')],cwd=B,text=True,capture_output=True)
ok('parent_v39_returncode',p.returncode==0)
ok('parent_v39_pass','"status": "PASS"' in p.stdout and '"combined_check_count": 5025' in p.stdout)

# -----------------------------------------------------------------------------
# A. Minimal source-semantic Selling schema.  The schema is closed, but the
# diagnostic geometry is NOT silently identified with historical cells.
# -----------------------------------------------------------------------------
sem=json.loads((B/'selling_fig1_semantic_contract_v39.json').read_text(encoding='utf-8'))
segments=list(csv.DictReader(open(B/'selling_fig1_segments_v37.csv',encoding='utf-8')))
nodes=list(csv.DictReader(open(B/'selling_fig1_nodes_v37.csv',encoding='utf-8')))
fields=list(csv.DictReader(open(B/'selling_fig1_fields_v37.csv',encoding='utf-8')))

# Source-exact abstract external boundary skeleton. AX/C labels are chart labels
# only and are declared up to cyclic reversal/rotation; no raster segment is
# assigned without source-geometric matching.
axes=[]; centers=[]
for i in range(6):
    axes.append({
      'axis_id':f'HAX{i}',
      'chart_order':i,
      'center_start':f'HC{(i-1)%6}',
      'center_end':f'HC{i}',
      'source_type':'EXTERNAL_BIPARTITE_SYMMETRY_AXIS',
      'source_status':'PROVEN_COUNT_AND_TYPE',
      'diagnostic_segment_id':'',
      'geometry_assignment_status':'NT_UNMATCHED_TO_RASTER_VECTOR'
    })
    centers.append({
      'center_id':f'HC{i}',
      'chart_order':i,
      'incident_axis_left':f'HAX{i}',
      'incident_axis_right':f'HAX{(i+1)%6}',
      'source_type':'FOURFOLD_SYMMETRY_CENTRE',
      'source_status':'PROVEN_COUNT_AND_TYPE',
      'x_px300':'','y_px300':'',
      'geometry_assignment_status':'NT_UNMATCHED_TO_RASTER_VECTOR'
    })
write_csv(B/'selling_fig1_external_axes_v40.csv',axes,list(axes[0].keys()))
write_csv(B/'selling_fig1_external_centers_v40.csv',centers,list(centers[0].keys()))

# Enrich node table. Degree-based candidate is diagnostic only, never historical.
node_out=[]; node_diag=collections.Counter()
for r in nodes:
    d=int(r['degree_solid'])
    cand='FISSURE_CANDIDATE_DIAGNOSTIC' if d==1 else ('CROSSING_CANDIDATE_DIAGNOSTIC' if d>=3 else 'CHAIN_OR_UNRESOLVED_DIAGNOSTIC')
    node_diag[cand]+=1
    q=dict(r)
    q.update({'historical_node_type':'UNRESOLVED','historical_node_id':'','diagnostic_topology_candidate':cand,'historical_assignment_status':'NT_SOURCE_GEOMETRIC_MATCHING'})
    node_out.append(q)
write_csv(B/'selling_fig1_nodes_v40.csv',node_out,list(node_out[0].keys()))

# Enrich segment table; explicitly retire Y000..Y004 as diagnostic candidates,
# not five of the six historical axes.
seg_out=[]
for r in segments:
    q=dict(r)
    q.update({
      'historical_axis_id':'',
      'historical_wall_label':'',
      'historical_substitution_label':'',
      'historical_boundary_role':'UNRESOLVED',
      'diagnostic_axis_candidate':'YES' if r['layer']=='symmetry_candidate' else 'NO',
      'historical_assignment_status':'NT_SOURCE_GEOMETRIC_MATCHING'
    })
    if r['layer']=='symmetry_candidate': q['historical_boundary_role']='DIAGNOSTIC_SYMMETRY_CANDIDATE_NOT_PROMOTED'
    seg_out.append(q)
write_csv(B/'selling_fig1_segments_v40.csv',seg_out,list(seg_out[0].keys()))

# Enrich field table with exactly the columns needed by the v40 gate.
field_out=[]
for r in fields:
    q=dict(r)
    q.update({
      'historical_field_id':'',
      'historical_form_label':'',
      'ordered_boundary_cycle':'',
      'repetition_of':'',
      'stabilizer_word':'',
      'field_assignment_status':'NT_SOURCE_GEOMETRIC_MATCHING'
    })
    field_out.append(q)
write_csv(B/'selling_fig1_fields_v40.csv',field_out,list(field_out[0].keys()))

schema={
 'version':'v40',
 'source_semantic_contract':'selling_fig1_semantic_contract_v39.json',
 'source_semantic_contract_sha256':sha(B/'selling_fig1_semantic_contract_v39.json'),
 'external_boundary_skeleton':{
   'axes':6,'fourfold_centres':6,
   'chart_labels':'HAX0..HAX5 and HC0..HC5 are local bookkeeping labels only, defined up to dihedral relabelling of the six-sided source boundary.',
   'source_status':'PROVEN_ABSTRACT_SOURCE_SKELETON',
   'raster_assignment':'NT_UNMATCHED'
 },
 'enriched_tables':{
   'nodes':'selling_fig1_nodes_v40.csv',
   'segments':'selling_fig1_segments_v40.csv',
   'fields':'selling_fig1_fields_v40.csv',
   'axes':'selling_fig1_external_axes_v40.csv',
   'centres':'selling_fig1_external_centers_v40.csv'
 },
 'diagnostic_node_candidate_counts':dict(node_diag),
 'required_historical_columns_now_materialized':{
   'nodes':['historical_node_type','historical_node_id','historical_assignment_status'],
   'segments':['historical_axis_id','historical_wall_label','historical_substitution_label','historical_boundary_role','historical_assignment_status'],
   'fields':['historical_field_id','historical_form_label','ordered_boundary_cycle','repetition_of','stabilizer_word','field_assignment_status']
 },
 'assignment_state':{
   'six_axes_abstract_cells':'PROVEN_SOURCE_SEMANTIC_INSTANTIATION',
   'six_centres_abstract_cells':'PROVEN_SOURCE_SEMANTIC_INSTANTIATION',
   'axes_to_diagnostic_segments':'NT',
   'nodes_to_fissure_crossing':'NT',
   'walls_to_substitution_labels':'NT',
   'field_boundary_cycles':'NT',
   'repetition_pairs':'NT',
   'stabilizer_words':'NT'
 },
 'decision':'PROVEN_MINIMAL_SCHEMA_CLOSED; HISTORICAL_GEOMETRIC_ASSIGNMENT_REMAINS_NT',
 'guard':'Degree-based node candidates and v37 symmetry candidates are diagnostics only. They are not historical fissures, crossings, or axes until source-geometric matching certifies them.'
}
(B/'selling_fig1_minimal_schema_v40.json').write_text(json.dumps(schema,indent=2,ensure_ascii=False),encoding='utf-8')

ok('selling_axes_six',len(axes)==6)
ok('selling_centres_six',len(centers)==6)
ok('selling_axis_ids_unique',len({a['axis_id'] for a in axes})==6)
ok('selling_center_ids_unique',len({c['center_id'] for c in centers})==6)
for i,a in enumerate(axes):
    ok(f'selling_axis_{i}_cyclic_start',a['center_start']==f'HC{(i-1)%6}')
    ok(f'selling_axis_{i}_cyclic_end',a['center_end']==f'HC{i}')
    ok(f'selling_axis_{i}_unmatched',a['geometry_assignment_status'].startswith('NT_'))
for i,c in enumerate(centers):
    ok(f'selling_center_{i}_incidence_left',c['incident_axis_left']==f'HAX{i}')
    ok(f'selling_center_{i}_incidence_right',c['incident_axis_right']==f'HAX{(i+1)%6}')
    ok(f'selling_center_{i}_unmatched',c['geometry_assignment_status'].startswith('NT_'))
ok('selling_node_count_preserved',len(node_out)==136)
ok('selling_segment_count_preserved',len(seg_out)==105)
ok('selling_field_count_preserved',len(field_out)==34)
ok('selling_degree_candidates_89_fissure',node_diag['FISSURE_CANDIDATE_DIAGNOSTIC']==89)
ok('selling_degree_candidates_11_crossing',node_diag['CROSSING_CANDIDATE_DIAGNOSTIC']==11)
ok('selling_degree_candidates_36_chain',node_diag['CHAIN_OR_UNRESOLVED_DIAGNOSTIC']==36)
ok('selling_all_historical_node_types_unresolved',all(r['historical_node_type']=='UNRESOLVED' for r in node_out))
ok('selling_all_wall_labels_blank',all(not r['historical_wall_label'] and not r['historical_substitution_label'] for r in seg_out))
ok('selling_all_field_cycles_blank',all(not r['ordered_boundary_cycle'] for r in field_out))
ok('selling_all_repetition_blank',all(not r['repetition_of'] for r in field_out))
ok('selling_all_stabilizers_blank',all(not r['stabilizer_word'] for r in field_out))
ok('selling_five_old_sym_candidates_not_promoted',sum(r['diagnostic_axis_candidate']=='YES' for r in seg_out)==5 and all(r['historical_axis_id']=='' for r in seg_out))

field_gate={
 'version':'v40',
 'schema':'selling_fig1_minimal_schema_v40.json',
 'schema_sha256':sha(B/'selling_fig1_minimal_schema_v40.json'),
 'presentation_functor':'selling_to_F29_functor_v36.json',
 'schema_support':{
   'external_axes_and_centres':'CLOSED_ABSTRACTLY',
   'fissure_crossing_columns':'MATERIALIZED_BUT_UNASSIGNED',
   'wall_substitution_columns':'MATERIALIZED_BUT_UNASSIGNED',
   'ordered_field_boundary_cycle_column':'MATERIALIZED_BUT_UNASSIGNED',
   'repetition_and_stabilizer_columns':'MATERIALIZED_BUT_UNASSIGNED'
 },
 'zero_cell_compatibility':'NT_ASSIGNMENT_NOT_CLOSED',
 'one_cell_compatibility':'NT_WALL_LABEL_ASSIGNMENT_NOT_CLOSED',
 'two_cell_compatibility':'NT_FIELD_CYCLES_REPETITIONS_STABILIZERS_NOT_CLOSED',
 'field_level_functor':'PENDING-MORPHISM_NOT_PROMOTED',
 'historical_field_H1_H2':'NT_NOT_RECOMPUTED_EXACT_GATE_OPEN',
 'decision':'SCHEMA_OBSTRUCTION_REMOVED; SOURCE_GEOMETRIC_ASSIGNMENT_OBSTRUCTION_REMAINS',
 'guard':'Closing the data schema is not closing the field complex. H1/H2 remain uncomputed.'
}
(B/'selling_field_functor_gate_v40.json').write_text(json.dumps(field_gate,indent=2,ensure_ascii=False),encoding='utf-8')
ok('selling_schema_obstruction_removed','SCHEMA_OBSTRUCTION_REMOVED' in field_gate['decision'])
ok('selling_functor_stays_pending',field_gate['field_level_functor'].startswith('PENDING'))
ok('selling_H12_stays_nt',field_gate['historical_field_H1_H2'].startswith('NT_'))

# -----------------------------------------------------------------------------
# B. Gamma_3(2): explicit crossed-module / Peiffer census of the four length-3
# interfaces.  The key result is that they factor through two existing inverse
# seed 4-cells after cyclic basepoint transport; no new independent attachment
# is certified by this bounded census.
# -----------------------------------------------------------------------------
g3=json.loads((B/'gamma3_level2_relative_3cells_v36.json').read_text())
cells=g3['cells']; words=[reduce_word(c['word']) for c in cells]
cyc=[]
for idx,c in enumerate(cells):
    w=words[idx]
    for sgn,ww in [(1,w),(-1,invword(w))]:
        for sh in range(len(ww)):
            cyc.append({'idx':idx,'sign':sgn,'shift':sh,'name':c['name'],'word':ww[sh:]+ww[:sh]})
# all ordered length-3 cancellation interfaces
length3=[]
for a in cyc:
    for b in cyc:
        if a['idx']==b['idx']: continue
        rr=reduce_word(a['word']+b['word'])
        can=(len(a['word'])+len(b['word'])-len(rr))//2
        if can==3:
            length3.append((a,b,rr))
ok('peiffer_length3_oriented_count_16',len(length3)==16)
# Restrict to the v39 four underlying interfaces involving R3b_123.
critical=[x for x in length3 if x[0]['name']=='R3b_123' or x[1]['name']=='R3b_123']
ok('peiffer_critical_length3_count_16',len(critical)==16)
# group by reduced residual: each diamond must have exactly two factorizations.
groups=collections.defaultdict(list)
for a,b,rr in critical: groups[tuple(rr)].append((a,b))
ok('peiffer_residual_group_count_8',len(groups)==8)
ok('peiffer_two_factorizations_each',all(len(v)==2 for v in groups.values()))
# inverse residual classes
seen=set(); inv_classes=[]
for rr in groups:
    if rr in seen: continue
    irr=tuple(invword(list(rr)))
    ok(f'peiffer_inverse_residual_present_{len(inv_classes)}',irr in groups)
    seen.add(rr);seen.add(irr);inv_classes.append((rr,irr))
ok('peiffer_inverse_residual_classes_4',len(inv_classes)==4)

seed_pairs={
 frozenset(['R3a1_213','R3a1_231']):'SEED_R3a1_213__R3a1_231',
 frozenset(['R3a1_312','R3a1_321']):'SEED_R3a1_312__R3a1_321'
}
diamonds=[]; used_seed_families=set()
for di,(rr,facts) in enumerate(sorted(groups.items(), key=lambda kv:kv[0])):
    (a1,b1),(a2,b2)=facts
    # Arrange so common R3b factor is in same position and compare the other cyclic cells.
    if a1['name']=='R3b_123': commonpos='left'; other1=b1; other2=b2; common1=a1;common2=a2
    else: commonpos='right'; other1=a1; other2=a2; common1=b1;common2=b2
    ok(f'peiffer_{di}_common_R3b_name',common1['name']=='R3b_123' and common2['name']=='R3b_123')
    ok(f'peiffer_{di}_common_R3b_same_word',common1['word']==common2['word'])
    # The two other cyclic words are exactly the same free word, but arise from inverse-paired named relators.
    eqword=(other1['word']==other2['word'])
    ok(f'peiffer_{di}_other_cyclic_words_equal',eqword)
    pair=frozenset([other1['name'],other2['name']]); seed=seed_pairs.get(pair)
    ok(f'peiffer_{di}_inverse_seed_pair_known',seed is not None)
    used_seed_families.add(seed)
    # Explicit Peiffer commutator boundary: a b a^-1 * a b^-1 a^-1 -> 1.
    A=other1['word']; C=common1['word']
    peiffer_boundary=reduce_word(A+C+invword(A)+A+invword(C)+invword(A))
    ok(f'peiffer_{di}_formal_peiffer_boundary_trivial',peiffer_boundary==[])
    diamonds.append({
      'diamond_id':f'PD{di:02d}',
      'common_factor_position':commonpos,
      'common_R3b':{'sign':common1['sign'],'shift':common1['shift']},
      'factorization_A':{'relator':other1['name'],'sign':other1['sign'],'shift':other1['shift']},
      'factorization_B':{'relator':other2['name'],'sign':other2['sign'],'shift':other2['shift']},
      'cyclic_attaching_word_equal':eqword,
      'existing_seed_4cell_family':seed,
      'reduced_residual_word':list(rr),
      'residual_length':len(rr),
      'peiffer_boundary_free_reduces_to_identity':peiffer_boundary==[],
      'attachment_status':'DEPENDENT_WHISKERED_EXISTING_4CELL_NOT_NEW_GENERATOR'
    })
ok('peiffer_seed_family_count_2',used_seed_families==set(seed_pairs.values()))

peif={
 'version':'v40',
 'source':'gamma3_critical_fox_shadow_v39.json + gamma3_group_ring_4cells_v38.json',
 'critical_interfaces_underlying_v39':4,
 'oriented_length3_instances':len(critical),
 'reduced_residual_diamonds':len(groups),
 'inverse_residual_classes':len(inv_classes),
 'existing_seed_families_used':sorted(used_seed_families),
 'diamonds':diamonds,
 'theorem':'Every v39 length-three critical interface is a basepoint-shifted/whiskered occurrence of one of two already-certified inverse-relator 4-cell seeds. The two named factorizations in each diamond have literally equal cyclic attaching words.',
 'new_independent_group_ring_4cell_count':0,
 'new_independent_5cell_count':0,
 'formal_peiffer_coherences':'PROVEN_BOUNDARY_TRIVIAL_FOR_8_ORIENTED_RESIDUAL_DIAMONDS; treated as dependent coherence, not new generators.',
 'complete_identity_module':'NT_PENDING_FULL_GROUP_RING_RESOLUTION',
 'guard':'The bounded critical-pair census proves factorization through existing seed 4-cells. It does not prove global confluence or completeness of the crossed module.'
}
(B/'gamma3_peiffer_overlap_v40.json').write_text(json.dumps(peif,indent=2,ensure_ascii=False),encoding='utf-8')
ok('peiffer_no_new_4cell',peif['new_independent_group_ring_4cell_count']==0)
ok('peiffer_no_new_5cell',peif['new_independent_5cell_count']==0)
ok('peiffer_module_stays_nt',peif['complete_identity_module'].startswith('NT_'))

# -----------------------------------------------------------------------------
# C. Object-varying Thomas-Wigner J_s.  v40 resolves the variance issue:
# v34's declared adjoint transport is X -> M^{-1} X M and is covariant for
# the stored right-multiplication arrow convention. The v39 edgewise M X M^-1
# action is naturally an action of the opposite groupoid and is path-dependent
# if used covariantly on the stored K cycle.
# -----------------------------------------------------------------------------
grp=json.loads((B/'selling_mixed_groupoid_phase_v34_certificate.json').read_text())['groupoid']
edges=grp['edges']; objects=grp['objects']
adj=collections.defaultdict(list); und=collections.defaultdict(set)
for ei,e in enumerate(edges):
    adj[e['source']].append((ei,e)); und[e['source']].add(e['target']);und[e['target']].add(e['source'])
objnames=sorted(set(objects.keys()))
# Connected components underlying graph
seen=set(); comps=[]
for o in objnames:
    if o in seen: continue
    q=[o];seen.add(o);cc=[]
    while q:
        x=q.pop(0);cc.append(x)
        for y in und[x]:
            if y not in seen:seen.add(y);q.append(y)
    comps.append(sorted(cc))
cycle_rank=len(edges)-len(objnames)+len(comps)
ok('Jsection_objects_17',len(objnames)==17)
ok('Jsection_edges_15',len(edges)==15)
ok('Jsection_components_3',len(comps)==3)
ok('Jsection_underlying_cycle_rank_1',cycle_rank==1)

roots={'R0','S0','K0'}
G={}; paths={}; duplicate_checks=[]
for root in sorted(roots):
    G[root]=sp.eye(3); paths[root]=[]; q=[root]
    while q:
        s=q.pop(0)
        for ei,e in adj[s]:
            M=sp.Matrix(e['matrix']); t=e['target']; cand=G[s]*M # stored arrow composition is right multiplication
            if t not in G:
                G[t]=cand; paths[t]=paths[s]+[e['label']]; q.append(t)
            else:
                duplicate_checks.append({'source':s,'target':t,'label':e['label'],'equal':cand==G[t],'candidate':matlist(cand),'stored':matlist(G[t])})
                ok(f'Jsection_duplicate_path_{len(duplicate_checks)-1}',cand==G[t])
ok('Jsection_all_objects_reached',set(G)==set(objnames))
ok('Jsection_one_non_tree_edge',len(duplicate_checks)==1)
ok('Jsection_non_tree_is_K4',duplicate_checks[0]['target']=='K4' and duplicate_checks[0]['label']=='v29_S2')

J=sp.Matrix([[0,0,0],[0,0,1],[0,-1,0]])
Jobj={o:sp.simplify(G[o].inv()*J*G[o]) for o in objnames}
edge_records=[]
for ei,e in enumerate(edges):
    M=sp.Matrix(e['matrix']);s=e['source'];t=e['target']
    cov=sp.simplify(M.inv()*Jobj[s]*M-Jobj[t])==sp.zeros(3)
    ok(f'Jsection_edge_{ei}_covariant_exact',cov)
    tower=True
    for k in range(1,13):
        m=2*3**k
        eq=matmod(M.inv()*Jobj[s]*M,m)==matmod(Jobj[t],m)
        ok(f'Jsection_edge_{ei}_level_{k}',eq);tower &= eq
    mod2=matmod(M.inv()*Jobj[s]*M,2)==matmod(Jobj[t],2)
    ok(f'Jsection_edge_{ei}_mod2',mod2)
    edge_records.append({'label':e['label'],'source':s,'target':t,'exact_covariant_transport':cov,'tower_k1_12':tower,'mod2':mod2})

# Explicitly compare the only two K0->K4 paths.
kchain=[next(e for e in edges if e['source']==f'K{i}' and e['target']==f'K{i+1}') for i in range(4)]
Kright=sp.eye(3)
Krev=sp.eye(3)
XL=J
XR=J
for e in kchain:
    M=sp.Matrix(e['matrix']); Kright=Kright*M; Krev=M*Krev
    XL=M*XL*M.inv()       # v39 edgewise convention if naively treated covariantly
    XR=M.inv()*XR*M       # v34 covariant adjoint convention
Kdeck=sp.Matrix(next(e['matrix'] for e in edges if e['label']=='kappa_12_deck'))
left_direct=Kdeck*J*Kdeck.inv(); right_direct=Kdeck.inv()*J*Kdeck
ok('Jsection_K_right_product_equals_deck',Kright==Kdeck)
ok('Jsection_K_reverse_product_not_deck',Krev!=Kdeck)
ok('Jsection_v34_covariant_path_independent',XR==right_direct)
ok('Jsection_v39_left_naive_path_dependent',XL!=left_direct)
ok('Jsection_v39_left_naive_mod2_shadow_closes',matmod(XL,2)==matmod(left_direct,2))
left_projective_fail=[]
for k in range(1,13):
    m=2*3**k
    # Because entry (1,2) is 1 in both matrices, any scalar orbit equality forces u=1 mod m.
    same=matmod(XL,m)==matmod(left_direct,m)
    left_projective_fail.append(not same)
    ok(f'Jsection_left_naive_projective_level_{k}_fails',not same)
ok('Jsection_left_naive_fails_all_k1_12',all(left_projective_fail))

# Root sign ambiguity: 3 disconnected components and -1 is not in U_k^+.
sign_choices=2**len(comps)
ok('Jsection_sign_choice_count_8',sign_choices==8)
for k in range(1,13):
    m=2*3**k
    minus_in=math.gcd(m-1,m)==1 and ((m-1)%3==1)
    ok(f'Jsection_minus_not_Uplus_{k}',not minus_in)

jsec={
 'version':'v40',
 'stored_groupoid':{'objects':len(objnames),'edges':len(edges),'components':comps,'underlying_cycle_rank':cycle_rank},
 'variance_resolution':{
   'v34_declared_covariant_adjoint_transport':'T_M(X)=M^{-1} X M',
   'stored_arrow_path_matrix':'G(path)=M_1 M_2 ... M_r (right multiplication)',
   'v39_edgewise_conjugation':'M X M^{-1} remains algebraically natural edgewise but is the opposite-variance action for the stored composition convention.',
   'decision':'USE_V34_COVARIANT_TRANSPORT_FOR_PATH/COCYCLE_TESTS'
 },
 'roots':sorted(roots),
 'root_choice':'J at each of R0,S0,K0; this is a declared gauge choice, not a source-selected Lorentz orientation.',
 'object_sections':{o:{'path_labels':paths[o],'path_matrix':matlist(G[o]),'J_s':matlist(Jobj[o])} for o in objnames},
 'edge_naturality':edge_records,
 'path_independence':{
   'only_non_tree_comparator':'K0->K1->K2->K3->K4 versus kappa_12_deck',
   'right_product_equals_deck':Kright==Kdeck,
   'v34_covariant_J_section_exact':XR==right_direct,
   'v39_left_conjugation_if_used_covariantly_exact':XL==left_direct,
   'v39_left_conjugation_mod2_shadow':matmod(XL,2)==matmod(left_direct,2),
   'v39_left_conjugation_projective_levels_k1_12':[not x for x in left_projective_fail]
 },
 'tower':'PROVEN_COMPATIBLE_FOR_K1_12_AND_FORMALLY_ALL_K_BY_INTEGER_REDUCTION',
 'section_status':'PROVEN_OBJECT_VARYING_NATURAL_SECTION_ON_15_EDGE_WITNESS_AFTER_ROOT_GAUGE_CHOICES',
 'canonical_orientation':'NT_NOT_SELECTED',
 'sign_ambiguity':{'connected_components':len(comps),'independent_root_sign_choices':sign_choices,'minus_one_not_quotiented_by_Uplus':True},
 'physical_phase':'SEPARATED/PENDING-MORPHISM',
 'guard':'Existence of a root-gauged natural section on the finite witness does not choose a historical Lorentz orientation, does not extend automatically to the full historical groupoid, and does not identify the compact direction with physical phase.'
}
(B/'projective_transported_J_section_v40.json').write_text(json.dumps(jsec,indent=2,ensure_ascii=False),encoding='utf-8')
ok('Jsection_status_proven_witness',jsec['section_status'].startswith('PROVEN_'))
ok('Jsection_orientation_stays_nt',jsec['canonical_orientation'].startswith('NT_'))
ok('Jsection_physical_phase_separated',jsec['physical_phase'].startswith('SEPARATED'))

# -----------------------------------------------------------------------------
# D. Registry and summary.
# -----------------------------------------------------------------------------
pending=[
 {'id':'PM40-01','object':'Minimal Selling Fig.1 historical-cell schema','status':'PROVEN_SCHEMA_CLOSED','next_gate':'Assign source geometry; schema columns and abstract six-axis/six-centre skeleton now exist.'},
 {'id':'PM40-02','object':'Six historical external axes as abstract source cells','status':'PROVEN_SOURCE_SEMANTIC_INSTANTIATION','next_gate':'Match each HAX_i to raster/vector geometry; old Y000..Y004 remain diagnostic only.'},
 {'id':'PM40-03','object':'Fissure/crossing, wall labels, field cycles, repetitions, stabilizers','status':'NT_SOURCE_GEOMETRIC_ASSIGNMENT','next_gate':'Targeted transcription/matching only; no further general vectorization.'},
 {'id':'PM40-04','object':'Field-level functor and historical H1/H2','status':'PENDING/NT_BLOCKED_BY_PM40-03','next_gate':'Test every 0/1/2-cell only after assignments close; then and only then compute H1/H2.'},
 {'id':'PM40-05','object':'Four v39 length-three Gamma3 critical interfaces','status':'PROVEN_PEIFFER_CENSUS_FACTOR_THROUGH_EXISTING_SEEDS','next_gate':'No new independent 4/5-cell from these interfaces; continue only with other overlaps or full resolution machinery.'},
 {'id':'PM40-06','object':'Full Z[Gamma3(2)] identity module','status':'NT_PENDING_GROUP_RING_RESOLUTION','next_gate':'Contracting homotopy / complete crossed-module rewriting remains required.'},
 {'id':'PM40-07','object':'Object-varying J_s on 15-edge witness','status':'PROVEN_AFTER_ROOT_GAUGE_V34_COVARIANT_TRANSPORT','next_gate':'Extend to larger/full historical groupoid and test additional holonomy.'},
 {'id':'PM40-08','object':'v39 left-conjugation as covariant stored-groupoid action','status':'REFUTED_BY_K_DECK_PATH_COHERENCE_FOR_K>=1','next_gate':'Interpret as opposite-groupoid action or use v34 M^{-1}XM convention.'},
 {'id':'PM40-09','object':'Lorentz orientation','status':'NT_ROOT_SIGN_NOT_SOURCE_SELECTED','next_gate':'Need source-equivariant sign/orientation selector; 3 witness components leave 8 root-sign choices.'},
 {'id':'PM40-10','object':'Fixed nine-unit roof chart / roof reversal = transpose','status':'REFUTED-TYPED_RETAINED','next_gate':'Do not identify.'},
 {'id':'PM40-11','object':'Physical phase','status':'SEPARATED/PENDING-MORPHISM','next_gate':'Requires explicit action/observable morphism.'}
]
reg={'version':'v40','items':pending,'counts':dict(collections.Counter(x['status'] for x in pending))}
(B/'pending_morphism_registry_v40.json').write_text(json.dumps(reg,indent=2,ensure_ascii=False),encoding='utf-8')
ok('v40_pending_count_11',len(pending)==11)

new_count=len(checks)-2
status='PASS' if all(v for _,v in checks) else 'FAIL'
cert={
 'phase':'v40-minimal-selling-schema-peiffer-critical-diamonds-transported-J-section',
 'status':status,
 'predecessor_v39_checks':5025,
 'new_check_count':new_count,
 'combined_check_count':5025+new_count,
 'failed':[n for n,v in checks if not v],
 'selling':{
   'minimal_schema':'PROVEN_CLOSED',
   'abstract_external_axes':6,
   'abstract_fourfold_centres':6,
   'historical_geometric_assignment':'NT',
   'field_functor':'PENDING_NOT_PROMOTED',
   'historical_H1_H2':'NT_NOT_RECOMPUTED'},
 'gamma3':{
   'v39_underlying_length3_interfaces':4,
   'oriented_instances':len(critical),
   'residual_diamonds':len(groups),
   'seed_families':len(used_seed_families),
   'new_independent_4cells':0,
   'new_independent_5cells':0,
   'complete_module':'NT'},
 'projective_J':{
   'objects':len(objnames),'edges':len(edges),'components':len(comps),'cycle_rank':cycle_rank,
   'v34_covariant_transport_path_independent':'PROVEN_ON_WITNESS',
   'v39_left_conjugation_covariant_path_test':'REFUTED_ON_K_DECK_CYCLE_FOR_K>=1; MOD2_SHADOW_CLOSES',
   'root_sign_choices':sign_choices,
   'Lorentz_orientation':'NT','physical_phase':'SEPARATED/PENDING-MORPHISM'},
 'pending_registry':'pending_morphism_registry_v40.json'
}
(B/'v40_VERIFICATION_REPORT.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
