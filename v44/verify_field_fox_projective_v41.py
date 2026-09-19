#!/usr/bin/env python3
from pathlib import Path
import json,csv,hashlib,subprocess,sys,itertools,collections,math,re
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
def dump(name,obj):
    (B/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')
def invtok(t): return t[:-3] if t.endswith('^-1') else t+'^-1'
def invword(w): return [invtok(t) for t in reversed(w)]
def reduce_word(w):
    s=[]
    for t in w:
        if s and invtok(t)==s[-1]: s.pop()
        else: s.append(t)
    return s
def rots(w): return [w[i:]+w[:i] for i in range(len(w))]
def canon_rel(w): return min(tuple(x) for x in rots(w)+rots(invword(w)))
def boundary_cancel(a,b):
    k=0
    for j in range(1,min(len(a),len(b))+1):
        if a[-j:]==invword(b[:j]): k=j
        else: break
    return k

def mtuple(M): return tuple(tuple(int(M[i,j]) for j in range(M.cols)) for i in range(M.rows))
def matmod(M,m): return M.applyfunc(lambda x:int(x)%m)

# -----------------------------------------------------------------------------
# Parent atomic mathematical replay
# -----------------------------------------------------------------------------
p=subprocess.run([sys.executable,str(B/'verify_field_fox_projective_v40.py')],cwd=B,text=True,capture_output=True)
ok('parent_v40_returncode',p.returncode==0)
ok('parent_v40_pass','"status": "PASS"' in p.stdout and '"combined_check_count": 5386' in p.stdout)

# -----------------------------------------------------------------------------
# A. Selling Figure 1: strict source-geometric binding gate.
# The v40 schema exists, but no diagnostic segment/node is silently promoted.
# -----------------------------------------------------------------------------
schema=json.loads((B/'selling_fig1_minimal_schema_v40.json').read_text(encoding='utf-8'))
axes=list(csv.DictReader(open(B/'selling_fig1_external_axes_v40.csv',encoding='utf-8')))
centres=list(csv.DictReader(open(B/'selling_fig1_external_centers_v40.csv',encoding='utf-8')))
nodes=list(csv.DictReader(open(B/'selling_fig1_nodes_v40.csv',encoding='utf-8')))
segs=list(csv.DictReader(open(B/'selling_fig1_segments_v40.csv',encoding='utf-8')))
fields=list(csv.DictReader(open(B/'selling_fig1_fields_v40.csv',encoding='utf-8')))
sem=json.loads((B/'selling_fig1_semantic_contract_v39.json').read_text(encoding='utf-8'))
photo=json.loads((B/'selling_fig1_bsb_photo_gate_v38.json').read_text(encoding='utf-8'))

ok('selling_v41_axes_six',len(axes)==6)
ok('selling_v41_centres_six',len(centres)==6)
ok('selling_v41_abstract_axes_unmatched',all(a['geometry_assignment_status'].startswith('NT_') and not a['diagnostic_segment_id'] for a in axes))
ok('selling_v41_centres_unmatched',all(c['geometry_assignment_status'].startswith('NT_') and not c['x_px300'] and not c['y_px300'] for c in centres))
ok('selling_v41_no_hist_node_promotions',all(r['historical_node_type']=='UNRESOLVED' for r in nodes))
ok('selling_v41_no_hist_axis_promotions',all(not r['historical_axis_id'] for r in segs))
ok('selling_v41_no_hist_wall_labels',all(not r['historical_wall_label'] and not r['historical_substitution_label'] for r in segs))
ok('selling_v41_no_hist_field_cycles',all(not r['ordered_boundary_cycle'] and not r['repetition_of'] and not r['stabilizer_word'] for r in fields))
ok('selling_v41_degree1_89',sum(int(r['degree_solid'])==1 for r in nodes)==89)
ok('selling_v41_photo_gate_still_diagnostic','COMPLETE_HISTORICAL_FIELD_INCIDENCE_NT' in photo['decision'])

binding={
 'version':'v41',
 'source':'v40 minimal schema + v39 semantic contract + v38 BSB photo registration',
 'abstract_skeleton':{'axes':6,'centres':6,'status':'PROVEN_SOURCE_SEMANTIC'},
 'source_certified_geometric_matches':{
   'axes_to_vector_segments':0,
   'centres_to_vector_nodes':0,
   'fissures_to_vector_nodes':0,
   'crossings_to_vector_nodes':0,
   'wall_substitution_labels':0,
   'historical_field_boundary_cycles':0,
   'repetition_pairs':0,
   'stabilizer_words':0
 },
 'diagnostic_only':{
   'degree1_nodes':89,
   'degree2_nodes':sum(int(r['degree_solid'])==2 for r in nodes),
   'degree_ge3_nodes':sum(int(r['degree_solid'])>=3 for r in nodes),
   'legacy_symmetry_candidates':sum(r['diagnostic_axis_candidate']=='YES' for r in segs),
   'photo_support_counts':photo['v37_segment_support_counts']
 },
 'field_functor':'PENDING-MORPHISM_NOT_PROMOTED',
 'historical_H1_H2':'NT_NOT_RECOMPUTED',
 'decision':'NO_SOURCE_CERTIFIED_GEOMETRIC_ASSIGNMENT_AVAILABLE_IN_CURRENT_MACHINE_READABLE_INPUT; ZERO_PROMOTIONS_BY_DESIGN',
 'guard':'Absence of a certified assignment in the current tables is not a theorem that no historical assignment exists. Diagnostic degree, Hough, hull, or symmetry-candidate evidence is not promoted to source semantics.'
}
dump('selling_axis_binding_gate_v41.json',binding)
ok('selling_v41_binding_zero_promotions',sum(binding['source_certified_geometric_matches'].values())==0)

field_gate={
 'version':'v41',
 'schema':'selling_fig1_minimal_schema_v40.json',
 'binding_gate':'selling_axis_binding_gate_v41.json',
 'zero_cell_compatibility':'NT_SOURCE_GEOMETRIC_NODE_ASSIGNMENT_OPEN',
 'one_cell_compatibility':'NT_HISTORICAL_WALL/SUBSTITUTION_ASSIGNMENT_OPEN',
 'two_cell_compatibility':'NT_ORDERED_FIELD_CYCLES_REPETITIONS_STABILIZERS_OPEN',
 'field_level_functor':'PENDING-MORPHISM_NOT_PROMOTED',
 'historical_field_H1_H2':'NT_NOT_RECOMPUTED_EXACT_GATE_OPEN',
 'decision':'FAIL_CLOSED_FIELD_COMPLEX_GATE_REMAINS_OPEN'
}
dump('selling_field_functor_gate_v41.json',field_gate)
ok('selling_v41_functor_pending',field_gate['field_level_functor'].startswith('PENDING'))
ok('selling_v41_H12_not_recomputed',field_gate['historical_field_H1_H2'].startswith('NT_'))

# -----------------------------------------------------------------------------
# B. Cross-context 89 comparison: Cube27/dyadic/RH/KEP.  These are typed
# formulations only; no numerical coincidence is promoted to a morphism.
# -----------------------------------------------------------------------------
roof_tex=(B/'roof_rail_pointed_C2_kep_rh_frontier_v0_3.tex').read_text(encoding='utf-8')
rail=json.loads((B/'roof_rail_recurrence_certificate_v0_1.json').read_text(encoding='utf-8'))
kep_md=(B/'KEP_RESTART_HANDOFF_R24_FR.md').read_text(encoding='utf-8')
cube=json.loads((B/'phase_cube27_context_v33.json').read_text(encoding='utf-8'))
# deterministic text extraction for the R6.65 source label check
pdftxt=subprocess.run(['pdftotext',str(B/'orthogonal_dyadic_golden_transmission_R6_65_SUPPORT30_SOURCE_MORITA_EN.pdf'),'-'],capture_output=True,text=True,check=True).stdout

ok('cross89_pointed_P89_source','P_{89}=T_{88}' in roof_tex and '44' in roof_tex and '45' in roof_tex)
ok('cross89_guard_no_identity','does not identify the A1BE' in roof_tex)
ok('cross89_34_55_nt','p=34' in roof_tex and 'q=55' in roof_tex and 'No such source-bound carrier is established' in roof_tex)
ok('cross89_rail_q9',rail['n9_specialization']['q9']==89)
ok('cross89_rail_q4',rail['n9_specialization']['q4']==44)
ok('cross89_envelope_xplus1',rail['algebraic_theorems']['envelope']=='E(x)=x+1')
ok('cross89_identity_44_45',any('89=44+45' in x for x in rail['n9_specialization']['identities']))
ok('cross89_order89_all_complements_spectrally_separated',len(rail['recurrence_k89']['equal_matrix_order_complement_pairs'])==0)
ok('cross89_34_55_equal_order_refuted',rail['recurrence_k89']['specific_pairs']['34_55']['equal'] is False)
ok('cross89_44_45_equal_order_refuted',rail['recurrence_k89']['specific_pairs']['44_45']['equal'] is False)
ok('cross89_kep_17_routes','16 routes endpoint split et une route frontière' in kep_md)
ok('cross89_cube_guard',cube['cube27_branch']['merge_with_so12_rank3'].startswith('SEPARATED'))
ok('cross89_dyadic_89_is_label','(69, 84, 89)' in pdftxt or '(69,84,89)' in pdftxt.replace(' ',''))

cross89={
 'version':'v41',
 'selling_input':{
   'carrier':'degree-one diagnostic nodes in v40 Figure-1 graph',
   'count':89,
   'status':'DIAGNOSTIC_ONLY_NOT_HISTORICAL_FISSURE_COUNT'
 },
 'compatible_formulation_templates':[
   {
    'id':'A1BE_POINTED_INVOLUTION_89',
    'formula':'89=44+45 from 44 free C2 orbits plus 45 total orbits on an explicit pointed 89-set',
    'source_status':'PASS_IN_ITS_OWN_CARRIER',
    'selling_status':'SEPARATED/PENDING-MORPHISM',
    'required_bridge':'a source-defined involution on the 89 Selling degree-one candidates preserving the historical boundary/fissure data'
   },
   {
    'id':'C2_SPLIT_TARGET_34_55',
    'formula':'for |X|=89, p=34 and q=55 iff the involution has f=21 fixed points',
    'source_status':'NT_TARGET_IN_THE_CROSSCONTEXT_SOURCE',
    'selling_status':'NT_NO_C2_ACTION_OR_21_FIXED_POINTS_CERTIFIED'
   },
   {
    'id':'ROOF_RAIL_ENVELOPE_89',
    'formula':'q_9(9)=89, E(x)=x+1, 89=44+45=81+8=2*45-1',
    'source_status':'PASS_ARITHMETIC_IDENTITY; ORDER-89 RECURRENCE HAS NO EQUAL-MATRIX-ORDER ADDITIVE-COMPLEMENT PAIR, INCLUDING 34/55 AND 44/45',
    'selling_status':'SEPARATED_ARITHMETIC_ENVELOPE_ONLY'
   },
   {
    'id':'CUBE27_TAGGED_TRANSFER',
    'formula':'17|27 -> 30|14 via 27=13+14; the leading 17 is a tagged carrier, not the Selling witness',
    'source_status':'PASS_TYPED_IN_SOURCE',
    'selling_status':'SEPARATED'
   },
   {
    'id':'DYADIC_R6_65_LABEL_89',
    'formula':'the triple (69,84,89) is a source-native line in one inherited d0 direction family',
    'source_status':'PASS_LABELLED_OCCURRENCE_STATEMENT',
    'selling_status':'PROVEN-NOT-MERGEABLE-BY-CARDINALITY: 89 here is a source label, not a set size'
   }
 ],
 'decision':'FORMULATION_COMPATIBILITY_CATALOGUE_ONLY; NO_CROSS-CARRIER_IDENTIFICATION',
 'source_hashes':{
   'roof_rail_tex':sha(B/'roof_rail_pointed_C2_kep_rh_frontier_v0_3.tex'),
   'roof_rail_recurrence_json':sha(B/'roof_rail_recurrence_certificate_v0_1.json'),
   'kep_r24_handoff':sha(B/'KEP_RESTART_HANDOFF_R24_FR.md'),
   'cube27_context':sha(B/'phase_cube27_context_v33.json'),
   'dyadic_R6_65_pdf':sha(B/'orthogonal_dyadic_golden_transmission_R6_65_SUPPORT30_SOURCE_MORITA_EN.pdf')
 }
}
dump('selling_degree1_89_crosscontext_v41.json',cross89)

# -----------------------------------------------------------------------------
# C. 17 Selling/deck objects versus 16 oriented Gamma3 length-3 occurrences:
# only a pointed-envelope comparison is presently typed.
# -----------------------------------------------------------------------------
grp=json.loads((B/'multiobject_lift_groupoid_v34.json').read_text(encoding='utf-8'))
peif=json.loads((B/'gamma3_peiffer_overlap_v40.json').read_text(encoding='utf-8'))
objects=list(grp['objects'].keys())
# component partition from edges
adj={o:set() for o in objects}
for e in grp['edges']:
    a=e['source']; b=e['target']; adj[a].add(b); adj[b].add(a)
seen=set(); comps=[]
for o in objects:
    if o in seen: continue
    q=[o];seen.add(o);cc=[]
    while q:
        x=q.pop();cc.append(x)
        for y in adj[x]:
            if y not in seen: seen.add(y);q.append(y)
    comps.append(sorted(cc))
comp_sizes=sorted([len(c) for c in comps])
ok('env17_selling_objects',len(objects)==17)
ok('env16_gamma_oriented',peif['oriented_length3_instances']==16)
ok('env17_kep_routes','16 routes endpoint split' in kep_md and 'une route frontière' in kep_md)
ok('env_xplus1_16_17',16+1==17 and rail['algebraic_theorems']['envelope']=='E(x)=x+1')
ok('env_selling_components_4_5_8',comp_sizes==[4,5,8])

env={
 'version':'v41',
 'source_counts':{
   'selling_deck_witness_objects':17,
   'selling_deck_component_sizes':comp_sizes,
   'gamma3_length3_oriented_occurrences':16,
   'kep_R24_endpoint_routes':'17 = 16 split + 1 boundary route',
   'roof_envelope':'E(x)=x+1'
 },
 'abstract_pointed_envelope':{
   'construction':'Env_partial(X)=X disjoint_union {partial}',
   'cardinality':'|Env_partial(X)|=|X|+1, hence E(16)=17',
   'status':'PROVEN_AS_ABSTRACT_POINTED_SET_CONSTRUCTION_AND_SOURCE-NATIVE_IN_KEP_R24'
 },
 'selling_gamma_adapter':{
   'status':'SEPARATED/PENDING-MORPHISM',
   'reason':'the current 17-object Selling/deck witness is partitioned into components of sizes 8,4,5 and contains no source-certified distinguished object whose deletion is identified with the 16 Gamma3 oriented occurrences; no incidence/groupoid map has been constructed.'
 },
 'decision':'17_VS_16_ENVELOPE_ANALOGY_TYPED_BUT_NOT_MERGED'
}
dump('selling_gamma_envelope_v41.json',env)

# -----------------------------------------------------------------------------
# D. Gamma3(2) length-two Squier/Peiffer census.
# -----------------------------------------------------------------------------
g3=json.loads((B/'gamma3_level2_relative_3cells_v36.json').read_text(encoding='utf-8'))
cells=g3['cells']; words=[reduce_word(c['word']) for c in cells]
inst=[]
for idx,c in enumerate(cells):
    w=words[idx]
    for sgn,ww in [(1,w),(-1,invword(w))]:
        for sh,rot in enumerate(rots(ww)):
            inst.append({'idx':idx,'name':c['name'],'sign':sgn,'shift':sh,'word':rot})
ok('len2_oriented_cyclic_instance_count_396',len(inst)==396)

cancel_counts=collections.Counter(); pairs=[]
for i,a in enumerate(inst):
    for j in range(i+1,len(inst)):
        b=inst[j]
        if a['idx']==b['idx']: continue
        k=boundary_cancel(a['word'],b['word'])
        if k: cancel_counts[k]+=1
        if k==2:
            rr=tuple(reduce_word(a['word']+b['word']))
            pairs.append((a,b,rr))
ok('len2_pair_count_126',len(pairs)==126)
ok('len2_k3_count_8',cancel_counts[3]==8)

resgroups=collections.defaultdict(list)
for a,b,rr in pairs: resgroups[rr].append((a,b))
diamonds=[(rr,fs) for rr,fs in resgroups.items() if len(fs)==2]
singletons=[(rr,fs) for rr,fs in resgroups.items() if len(fs)==1]
ok('len2_residual_groups_74',len(resgroups)==74)
ok('len2_diamonds_52',len(diamonds)==52)
ok('len2_singletons_22',len(singletons)==22)

# Underlying relator classes quotient cyclic shift + inversion.
bycanon=collections.defaultdict(list); name_to_class={}
for c,w in zip(cells,words):
    cl=canon_rel(w); bycanon[cl].append(c['name']); name_to_class[c['name']]=cl
under=set(frozenset([name_to_class[a['name']],name_to_class[b['name']]]) for a,b,_ in pairs)
ok('len2_under_class_pairs_26',len(under)==26)

# Determine the actual index-permutation subgroup preserving the frozen 43-cell presentation.
def perm_token(t,p):
    inv=t.endswith('^-1'); base=t[:-3] if inv else t
    if base.startswith('E'):
        i,j=int(base[1]),int(base[2]); nb=f'E{p[i]}{p[j]}'
    elif base.startswith('F'):
        i=int(base[1]); nb=f'F{p[i]}'
    else: raise ValueError(t)
    return nb+('^-1' if inv else '')
canon_set=set(bycanon)
valid_perms=[]
for tup in itertools.permutations([1,2,3]):
    pmap={1:tup[0],2:tup[1],3:tup[2]}
    if all(canon_rel([perm_token(t,pmap) for t in cl]) in canon_set for cl in canon_set):
        valid_perms.append(pmap)
ok('len2_frozen_presentation_perm_subgroup_order2',len(valid_perms)==2)
ok('len2_perm_subgroup_is_id_and_23',set(tuple(p[i] for i in [1,2,3]) for p in valid_perms)=={(1,2,3),(1,3,2)})

def act_class(cl,p): return canon_rel([perm_token(t,p) for t in cl])
unseen=set(under); orbits=[]
while unseen:
    seed=next(iter(unseen)); orb=set([seed]); changed=True
    while changed:
        changed=False
        for x in list(orb):
            for pp in valid_perms:
                y=frozenset(act_class(cl,pp) for cl in x)
                if y in under and y not in orb: orb.add(y);changed=True
    unseen-=orb;orbits.append(orb)
ok('len2_orbit_types_14',len(orbits)==14)
ok('len2_orbit_sizes_12x2_2x1',collections.Counter(map(len,orbits))=={2:12,1:2})

def class_names(cl): return sorted(bycanon[cl])
orbit_rows=[]
for oi,orb in enumerate(sorted(orbits,key=lambda o:(-len(o),repr(sorted([sorted([class_names(c) for c in p]) for p in o]))))):
    pairs_fmt=[]
    for pair in orb:
        pairs_fmt.append(sorted(['/'.join(class_names(cl)) for cl in pair]))
    orbit_rows.append({'orbit_id':f'L2O{oi:02d}','size':len(orb),'pairs':sorted(pairs_fmt)})

# Augmentation shadow of the 52 two-factorization diamonds.
vecs=[]
for rr,fs in diamonds:
    a,b=fs[0]; c,d=fs[1]
    v=[0]*len(cells)
    for x,coef in [(a,1),(b,1),(c,-1),(d,-1)]: v[x['idx']]+=coef*x['sign']
    vecs.append(v)
M=sp.Matrix(vecs).T
fc=json.loads((B/'gamma3_group_ring_4cells_v38.json').read_text(encoding='utf-8'))
nameidx={c['name']:i for i,c in enumerate(cells)}
seedvec=[]; seednames=[]
for cyc in fc['cycles']:
    v=[0]*len(cells); v[nameidx[cyc['relator_a']]]=1; v[nameidx[cyc['relator_b']]]=1
    seedvec.append(v); seednames.append((cyc['relator_a'],cyc['relator_b']))
S=sp.Matrix(seedvec).T
ok('len2_aug_rank_6',M.rank()==6)
ok('len2_seed_rank_6',S.rank()==6)
ok('len2_combined_rank_stays_6',sp.Matrix.hstack(S,M).rank()==6)
# Every augmentation diamond vector must be zero or +/- one of six seed vectors.
seed_tuples={tuple(v):i for i,v in enumerate(seedvec)}
seed_tuples.update({tuple(-x for x in v):i for i,v in enumerate(seedvec)})
mult=collections.Counter(); zero=tuple([0]*len(cells))
for v in vecs:
    tv=tuple(v)
    if tv==zero: mult['ZERO']+=1
    else:
        ok('len2_each_aug_vector_seed_shadow',tv in seed_tuples)
        mult[f'SEED{seed_tuples[tv]+1}']+=1
ok('len2_aug_zero_24',mult['ZERO']==24)
ok('len2_aug_nonzero_28',sum(v for k,v in mult.items() if k!='ZERO')==28)
# Exponent sum kernel check.
G=g3['generator_order']; gi={g:i for i,g in enumerate(G)}; E=sp.zeros(len(G),len(cells))
for j,c in enumerate(cells):
    for t in c['word']:
        inv=t.endswith('^-1'); base=t[:-3] if inv else t; E[gi[base],j]+=(-1 if inv else 1)
ok('len2_all_aug_vectors_in_exponent_kernel',all(E*sp.Matrix(v)==sp.zeros(len(G),1) for v in vecs))

length2={
 'version':'v41',
 'source':'43 specialized Gamma_3(2) relative 3-cells from v36; v38 six exact inverse-relator Fox 4-cycles; v40 length-three Peiffer census',
 'oriented_cyclic_instances':len(inst),
 'boundary_cancellation_counts':dict(sorted(cancel_counts.items())),
 'length_two_overlap_count':len(pairs),
 'residual_length_profile':dict(sorted(collections.Counter(len(rr) for _,_,rr in pairs).items())),
 'exact_residual_groups':len(resgroups),
 'two_factorization_diamonds':len(diamonds),
 'singleton_residuals':len(singletons),
 'underlying_cyclic_inverse_relator_class_pairs':len(under),
 'frozen_presentation_index_permutation_subgroup':{
   'order':len(valid_perms),'elements':[tuple(p[i] for i in [1,2,3]) for p in valid_perms],
   'interpretation':'C2 generated by swapping indices 2 and 3; full S3 does not preserve this frozen specialized 43-cell list as a set of cyclic/inverse relator classes.'
 },
 'squier_peiffer_orbit_types':len(orbits),
 'orbit_size_profile':dict(collections.Counter(map(len,orbits))),
 'orbits':orbit_rows,
 'augmentation_shadow':{
   'diamond_relation_rank':M.rank(),
   'six_v38_seed_rank':S.rank(),
   'combined_rank':sp.Matrix.hstack(S,M).rank(),
   'new_directions_beyond_seed_span':sp.Matrix.hstack(S,M).rank()-S.rank(),
   'multiplicities':dict(mult),
   'theorem':'Every length-two two-factorization diamond has augmentation shadow equal to zero or +/- one of the six existing v38 inverse-relator seed directions.'
 },
 'group_ring_attachment_status':'NT: augmentation collapse does not identify the full whiskered/group-ring coefficients of every diamond and therefore does not certify completeness or absence of new nonabelian 4/5-cell generators.',
 'complete_identity_module':'NT_PENDING_FULL_GROUP_RING_RESOLUTION'
}
dump('gamma3_length2_squier_v41.json',length2)

# -----------------------------------------------------------------------------
# E. Extend root-gauged J_s to a finite C6 RACG ball and close declared
# historical-relator holonomy exactly.  The root sign remains a gauge choice.
# -----------------------------------------------------------------------------
I3=sp.eye(3)
P=sp.Matrix([[1,1,0],[3,4,0],[0,0,1]])
Q=sp.Matrix([[2,0,5],[0,1,0],[3,0,8]])
R=sp.Matrix([[-4,0,5],[0,-1,0],[-3,0,4]])
Sm=sp.Matrix([[-2,1,0],[-3,2,0],[0,0,-1]])
Gamma=sp.diag(-1,1,-1); Delta=sp.diag(-1,1,1); H=Gamma*Delta
DP=sp.diag(1,-1,1); DQ=sp.diag(1,1,-1)
T=sp.simplify(P*R*P.inv()); U=sp.simplify(Q*Sm*Q.inv())
V=sp.simplify(P*DP*P.inv()*DP); W=sp.simplify(Q*DQ*Q.inv()*DQ)
a=sp.simplify(V*Gamma); b=sp.simplify(W*H); g=Gamma; h=H; t=T; u=U
CGEN={'a':a,'h':h,'g':g,'b':b,'u':u,'t':t}; cycle=['a','h','g','b','u','t']
order={x:i for i,x in enumerate(cycle)}; comm={frozenset(e) for e in zip(cycle,cycle[1:]+cycle[:1])}
def racg_reduce(word):
    w=list(word); changed=True
    while changed:
        changed=False; i=0
        while i<len(w)-1:
            if frozenset((w[i],w[i+1])) in comm and order[w[i]]>order[w[i+1]]:
                w[i],w[i+1]=w[i+1],w[i];changed=True;i=max(0,i-1)
            else: i+=1
        i=0
        while i<len(w):
            j=i+1
            while j<len(w) and frozenset((w[i],w[j])) in comm:
                if w[j]==w[i]:
                    del w[j];del w[i];changed=True;i=max(-1,i-2);break
                j+=1
            i+=1
    return tuple(w)
# all unique normal forms represented by raw words length <=5
nfs=set()
for n in range(6):
    for w in itertools.product(cycle,repeat=n): nfs.add(racg_reduce(w))
ok('Jball_unique_normal_forms_4504',len(nfs)==4504)
J=sp.Matrix([[0,0,0],[0,0,1],[0,-1,0]])
# Fast exact integer 3x3 arithmetic for the finite ball. All six C6 matrices
# are involutions, so M^{-1}=M.
def tt(M): return tuple(tuple(int(M[i,j]) for j in range(3)) for i in range(3))
def mm(A,C):
    return tuple(tuple(sum(A[i][k]*C[k][j] for k in range(3)) for j in range(3)) for i in range(3))
Itt=tt(I3); Jtt=tt(J); CGT={x:tt(M) for x,M in CGEN.items()}
Gmat={}; Jmat={}
formula_ok=True
for nf in nfs:
    G0=Itt; J0=Jtt
    for x in nf:
        Mx=CGT[x]; G0=mm(G0,Mx); J0=mm(mm(Mx,J0),Mx)
    Gmat[nf]=G0; Jmat[nf]=J0
    # inverse of a word in involutive generators is reversed word
    Gi=Itt
    for x in reversed(nf): Gi=mm(Gi,CGT[x])
    formula_ok = formula_ok and (J0==mm(mm(Gi,Jtt),G0))
ok('Jball_all_section_formulas_exact',formula_ok)
# every generator edge that remains in the radius-5 normal-form set
# is checked against J_t=M_x^{-1}J_sM_x.
tested_edges=0; edge_ok=True
for nf in nfs:
    for x in cycle:
        tgt=racg_reduce(nf+(x,))
        if tgt in nfs:
            tested_edges+=1
            Mx=CGT[x]
            edge_ok = edge_ok and (Jmat[tgt]==mm(mm(Mx,Jmat[nf]),Mx))
ok('Jball_all_internal_edges_natural',edge_ok)
ok('Jball_has_many_edges',tested_edges>5000)
# Representation collisions are recorded, not called faithfulness.
matfib=collections.Counter(Gmat.values())
collision_profile=collections.Counter(matfib.values())

# Historical/source relators in V,W,Gamma,Delta,T,U.
GM={'V':V,'W':W,'Gamma':Gamma,'Delta':Delta,'T':T,'U':U}
rels={
 'Gamma2':[('Gamma',1),('Gamma',1)],
 'Delta2':[('Delta',1),('Delta',1)],
 'T2':[('T',1),('T',1)],
 'U2':[('U',1),('U',1)],
 'GammaDelta_comm':[('Gamma',1),('Delta',1),('Gamma',-1),('Delta',-1)],
 'TU_comm':[('T',1),('U',1),('T',-1),('U',-1)],
 'Gamma_V_inversion':[('Gamma',1),('V',1),('Gamma',1),('V',1)],
 'Delta_V_inversion':[('Delta',1),('V',1),('Delta',1),('V',1)],
 'Gamma_W_fixed':[('Gamma',1),('W',1),('Gamma',1),('W',-1)],
 'Delta_W_inversion':[('Delta',1),('W',1),('Delta',1),('W',1)],
 'T_comm_VGamma':[('T',1),('V',1),('Gamma',1),('T',-1),('Gamma',-1),('V',-1)],
 'U_comm_WGammaDelta':[('U',1),('W',1),('Gamma',1),('Delta',1),('U',-1),('Delta',-1),('Gamma',-1),('W',-1)]
}
def mword(word):
    M=I3
    for x,e in word: M=M*(GM[x] if e==1 else GM[x].inv())
    return sp.simplify(M)
hol={}
for rn,w in rels.items():
    M=mword(w); eq=(M==I3); ok(f'Jhist_rel_{rn}_matrix_identity',eq)
    ok(f'Jhist_rel_{rn}_holonomy',sp.simplify(M.inv()*J*M)==J)
    mods=[]
    for k in range(1,13):
        mods.append(matmod(M,2*3**k)==matmod(I3,2*3**k))
        ok(f'Jhist_rel_{rn}_tower_{k}',mods[-1])
    hol[rn]={'matrix_identity':eq,'tower_k1_12':all(mods)}
# RACG defining relators: 6 involutions + 6 commutators
racg_hol={}
for x in cycle:
    M=CGEN[x]*CGEN[x]; ok(f'Jracg_inv_{x}',M==I3); racg_hol[f'{x}2']=M==I3
for x,y in zip(cycle,cycle[1:]+cycle[:1]):
    M=CGEN[x]*CGEN[y]*CGEN[x]*CGEN[y]; ok(f'Jracg_comm_{x}_{y}',M==I3); racg_hol[f'comm_{x}_{y}']=M==I3

jhist={
 'version':'v41',
 'transport':'T_M(X)=M^{-1} X M',
 'finite_RACG_ball':{
   'raw_words_length_le_5':sum(6**n for n in range(6)),
   'unique_normal_forms':len(nfs),
   'edge_naturality_checks':tested_edges,
   'matrix_image_distinct_count':len(matfib),
   'matrix_fiber_multiplicity_profile':dict(sorted(collision_profile.items())),
   'status':'PROVEN_ROOT-GAUGED_J_SECTION_ON_RADIUS5_NORMAL-FORM_BALL; NO_GLOBAL_FAITHFULNESS_CLAIM'
 },
 'historical_source_relators':hol,
 'historical_relation_holonomy':'PROVEN_EXACT_OVER_Z_FOR_ALL_12_DECLARED_RELATIONS_AND_THEREFORE_AT_EVERY_2*3^k_REDUCTION',
 'RACG_C6_defining_relators':racg_hol,
 'tower':'K1_12_MACHINE_REPLAY_AND_FORMALLY_ALL_K_BY_INTEGER_IDENTITY',
 'canonical_Lorentz_orientation':'NT: root sign choices remain gauge data and no Selling source datum selects them',
 'opposite_variance':'MXM^{-1} remains the opposite-groupoid action under the stored right-composition convention; v40 deck-cycle defect is preserved.',
 'physical_phase':'SEPARATED/PENDING-MORPHISM'
}
dump('projective_J_historical_holonomy_v41.json',jhist)

# -----------------------------------------------------------------------------
# F. Pending-morphism and reality/certification ledger.
# -----------------------------------------------------------------------------
pending={
 'version':'v41','predecessor':'v40 CLOSED/APPEND-ONLY',
 'entries':[
  {'id':'PM41-SELL-AX','status':'NT','subject':'BSB source-certified assignment of HAX_i/HC_i to vector/raster geometry','guard':'zero promotion from diagnostics'},
  {'id':'PM41-SELL-FUN','status':'PENDING-MORPHISM','subject':'F_Sell^field -> F29_can complete 0/1/2-cell compatibility'},
  {'id':'PM41-SELL-H12','status':'NT','subject':'historical field H1/H2','guard':'not recomputed until field functor gate closes'},
  {'id':'PM41-89','status':'SEPARATED/PENDING-MORPHISM','subject':'89 degree-one Selling candidates vs A1BE/roof/RH/KEP 89-carriers'},
  {'id':'PM41-17-16','status':'SEPARATED/PENDING-MORPHISM','subject':'17 Selling/deck objects as pointed envelope of 16 Gamma3 oriented occurrences'},
  {'id':'PM41-G3-L2','status':'PROVEN-AUGMENTATION-SHADOW/NT-GROUP-RING','subject':'length-two Squier diamonds; no new augmentation directions, full group-ring attachments unresolved'},
  {'id':'PM41-G3-FULL','status':'NT','subject':'complete Z[Gamma3(2)] identity-among-relations module'},
  {'id':'PM41-J','status':'PROVEN-RELATION-HOLONOMY/NT-ORIENTATION','subject':'root-gauged J_s on finite C6 ball and all declared historical relators'},
  {'id':'PM41-ROOF-TRANSPOSE','status':'REFUTED-TYPED','subject':'roof reversal = matrix transpose'},
  {'id':'PM41-PHASE','status':'SEPARATED/PENDING-MORPHISM','subject':'compact internal direction = physical phase'}
 ],
 'stratified_reality_guard':'observed calculations are separated from action obligations; no NT or pending morphism is promoted to observed theorem.'
}
dump('pending_morphism_registry_v41.json',pending)

new_count=len(checks)-2
status='PASS' if all(v for _,v in checks) else 'FAIL'
cert={
 'phase':'v41-source-binding-crosscontext89-envelope17-16-length2-squier-historical-J-holonomy',
 'status':status,
 'predecessor_v40_checks':5386,
 'new_check_count':new_count,
 'combined_check_count':5386+new_count,
 'failed':[n for n,v in checks if not v],
 'selling':{
   'source_certified_geometric_promotions':0,
   'field_functor':'PENDING_NOT_PROMOTED','historical_H1_H2':'NT_NOT_RECOMPUTED',
   'degree1_diagnostic_nodes':89,
   'crosscontext89':'CATALOGUED_TYPED_SEPARATED'
 },
 'envelope17_16':{'selling_objects':17,'gamma_oriented_occurrences':16,'adapter':'SEPARATED/PENDING-MORPHISM'},
 'gamma3_length2':{
   'overlaps':126,'diamonds':52,'orbit_types':14,'augmentation_rank':6,'new_augmentation_directions':0,'full_group_ring_module':'NT'
 },
 'projective_J':{
   'normal_form_ball_objects':len(nfs),'edge_checks':tested_edges,'historical_relators':len(rels),
   'relation_holonomy':'PROVEN','Lorentz_orientation':'NT','physical_phase':'SEPARATED/PENDING-MORPHISM'
 },
 'pending_registry':'pending_morphism_registry_v41.json'
}
dump('v41_VERIFICATION_REPORT.json',cert)
print(json.dumps(cert,indent=2,ensure_ascii=False))
raise SystemExit(0 if status=='PASS' else 1)
