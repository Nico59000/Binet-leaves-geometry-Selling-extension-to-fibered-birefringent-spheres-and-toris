#!/usr/bin/env python3
from pathlib import Path
import json,csv,hashlib,itertools,collections,re,sys
import sympy as sp
B=Path(__file__).resolve().parent
checks=[]
def ok(name,cond):
    checks.append((name,bool(cond)))
    if not cond: print('FAIL',name,file=sys.stderr)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def dump(name,obj): (B/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,sort_keys=True),encoding='utf-8')
def invtok(t): return t[:-3] if t.endswith('^-1') else t+'^-1'
def invword(w): return tuple(invtok(t) for t in reversed(w))
def reduce_word(w):
    s=[]
    for t in w:
        if s and invtok(t)==s[-1]: s.pop()
        else: s.append(t)
    return tuple(s)
def boundary_cancel(a,b):
    k=0
    for j in range(1,min(len(a),len(b))+1):
        if tuple(a[-j:])==invword(b[:j]): k=j
        else: break
    return k
def grterm(d,g,c):
    if c:
        d[g]=d.get(g,0)+c
        if d[g]==0: del d[g]
def vterm(v,k,g,c):
    d=v.setdefault(k,{}); grterm(d,g,c)
    if not d:v.pop(k,None)
def vsub(A,C):
    out={k:dict(d) for k,d in A.items()}
    for k,d in C.items():
        for g,c in d.items(): vterm(out,k,g,-c)
    return out

# -----------------------------------------------------------------------------
# Parent v43 append-only integrity gate.
# -----------------------------------------------------------------------------
tok=json.loads((B/'v43_ATOMIC_COMPLETION_TOKEN.json').read_text(encoding='utf-8'))
ok('parent_v43_closed',tok.get('publication_state')=='CLOSED/APPEND-ONLY')
ok('parent_v43_atomic',tok.get('atomic_completion') is True)
ok('parent_v43_checks_7566',tok.get('verification',{}).get('combined')==7566)
ok('parent_v43_atomic_gates',all(v=='PASS' for v in tok.get('atomic_conjunction',{}).values()))
entries=[];bad=[]
for line in (B/'v43_SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
    if not line.strip():continue
    h,n=line.split('  ',1);entries.append(n);p=B/n
    if (not p.exists()) or sha(p)!=h:bad.append(n)
ok('parent_v43_ledger_195_entries',len(entries)==195)
ok('parent_v43_ledger_byte_exact',not bad)
ok('parent_v43_pdf_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v43.pdf')==tok['pdf']['sha256'])
ok('parent_v43_tex_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v43.tex')==tok['tex']['sha256'])
ok('parent_LOG0016_exact',sha(B/'LOG_0016.pdf')=='423da98bd897c5dcb2b3556f4c9dc1367a4f12263d247608cfb2ba7a00f2e47b' and (B/'LOG_0016.pdf').stat().st_size==1991876)

# -----------------------------------------------------------------------------
# A. Page-210 overflow repair in v44 source.
# -----------------------------------------------------------------------------
tex=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v44.tex').read_text(encoding='utf-8')
ok('layout44_old_unbreakable_line_removed',r'\\texttt{historical\\_substitution\\_label} et le rôle de frontière' not in tex)
ok('layout44_breakable_axis_path',r'\path{historical_axis_id}' in tex)
ok('layout44_breakable_wall_path',r'\path{historical_wall_label}' in tex)
ok('layout44_breakable_substitution_path',r'\path{historical_substitution_label}' in tex)
ok('layout44_explicit_linebreak',r'\path{historical_wall_label},\\\path{historical_substitution_label}' in tex)
layout={
 'version':'v44','source_page_in_v43':210,
 'defect':'The bullet listing historical_axis_id, historical_wall_label and historical_substitution_label was set as one unbreakable monospace run and crossed the right text border.',
 'repair':'Replace the three long texttt identifiers by breakable path tokens and insert an explicit line break before historical_substitution_label.',
 'source_fix_status':'PROVEN_IN_TEX; FINAL_PDF_RENDER_QA_REQUIRED_BEFORE_PUBLICATION',
 'guard':'No mathematical content or historical status is changed by the typographic repair.'}
dump('v44_page210_overflow_fix.json',layout)

# -----------------------------------------------------------------------------
# B. Selling/BSB source-certified historical C2: explicit acquisition protocol.
# -----------------------------------------------------------------------------
sem=json.loads((B/'selling_fig1_semantic_contract_v39.json').read_text(encoding='utf-8'))
axes=list(csv.DictReader(open(B/'selling_fig1_external_axes_v40.csv',encoding='utf-8')))
centres=list(csv.DictReader(open(B/'selling_fig1_external_centers_v40.csv',encoding='utf-8')))
nodes=list(csv.DictReader(open(B/'selling_fig1_nodes_v40.csv',encoding='utf-8')))
segs=list(csv.DictReader(open(B/'selling_fig1_segments_v40.csv',encoding='utf-8')))
fields=list(csv.DictReader(open(B/'selling_fig1_fields_v40.csv',encoding='utf-8')))
photo=json.loads((B/'selling_fig1_bsb_photo_gate_v38.json').read_text(encoding='utf-8'))
ok('sell44_source_six_axes',len(axes)==6 and all(a['source_type']=='EXTERNAL_BIPARTITE_SYMMETRY_AXIS' for a in axes))
ok('sell44_source_six_centres',len(centres)==6 and all(c['source_type']=='FOURFOLD_SYMMETRY_CENTRE' for c in centres))
ok('sell44_89_degree1_diagnostic',sum(int(r['degree_solid'])==1 for r in nodes)==89)
ok('sell44_axes_still_unmatched',all(a['geometry_assignment_status'].startswith('NT_') for a in axes))
ok('sell44_centres_still_unmatched',all(c['geometry_assignment_status'].startswith('NT_') for c in centres))
ok('sell44_no_historical_node_types',all(r['historical_node_type']=='UNRESOLVED' for r in nodes))
ok('sell44_photo_registration_exists','registration' in photo or 'homography' in str(photo).lower())
rules={r['id']:r for r in sem['source_semantic_rules']}
for rid in ['S39-R1','S39-R2','S39-R3','S39-R4','S39-R5','S39-R6','S39-R8','S39-R9','S39-R10']:
    ok('sell44_semantic_rule_'+rid,rid in rules and str(rules[rid]['status']).startswith('PROVEN'))
protocol={
 'version':'v44','status':'PROVEN_MINIMAL_CERTIFICATION_PROTOCOL; HISTORICAL_ACTION_NOT_YET_MATERIALIZED',
 'source_basis':{
   'semantic_rules':['S39-R1','S39-R2','S39-R3','S39-R4','S39-R5','S39-R6','S39-R8','S39-R9','S39-R10'],
   'primary_plate':'LOG_0016.pdf / Tafel III Fig.1 plus BSB photo witness',
   'photo_registration':'v38 hash-locked homography/cross-registration'},
 'required_axis_binding':[ 
   'For each HAX_i, certify two source-visible fourfold centres HC_{i-1},HC_i on the primary plate and their registered canonical coordinates.',
   'The line through those two certified centres is the only admissible geometric reflection axis; Hough/vectorizer candidates may assist but cannot certify it.',
   'Verify that the reflected source strokes agree with the historical bipartite symmetry statement on a neighbourhood on both sides of the line.'
 ],
 'reflection_action_formula':'For a certified line n.x=d, s_i(x)=x-2 n (n.x-d)/(n.n).',
 'node_action_gate':[ 
   'Apply s_i only after the axis is certified.',
   'Match reflected diagnostic nodes by a unique bijection inside a tolerance derived from v38 registration error plus source-stroke uncertainty, not a hand-chosen tolerance.',
   'Require involutivity pi_i^2=id, degree preservation, and preservation of incident segment layer and source-certified wall/field incidence.',
   'Only if the 89 degree-one set is invariant may it be called an 89-node C2 carrier; otherwise the historical carrier is the certified invariant subset instead.'
 ],
 'six_axis_conjugacy_gate':'The printed sixfold central symmetry may be used to conjugate the six reflections only after the central rotation itself is bound to the raster/source incidence.',
 'current_result':{'axis_bindings':0,'C2_actions_on_89':0,'orbit_census':'NOT_ACTIVATED_FAIL_CLOSED'},
 'carrier_comparisons':{'44|45':'SEPARATED_NOT_TESTED','34|55':'SEPARATED_NOT_TESTED'},
 'guard':'A visually plausible reflection or a degree-preserving permutation is not historical evidence.'}
dump('selling_bsb_C2_action_protocol_v44.json',protocol)

# Data needed for exact historical H1/H2.
required={
 'version':'v44','status':'PROVEN_REQUIREMENT_SCHEMA; RECOMPUTATION_BLOCKED',
 'historical_cell_complex':{
   'C0':'source-certified fissure/crossing/intersection nodes with unique historical IDs',
   'C1':'oriented source-certified boundary arcs with source/target node IDs and historical substitution matrix M_e',
   'C2':'source-certified fields with ordered oriented boundary cycles',
   'quotient':'explicit repeated-field identifications for the finite domain modulo repetition',
   'isotropy':'per-cell stabilizer/self-transformation words for repeated fields and symmetric cells'},
 'local_systems':{
   'Binet_sign_line':'edge transport epsilon_e induced by the historical congruence arrow on the selected time-cone component',
   'Adjoint_soQ':'edge transport Ad(M_e):X -> M_e^{-1} X M_e in compatible bases of so(Q_v)'},
 'cochain_matrices':{
   'd0':'for edge e:v->w, (d0 s)_e = s_w - T_e s_v',
   'd1':'transported signed sum of edge cochains around each ordered field boundary; terminal transport/stabilizer must agree with the repeated-field identification',
   'chain_condition':'d1*d0=0 exactly; face holonomies must equal the declared stabilizer actions'},
 'F29_functor_gate':{
   'objects':'map every historical 0-cell/object to a canonical F29 object',
   'edges':'map each historical substitution edge to an explicit F29 path with matching endpoints',
   'faces':'supply an explicit F29 filling for every historical 2-cell and verify boundary equality'},
 'recompute_only_after':[
   'all historical IDs are total, not diagnostic placeholders',
   'all repeated fields and stabilizers are encoded',
   'all local-system transports are explicit',
   '0/1/2-cell functor compatibility closes'],
 'formula_after_gate':'H^1=ker(d1)/im(d0); H^2=C^2/im(d1) for the closed 2-dimensional carrier, with stabilizer-invariant cochains if the quotient is treated as a complex of groups.',
 'current_missing':{
   'historical_node_assignments':sum(r['historical_node_type']=='UNRESOLVED' for r in nodes),
   'historical_segment_assignments':sum(r['historical_assignment_status'].startswith('NT_') for r in segs),
   'historical_field_assignments':sum(r['field_assignment_status'].startswith('NT_') for r in fields),
   'repetition_pairs':sum(not r['repetition_of'] for r in fields),
   'stabilizer_words':sum(not r['stabilizer_word'] for r in fields)},
 'decision':'H1_H2_NT_NOT_RECOMPUTED'}
dump('selling_h12_recompute_requirements_v44.json',required)
ok('sell44_H12_nodes_missing_all',required['current_missing']['historical_node_assignments']==136)
ok('sell44_H12_fields_missing_all',required['current_missing']['historical_field_assignments']==34)
ok('sell44_H12_stabilizers_missing_all',required['current_missing']['stabilizer_words']==34)

# Preserve strict R5 adapter refutation.
r5=json.loads((B/'selling_gamma_R5_adapter_v43.json').read_text())
ok('sell44_R5_adapter_refuted',r5['strict_incidence_preserving_16_to_16_adapter'].startswith('REFUTED-TYPED'))

# -----------------------------------------------------------------------------
# C. Gamma3: enumerate minimal 3-relator comparison identities beyond adjacent pairwise Squier moves.
# -----------------------------------------------------------------------------
g3=json.loads((B/'gamma3_level2_relative_3cells_v36.json').read_text());cells=g3['cells'];gens=g3['generator_order'];ci={c['name']:i for i,c in enumerate(cells)}
STD={}
for i in range(1,4):
    F=sp.eye(3);F[i-1,i-1]=-1;STD[f'F{i}']=F
for i in range(1,4):
    for j in range(1,4):
        if i!=j:
            E=sp.eye(3);E[i-1,j-1]=2;STD[f'E{i}{j}']=E
def mt(M):return tuple(tuple(int(M[i,j]) for j in range(3)) for i in range(3))
def mm(A,C):return tuple(tuple(sum(A[i][k]*C[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def mi(A):return mt(sp.Matrix(A).inv())
I=mt(sp.eye(3));MAT={g:mt(M) for g,M in STD.items()};MATI={g:mi(MAT[g]) for g in MAT}
def wordmat(w):
    A=I
    for t in w:
        inv=t.endswith('^-1');b=t[:-3] if inv else t;A=mm(A,MATI[b] if inv else MAT[b])
    return A
inst=[]
for i,c in enumerate(cells):
    w=tuple(c['word'])
    for sign,sw in [(1,w),(-1,invword(w))]:
        for sh in range(len(sw)):
            q=sw[:sh]
            inst.append({'id':len(inst),'idx':i,'name':c['name'],'sign':sign,'shift':sh,'word':sw[sh:]+sw[:sh],'qinv':mi(wordmat(q))})
ok('g44_oriented_cyclic_instances_396',len(inst)==396)
def fvec(fs):
    v={}
    for x in fs:vterm(v,x['idx'],x['qinv'],x['sign'])
    return v
# M43 membership.
fc=json.loads((B/'gamma3_group_ring_4cells_v38.json').read_text())
seedpairs=[]
for cyc in fc['cycles']:
    a=ci[cyc['relator_a']];b=ci[cyc['relator_b']];seedpairs.append((a,b,wordmat(tuple(cells[a]['word']))))
seedcoords={i for a,b,_ in seedpairs for i in (a,b)}
square={}
for k,c in enumerate(cells):
    w=tuple(c['word'])
    if len(w)%2==0 and w[:len(w)//2]==w[len(w)//2:]:
        half=w[:len(w)//2];g=wordmat(half)
        if mm(g,g)==I:square[k]=(half,g)
ok('g44_order2_square_relators_18',len(square)==18)
def right_translate(d,r):
    out={}
    for h,c in d.items():grterm(out,mm(h,r),c)
    return out
def in_ideal_1minusg(d,g):
    seen=set()
    for h in list(d):
        if h in seen:continue
        hg=mm(h,g);seen|={h,hg}
        if d.get(h,0)+d.get(hg,0)!=0:return False
    return True
def inM43(v):
    allowed=seedcoords|set(square)
    if any(k not in allowed and d for k,d in v.items()):return False
    for a,b,r in seedpairs:
        if v.get(b,{})!=right_translate(v.get(a,{}),r):return False
    for k,(half,g) in square.items():
        if not in_ideal_1minusg(v.get(k,{}),g):return False
    return True
# Ordered pair products with positive boundary cancellation.
pairs=[];pairres={}
for ai,a in enumerate(inst):
    for bi,b in enumerate(inst):
        if a['idx']==b['idx']:continue
        k=boundary_cancel(a['word'],b['word'])
        if k>=1:
            r=reduce_word(a['word']+b['word']);pairs.append((ai,bi,k,r));pairres[(ai,bi)]=r
ok('g44_ordered_positive_cancel_pairs_7904',len(pairs)==7904)
ok('g44_ordered_cancel_profile',collections.Counter(k for _,_,k,_ in pairs)==collections.Counter({1:7540,2:252,3:16,4:96}))
# Sequential 3-relator factorizations, all three base relators distinct and positive cancellation at both joins.
groups=collections.defaultdict(list)
for ai,bi,k,r in pairs:
    a,b=inst[ai],inst[bi]
    for ci_,c in enumerate(inst):
        if c['idx'] in (a['idx'],b['idx']):continue
        k2=boundary_cancel(r,c['word'])
        if k2<1:continue
        final=reduce_word(r+c['word']);groups[final].append((ai,bi,ci_,k,k2))
triple_count=sum(len(v) for v in groups.values());prof=collections.Counter(len(v) for v in groups.values())
ok('g44_triples_150896',triple_count==150896)
ok('g44_triple_residual_groups_53572',len(groups)==53572)
ok('g44_triple_factor_profile',prof==collections.Counter({1:8192,2:25544,4:16816,8:2996,16:24}))
# Factorization graph generated only by adjacent pairwise Squier replacements.
disconnected=[]
for resid,fs in groups.items():
    if len(fs)<2:continue
    n=len(fs);parent=list(range(n))
    def find(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]];x=parent[x]
        return x
    def union(x,y):
        rx,ry=find(x),find(y)
        if rx!=ry:parent[ry]=rx
    b1=collections.defaultdict(list);b2=collections.defaultdict(list)
    for j,(a,b,c,k,k2) in enumerate(fs):
        b1[(c,pairres[(a,b)])].append(j)
        if (b,c) in pairres:b2[(a,pairres[(b,c)])].append(j)
    for bucket in list(b1.values())+list(b2.values()):
        for j in bucket[1:]:union(bucket[0],j)
    comps=collections.defaultdict(list)
    for j in range(n):comps[find(j)].append(j)
    if len(comps)>1:disconnected.append((resid,fs,list(comps.values())))
ok('g44_pairwise_squier_disconnected_triple_groups_248',len(disconnected)==248)
ok('g44_disconnected_are_two_singletons',all(len(fs)==2 and len(cs)==2 and sorted(map(len,cs))==[1,1] for _,fs,cs in disconnected))
triple_cycles=[]
for resid,fs,cs in disconnected:
    a=fs[cs[0][0]];b=fs[cs[1][0]]
    v=vsub(fvec([inst[i] for i in b[:3]]),fvec([inst[i] for i in a[:3]]))
    triple_cycles.append(v)
for j,v in enumerate(triple_cycles): ok(f'g44_triple_cycle_{j:03d}_absorbed_M43',inM43(v))
ok('g44_all_248_triple_cycles_absorbed',all(inM43(v) for v in triple_cycles))
triple_art={
 'version':'v44','scope':'ordered cyclic three-relator factorizations with three distinct base relators and positive free cancellation at both sequential joins',
 'oriented_cyclic_instances':396,'ordered_positive_cancellation_pairs':7904,'sequential_triples':triple_count,'residual_groups':len(groups),
 'factorization_profile':{str(k):v for k,v in sorted(prof.items())},
 'pairwise_Squier_adjacency':'replace either adjacent pair by another factorization with the same intermediate residual while holding the third/first occurrence fixed',
 'pairwise_Squier_disconnected_groups':248,'component_profile':'248 groups each consist of exactly two isolated triple factorizations',
 'comparison_cycles':248,'absorbed_by_M43':248,'new_group_ring_direction_mod_M43':0,
 'finite_character_residual_min_dimension_inherited':14,
 'decision':'NO_NEW_DIRECTION_FROM_THIS_MINIMAL_THREE_RELATOR_CENSUS; COMPLETE_RESOLUTION_REMAINS_NT',
 'guard':'This exhausts the declared sequential positive-cancellation triple census, not all possible higher identities among relations.'}
dump('gamma3_triple_peiffer_v44.json',triple_art)

# -----------------------------------------------------------------------------
# D. Relations among 18 second syzygies: exact 2-periodic order-two tails.
# -----------------------------------------------------------------------------
rows=[]
for j,(k,(half,g)) in enumerate(sorted(square.items()),1):
    # Check (1-g)(1+g)=0 and reverse product as exact group-ring identities.
    p1={};p2={}
    for a,ca in [(I,1),(g,-1)]:
        for b,cb in [(I,1),(g,1)]:grterm(p1,mm(a,b),ca*cb)
    for a,ca in [(I,1),(g,1)]:
        for b,cb in [(I,1),(g,-1)]:grterm(p2,mm(a,b),ca*cb)
    ok(f'g44_periodic_{j:02d}_g2_identity',mm(g,g)==I)
    ok(f'g44_periodic_{j:02d}_minus_plus_zero',not p1)
    ok(f'g44_periodic_{j:02d}_plus_minus_zero',not p2)
    rows.append({'sector_id':f'P2_{j:02d}','relator':cells[k]['name'],'half_word':' '.join(half),
                 'four_boundary':'(1-g)e_r','five_boundary':'(1+g)f','six_boundary':'(1-g)h',
                 'ker_d5_sector':'R(1-g) h ~= R/R(1+g)','next_kernel':'R(1+g) k ~= R/R(1-g)',
                 'status':'PROVEN_2_PERIODIC_ORDER2_TAIL'})
with open(B/'gamma3_periodic_5_6cells_v44.csv','w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
ok('g44_periodic_18_rows',len(rows)==18)
periodic={
 'version':'v44','ring':'R=Z[Gamma_3(2)]','sectors':18,
 'input_second_syzygies':'y_i=(1+g_i) f_i, g_i^2=1',
 'relations_among_second_syzygies':'ker(partial d5 on the 18 order-two sectors) = direct_sum_i R(1-g_i) h_i',
 'module_form':'direct_sum_i R(1-g_i) ~= direct_sum_i R/R(1+g_i)',
 'new_5cell_attachments':18,'new_6syzygy_directions':'18 canonical (1-g_i)h_i directions in the extended partial resolution',
 'periodicity':'after attaching k_i with d6 k_i=(1-g_i)h_i, the next kernel alternates back to R(1+g_i); each order-two sector has the standard 2-periodic (1-g_i),(1+g_i) tail',
 'cross_sector_relations':'ZERO in the declared 18-sector subcomplex because the f_i/h_i coordinates are distinct',
 'complete_5_6cell_resolution':'NT: no completeness beyond these certified sectors is claimed.'}
dump('gamma3_second_syzygy_relations_v44.json',periodic)

# -----------------------------------------------------------------------------
# E. Singer odd pseudoscalar -> historical orientation line: exact current obstruction.
# -----------------------------------------------------------------------------
v36=(B/'verify_selling_domain_roof_gamma3_v36.py').read_text(encoding='utf-8')
v22=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v22.tex').read_text(encoding='utf-8')
ok('singer44_Gamma_reverses_J','Gamma reverses J' in v36 or r'\\Gamma J\\Gamma^{-1}=-J' in tex)
ok('singer44_no_native_odd_mark','NO_NEW_GAMMA_BREAKING_MARK_FOUND_IN_SOURCE_TEXT_EXTRACTION' in v36)
ok('singer44_separate_pseudoscalar_present','pseudoscalare cubique orienté' in v22 and 'change de signe sous inversion du framing' in v22)
injection={
 'version':'v44','target_orientation_line':'span(J_s), with Gamma acting by -1 on the orientation sign',
 'current_unpointed_selling_deck_marks':'no serialized odd line or source mark; source extraction remains Gamma-symmetric',
 'nonzero_injection_from_current_trivial_scalar_carrier':'REFUTED-TYPED: Gamma-equivariance would require i(x)=-i(x), hence i(x)=0 over characteristic zero',
 'separate_Singer_domain':{'mark':'Omega_sigma(L) epsilon^3','framing_inversion':'acts by -1','status':'SEPARATED/CONDITIONAL-CANDIDATE'},
 'missing_equivariant_data':['a source-native assignment of a line L to each relevant historical object','a Singer framing sigma transported along every historical edge','an explicit identification of Singer-framing inversion with the historical Gamma reflection','path/face compatibility of that assignment'],
 'Singer_to_historical_injection':'NT/PENDING-MORPHISM until all four data are supplied',
 'Lorentz_orientation':'NT','physical_phase':'SEPARATED/PENDING-MORPHISM'}
dump('selling_singer_injection_gate_v44.json',injection)

# Preserve field functor fail closed.
fg=json.loads((B/'selling_field_functor_gate_v43.json').read_text())
ok('sell44_field_functor_not_promoted','PENDING' in str(fg).upper())
dump('selling_field_functor_gate_v44.json',{
 'version':'v44','predecessor':'selling_field_functor_gate_v43.json','zero_cell_gate':'NT','one_cell_gate':'NT','two_cell_gate':'NT',
 'source_certified_axis_action':'NOT_MATERIALIZED','field_functor':'PENDING-MORPHISM','historical_H1_H2':'NT_NOT_RECOMPUTED',
 'requirements_file':'selling_h12_recompute_requirements_v44.json','decision':'FAIL_CLOSED'})

# Registry.
registry={
 'version':'v44','items':[
  {'id':'PM44-01','object':'page 210 right-border overflow','status':'FIXED_IN_TEX_PENDING_FINAL_PDF_QA'},
  {'id':'PM44-02','object':'Selling/BSB six-axis source-certified reflection actions','status':'NT; PROTOCOL_MATERIALIZED'},
  {'id':'PM44-03','object':'89-node historical C2 orbit census','status':'NOT_ACTIVATED_FAIL_CLOSED'},
  {'id':'PM44-04','object':'strict R5-Gamma 16-to-16 critical-diamond incidence adapter','status':'REFUTED-TYPED'},
  {'id':'PM44-05','object':'minimal sequential three-relator Peiffer comparison cycles','status':'248/248 ABSORBED_BY_M43; NO_NEW_DIRECTION'},
  {'id':'PM44-06','object':'relations among 18 second syzygies','status':'PROVEN_18_2_PERIODIC_TAILS'},
  {'id':'PM44-07','object':'complete Gamma3 group-ring resolution','status':'NT; finite-character residual min dimension 14 retained'},
  {'id':'PM44-08','object':'nonzero injection from current unpointed Selling scalar marks to J-sign line','status':'REFUTED-TYPED'},
  {'id':'PM44-09','object':'Singer-framed odd pseudoscalar injection','status':'SEPARATED/NT/PENDING-MORPHISM'},
  {'id':'PM44-10','object':'historical field functor','status':'PENDING-MORPHISM'},
  {'id':'PM44-11','object':'historical H1/H2','status':'NT_NOT_RECOMPUTED'},
  {'id':'PM44-12','object':'roof reversal vs transpose','status':'REFUTED-TYPED'},
  {'id':'PM44-13','object':'physical phase','status':'SEPARATED/PENDING-MORPHISM'}]}
dump('pending_morphism_registry_v44.json',registry)

# Final report (math/verifier only; industrial gates supplied later).
passed=sum(v for _,v in checks); total=len(checks)
report={'version':'v44','predecessor_checks':7566,'new_checks':total,'combined_checks':7566+total,'new_passed':passed,'status':'PASS' if passed==total else 'FAIL','failed':[n for n,v in checks if not v]}
dump('v44_VERIFICATION_REPORT.json',report)
print(json.dumps(report,indent=2))
if passed!=total:sys.exit(1)
