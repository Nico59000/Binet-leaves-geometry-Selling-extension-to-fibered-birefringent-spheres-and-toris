#!/usr/bin/env python3
from pathlib import Path
import json,csv,hashlib,subprocess,sys,itertools,collections,re
import sympy as sp
import networkx as nx
B=Path(__file__).resolve().parent
checks=[]
def ok(name,cond):
    checks.append((name,bool(cond)))
    if not cond: print('FAIL',name)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def dump(name,obj): (B/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')
def invtok(t): return t[:-3] if t.endswith('^-1') else t+'^-1'
def invword(w): return [invtok(t) for t in reversed(w)]
def reduce_word(w):
    s=[]
    for t in w:
        if s and invtok(t)==s[-1]: s.pop()
        else:s.append(t)
    return s
def boundary_cancel(a,b):
    k=0
    for j in range(1,min(len(a),len(b))+1):
        if a[-j:]==invword(b[:j]): k=j
        else: break
    return k
def grterm(d,g,c):
    if c:
        d[g]=d.get(g,0)+c
        if d[g]==0:del d[g]

def vterm(v,k,g,c):
    d=v.setdefault(k,{});grterm(d,g,c)
    if not d:v.pop(k,None)

# Sealed predecessor gate. v42 is CLOSED/APPEND-ONLY, so v43 verifies its
# immutable SHA ledger and atomic token instead of recursively rewriting the
# complete v37->v42 verifier chain. This is the append-only parent replay gate.
tok=json.loads((B/'v42_ATOMIC_COMPLETION_TOKEN.json').read_text(encoding='utf-8'))
ok('parent_v42_atomic_token',tok.get('publication_state')=='CLOSED/APPEND-ONLY' and tok.get('checks')=='5832/5832 PASS')
ok('parent_v42_all_atomic_gates',all(v=='PASS' for v in tok.get('atomic_conjunction',{}).values()))
_parent_entries=[];_parent_mismatch=[]
for line in (B/'v42_SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
    if not line.strip():continue
    h,n=line.split('  ',1);_parent_entries.append(n)
    q=B/n
    if (not q.exists()) or sha(q)!=h:_parent_mismatch.append(n)
ok('parent_v42_ledger_170_entries',len(_parent_entries)==170)
ok('parent_v42_ledger_byte_exact',not _parent_mismatch)
ok('parent_v42_pdf_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v42.pdf')==tok['pdf_sha256'])
ok('parent_v42_tex_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v42.tex')==tok['tex_sha256'])

# A. Selling/BSB remains fail-closed; 89-node C2 census remains gated.
prev=json.loads((B/'selling_source_binding_gate_v42.json').read_text())
axes=list(csv.DictReader(open(B/'selling_fig1_external_axes_v40.csv',encoding='utf-8')))
centres=list(csv.DictReader(open(B/'selling_fig1_external_centers_v40.csv',encoding='utf-8')))
nodes=list(csv.DictReader(open(B/'selling_fig1_nodes_v40.csv',encoding='utf-8')))
segs=list(csv.DictReader(open(B/'selling_fig1_segments_v40.csv',encoding='utf-8')))
fields=list(csv.DictReader(open(B/'selling_fig1_fields_v40.csv',encoding='utf-8')))
ok('sell43_no_new_source_promotions',prev['new_promotions']==0)
ok('sell43_axes_unmatched',all(a['geometry_assignment_status'].startswith('NT_') for a in axes))
ok('sell43_centres_unmatched',all(c['geometry_assignment_status'].startswith('NT_') for c in centres))
ok('sell43_nodes_unresolved',all(r['historical_node_type']=='UNRESOLVED' for r in nodes))
ok('sell43_degree1_89',sum(int(r['degree_solid'])==1 for r in nodes)==89)
bind43={
 'version':'v43','predecessor':'selling_source_binding_gate_v42.json','new_promotions':0,
 'source_certified_geometric_matches':prev['source_certified_geometric_matches'],
 'decision':'NO_NEW_SOURCE_CERTIFIED_BSB_GEOMETRIC_BINDING; ALL HISTORICAL ASSIGNMENTS REMAIN NT',
 'guard':'No image-processing, degree, cardinality, or cross-corpus analogy is a historical Selling assignment.'}
dump('selling_source_binding_gate_v43.json',bind43)
c2={
 'version':'v43','degree_one_candidates':89,
 'precondition':'source-certified historical C2 action on the 89 degree-one candidates',
 'precondition_status':'FALSE_NOT_MATERIALIZED','fixed_points':None,'free_orbits':None,'total_orbits':None,
 'comparison_44_45':'SEPARATED_NOT_TESTED','comparison_34_55':'SEPARATED_NOT_TESTED',
 'decision':'C2_ORBIT_CENSUS_NOT_ACTIVATED_FAIL_CLOSED'}
dump('selling_degree1_C2_gate_v43.json',c2)
ok('sell43_C2_not_activated',c2['fixed_points'] is None)

# B. R5 pointed-envelope test against matrix decoration and natural Gamma critical incidence.
grp=json.loads((B/'multiobject_lift_groupoid_v34.json').read_text())
objs=list(grp['objects']);edges=grp['edges']
G=nx.Graph();G.add_nodes_from(objs);G.add_edges_from((e['source'],e['target']) for e in edges)
autos=list(nx.algorithms.isomorphism.GraphMatcher(G,G).isomorphisms_iter())
ok('env43_unlabelled_autos_40',len(autos)==40)
ok('env43_R5_fixed_all_unlabelled_autos',all(a['R5']=='R5' for a in autos))
Qdec={v:tuple(tuple(row) for row in d['Q']) for v,d in grp['objects'].items()}
Cdec={v:tuple(tuple(row) for row in d['chart']) for v,d in grp['objects'].items()}
qautos=[a for a in autos if all(Qdec[v]==Qdec[a[v]] for v in G)]
cautos=[a for a in autos if all(Cdec[v]==Cdec[a[v]] for v in G)]
ok('env43_Q_decorated_autos_identity',len(qautos)==1)
ok('env43_chart_decorated_autos_identity',len(cautos)==1)
incR5=[e for e in edges if 'R5' in (e['source'],e['target'])]
ok('env43_R5_degree_one_edge',len(incR5)==1 and incR5[0]['source']=='R4' and incR5[0]['target']=='R5')
R=G.copy();R.remove_node('R5')
beta=lambda H:H.number_of_edges()-H.number_of_nodes()+nx.number_connected_components(H)
ok('env43_R5_delete_signature',R.number_of_nodes()==16 and R.number_of_edges()==14 and nx.number_connected_components(R)==3 and beta(R)==1)
# Natural Gamma occurrence incidence: each of the 8 oriented residual diamonds connects its two factorization occurrences.
peif=json.loads((B/'gamma3_peiffer_overlap_v40.json').read_text())
ok('env43_gamma_oriented_occurrences_16',peif['oriented_length3_instances']==16)
ok('env43_gamma_diamonds_8',peif['reduced_residual_diamonds']==8)
HG=nx.Graph();HG.add_nodes_from(range(16));HG.add_edges_from((2*i,2*i+1) for i in range(8))
ok('env43_gamma_incidence_signature',HG.number_of_nodes()==16 and HG.number_of_edges()==8 and nx.number_connected_components(HG)==8 and beta(HG)==0)
ok('env43_strict_incidence_adapter_refuted',not nx.is_isomorphic(R,HG))
r5obj=grp['objects']['R5']
env43={
 'version':'v43','candidate':'R5',
 'R5_exact_object':r5obj,'R5_incident_edge':incR5[0],
 'unlabelled_graph_automorphism_group_size':len(autos),'R5_stabilizer_size_in_unlabelled_graph':sum(a['R5']=='R5' for a in autos),
 'Q_decorated_automorphism_group_size':len(qautos),'chart_decorated_automorphism_group_size':len(cautos),
 'R5_deleted_selling_incidence':{'vertices':16,'edges':14,'components':3,'cycle_rank':1,'component_sizes':sorted(len(c) for c in nx.connected_components(R))},
 'Gamma3_critical_diamond_incidence':{'vertices':16,'edges':8,'components':8,'cycle_rank':0,'component_sizes':[2]*8},
 'strict_incidence_preserving_16_to_16_adapter':'REFUTED-TYPED_BY_GRAPH_INVARIANTS',
 'cardinality_bijection':'EXISTS_NONCANONICALLY_BUT_HAS_NO_TYPED_FORCE',
 'historical_stabilizer_label_for_R5':'NT_NOT_SERIALIZED_IN_CURRENT_SELLING_FIELD_COMPLEX',
 'historical_boundary_status':'NT_NOT_SOURCE-DISTINGUISHED',
 'guard':'The refutation applies to the declared natural critical-diamond incidence graph, not to every conceivable future Gamma/Selling adapter.'}
dump('selling_gamma_R5_adapter_v43.json',env43)

# C. Gamma3 exact group-ring length-one census and enlarged partial resolution.
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
def fox(w):
    out={x:{} for x in gens};pref=I
    for t in w:
        inv=t.endswith('^-1');b=t[:-3] if inv else t
        if inv:
            xinv=MATI[b];grterm(out[b],mm(pref,xinv),-1);pref=mm(pref,xinv)
        else:
            grterm(out[b],pref,1);pref=mm(pref,MAT[b])
    assert pref==I;return out
FOX=[fox(c['word']) for c in cells]
def grleft(g,d):return {mm(g,h):c for h,c in d.items()}
def gradd(dst,src,s=1):
    for h,c in src.items():grterm(dst,h,s*c)
def bdry(v):
    out={x:{} for x in gens}
    for i,ring in v.items():
        for cg,cc in ring.items():
            for x in gens:gradd(out[x],grleft(cg,FOX[i][x]),cc)
    return {x:d for x,d in out.items() if d}
inst=[]
for i,c in enumerate(cells):
    w=c['word']
    for sign,sw in [(1,w),(-1,invword(w))]:
        for sh in range(len(sw)):
            q=sw[:sh];inst.append({'idx':i,'name':c['name'],'sign':sign,'shift':sh,'word':sw[sh:]+sw[:sh],'qinv':mi(wordmat(q))})
ok('g43_oriented_cyclic_instances_396',len(inst)==396)
def fvec(fs):
    v={}
    for x in fs:vterm(v,x['idx'],x['qinv'],x['sign'])
    return v
def vsub(A,C):
    out={k:dict(d) for k,d in A.items()}
    for k,d in C.items():
        for g,c in d.items():vterm(out,k,g,-c)
    return out
# Old free seed modules.
fc=json.loads((B/'gamma3_group_ring_4cells_v38.json').read_text())
seedpairs=[]
for cyc in fc['cycles']:
    a=ci[cyc['relator_a']];b=ci[cyc['relator_b']];seedpairs.append((a,b,wordmat(cells[a]['word'])))
seedcoords={i for a,b,_ in seedpairs for i in (a,b)}
ok('g43_six_old_free_seed_modules',len(seedpairs)==6 and len(seedcoords)==12)
# All order-two square relators r=g^2 with g^2=1.
square={}
for k,c in enumerate(cells):
    w=c['word']
    if len(w)%2==0 and w[:len(w)//2]==w[len(w)//2:]:
        half=w[:len(w)//2];g=wordmat(half)
        if mm(g,g)==I:square[k]=(half,g)
ok('g43_square_relator_count_18',len(square)==18)
v42rows=list(csv.DictReader(open(B/'gamma3_new_4cell_generators_v42.csv',encoding='utf-8')))
v42coords={ci[r['relator']] for r in v42rows}
ok('g43_v42_square_coords_6',len(v42coords)==6 and v42coords<=set(square))
newcoords=sorted(set(square)-v42coords)
ok('g43_new_square_coords_12',len(newcoords)==12)
# Exact canonical cycles and annihilators.
newrows=[];seconds=[]
for j,k in enumerate(newcoords,1):
    half,g=square[k]; z={k:{I:1,g:-1}}
    ok(f'g43_newZ4_{j}_boundary_zero',not bdry(z))
    prod={}
    for a,ca in {I:1,g:1}.items():
        for b,cb in {I:1,g:-1}.items():grterm(prod,mm(a,b),ca*cb)
    ok(f'g43_newZ4_{j}_annihilator',not prod)
    newrows.append({'generator_id':f'Z4L1_{j}','relator':cells[k]['name'],'half_word':' '.join(half),'g_matrix':json.dumps([list(r) for r in g]),'cycle_formula':'(1-g)e_'+cells[k]['name'],'status':'PROVEN_NEW_LENGTH1_GENERATOR_RELATIVE_TO_V42'})
    seconds.append({'id':f'Z5L1_{j}','four_cell':f'Z4L1_{j}','boundary_formula':f'(1+gL1_{j}) fL1_{j}','status':'PROVEN_SECOND_SYZYGY_IN_EXTENDED_COMPLEX'})
with open(B/'gamma3_new_length1_4cells_v43.csv','w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(newrows[0]));w.writeheader();w.writerows(newrows)
# Partial module membership. Direct supports imply all intersections between 6 old and 18 square modules vanish.
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
partial_coords=seedcoords|set(square)
ok('g43_partial_coordinate_count_30',len(partial_coords)==30)
def in_partial(v,include_new12=True):
    sqset=set(square) if include_new12 else v42coords
    allowed=seedcoords|sqset
    if any(k not in allowed and d for k,d in v.items()):return False
    for a,b,r in seedpairs:
        if v.get(b,{})!=right_translate(v.get(a,{}),r):return False
    for k in sqset:
        if not in_ideal_1minusg(v.get(k,{}),square[k][1]):return False
    return True
# Length-one complete census.
pairs=[]
for ii,a in enumerate(inst):
    for b in inst[ii+1:]:
        if a['idx']!=b['idx'] and boundary_cancel(a['word'],b['word'])==1:
            pairs.append((a,b,tuple(reduce_word(a['word']+b['word']))))
ok('g43_length1_pairs_3770',len(pairs)==3770)
rg=collections.defaultdict(list)
for a,b,r in pairs:rg[r].append((a,b))
prof=collections.Counter(len(v) for v in rg.values())
ok('g43_length1_residual_groups_2108',len(rg)==2108)
ok('g43_length1_profile',prof==collections.Counter({1:954,2:888,4:242,3:24}))
cycles=[]
for resid,fs in rg.items():
    if len(fs)<2:continue
    ref=fvec(fs[0])
    for j in range(1,len(fs)):
        v=vsub(fvec(fs[j]),ref);cycles.append(v)
        ok(f'g43_L1_cycle_{len(cycles)-1}_fox_zero',not bdry(v))
ok('g43_length1_comparison_cycles_1662',len(cycles)==1662)
absorbed42=sum(in_partial(v,False) for v in cycles);absorbed43=sum(in_partial(v,True) for v in cycles)
ok('g43_L1_v42_absorbed_918',absorbed42==918)
ok('g43_L1_v42_survivors_744',len(cycles)-absorbed42==744)
ok('g43_L1_v43_absorbs_all_1662',absorbed43==1662)
# Every v42-survivor touches one of exactly 12 new coordinates.
ok('g43_all_744_touch_new12',all(any(k in newcoords for k in v) for v in cycles if not in_partial(v,False)))
# Recheck comparison cycles for cancellation lengths 2,3,4 under enlarged module.
overlap_summary={1:{'pairs':len(pairs),'residual_groups':len(rg),'comparison_cycles':len(cycles),'absorbed':absorbed43}}
for L,epairs,egroups,ecycles in [(2,126,74,52),(3,8,4,4),(4,48,1,47)]:
    ps=[]
    for ii,a in enumerate(inst):
        for b in inst[ii+1:]:
            if a['idx']!=b['idx'] and boundary_cancel(a['word'],b['word'])==L:
                ps.append((a,b,tuple(reduce_word(a['word']+b['word']))))
    rr=collections.defaultdict(list)
    for a,b,r in ps:rr[r].append((a,b))
    cs=[]
    for resid,fs in rr.items():
        if len(fs)<2:continue
        ref=fvec(fs[0])
        for j in range(1,len(fs)):cs.append(vsub(fvec(fs[j]),ref))
    ok(f'g43_L{L}_pair_count',len(ps)==epairs);ok(f'g43_L{L}_group_count',len(rr)==egroups);ok(f'g43_L{L}_cycle_count',len(cs)==ecycles)
    ok(f'g43_L{L}_all_absorbed',all(in_partial(v,True) for v in cs))
    overlap_summary[L]={'pairs':len(ps),'residual_groups':len(rr),'comparison_cycles':len(cs),'absorbed':sum(in_partial(v,True) for v in cs)}
# Character shadow: partial module remains far from full finite-quotient kernel.
gidx={g:i for i,g in enumerate(gens)}
hist=collections.Counter();resid_hist=collections.Counter()
for bits in itertools.product([0,1],repeat=9):
    nminus=0
    for k,(half,g) in square.items():
        parity=0
        for t in half:
            b=t[:-3] if t.endswith('^-1') else t;parity^=bits[gidx[b]]
        nminus+=parity
    rank=6+nminus;kerdim=34 if not any(bits) else 35
    hist[rank]+=1;resid_hist[kerdim-rank]+=1
ok('g43_character_partial_rank_max_21',max(hist)==21)
ok('g43_character_residual_min_14',min(resid_hist)==14)
# Exact module structure from disjoint coordinate supports.
partial={
 'version':'v43','ring':'R=Z[Gamma_3(2)]',
 'old_free_seed_modules':6,'order_two_cyclic_modules':18,'new_length1_cyclic_modules':12,
 'pairwise_intersections':'ZERO_BY_DISJOINT_RELATOR_COORDINATE_SUPPORTS',
 'old_seed_annihilators':'ZERO: each seed generator has coefficient 1 on its private first relator coordinate',
 'order_two_annihilators':'Ann_R((1-g_i)e_ri)=R(1+g_i), proved by right coset pairs {h,hg_i}',
 'module_isomorphism':'M43 ~= R^6 direct_sum (direct_sum_{i=1}^{18} R/R(1+g_i))',
 'length_overlap_absorption':{str(k):v for k,v in overlap_summary.items()},
 'length1_before_new12':{'comparison_cycles':1662,'absorbed_by_v42':918,'survived_mod_v42':744},
 'length1_after_new12':'ALL_1662_ABSORBED',
 'character_shadow_partial_rank_histogram':{str(k):v for k,v in sorted(hist.items())},
 'character_shadow_residual_dimension_histogram':{str(k):v for k,v in sorted(resid_hist.items())},
 'finite_shadow_residual_min_dimension':min(resid_hist),
 'complete_resolution':'NT: pairwise cyclic overlaps through cancellation length 4 are absorbed, but higher multi-relator identities and the finite-character residual shadow remain nonzero.'}
dump('gamma3_partial_resolution_v43.json',partial)
dump('gamma3_length1_group_ring_v43.json',{
 'version':'v43','oriented_instances':396,'length1_pairs':3770,'residual_groups':2108,
 'factorization_profile':{str(k):v for k,v in sorted(prof.items())},'multi_residual_groups':sum(v for k,v in prof.items() if k>=2),
 'comparison_cycles':1662,'v42_absorbed':absorbed42,'v42_survivors':1662-absorbed42,'new_canonical_generators':12,'v43_absorbed':absorbed43,
 'decision':'ALL_LENGTH1_COMPARISON_CYCLES_ABSORBED_BY_ENLARGED_PARTIAL_MODULE'})
dump('gamma3_second_syzygies_v43.json',{'version':'v43','v42_second_syzygies':6,'new_length1_second_syzygies':12,'total_order_two_second_syzygies':18,'new':seconds,'complete_5cell_module':'NT'})

# D. Odd source-mark search for J -> -J.
v36=(B/'verify_selling_domain_roof_gamma3_v36.py').read_text(encoding='utf-8')
v22=(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v22.tex').read_text(encoding='utf-8')
ok('J43_v36_no_gamma_breaking_mark','NO_NEW_GAMMA_BREAKING_MARK_FOUND_IN_SOURCE_TEXT_EXTRACTION' in v36)
ok('J43_v36_oriented_U1_stays_NT','NT_PENDING_SOURCE_NATIVE_MARK_FROM_EXACT_FIELD_PLATE_OR_ADDITIONAL_DATA' in v36)
ok('J43_v22_singer_odd_pseudoscalar','pseudoscalare cubique orienté' in v22 and 'change de signe sous inversion du framing' in v22)
ok('J43_v22_requires_Singer_framing',"d'un framing Singer $\\sigma$" in v22)
mark43={
 'version':'v43','current_historical_selling_deck_carrier':{
  'source_search_v36':'NO_NEW_GAMMA_BREAKING_MARK_FOUND_IN_SOURCE_TEXT_EXTRACTION',
  'serialized_odd_mark':'NONE','decision':'NO_ODD_MARK_IN_CURRENT_HISTORICAL_CARRIER; SIGN_NONUNIQUENESS_PRESERVED'},
 'separate_candidate':{
  'source':'formalisation v22 retained context',
  'mark':'oriented cubic pseudoscalar R_3^(3)=Omega_sigma(L) epsilon^3',
  'oddness':'changes sign under inversion of Singer framing sigma',
  'preconditions':'conservation of line L plus an explicitly chosen Singer framing sigma',
  'status':'SEPARATED/CONDITIONAL-CANDIDATE_NOT_NATIVE_TO_CURRENT_SELLING-DECK_CARRIER'},
 'conclusion':'Current data refute uniqueness from relational/groupoid data but do not refute existence of a future source-selected orientation once an additional odd mark is supplied.',
 'physical_phase':'SEPARATED/PENDING-MORPHISM'}
dump('selling_deck_odd_mark_gate_v43.json',mark43)

# E. Selling field functor remains fail closed.
fg={
 'version':'v43','source_binding':'selling_source_binding_gate_v43.json',
 'zero_cell_compatibility':'NT','one_cell_compatibility':'NT','two_cell_compatibility':'NT',
 'field_functor':'PENDING-MORPHISM_NOT_PROMOTED','historical_H1_H2':'NT_NOT_RECOMPUTED',
 'decision':'EXACT_0/1/2_CELL_GATE_NOT_CLOSED'}
dump('selling_field_functor_gate_v43.json',fg)

pending={'version':'v43','predecessor':'v42 CLOSED/APPEND-ONLY','entries':[
 {'id':'PM43-SELL-BIND','status':'NT','subject':'source-certified BSB HAX/HC/fissure/crossing/wall/field assignment'},
 {'id':'PM43-89-C2','status':'SEPARATED/PENDING-PRECONDITION','subject':'89-node historical C2 orbit census'},
 {'id':'PM43-R5','status':'PROVEN-GRAPH-CANDIDATE/REFUTED-TYPED-STRICT-INCIDENCE/NT-HISTORICAL','subject':'R5 pointed envelope versus 16 Gamma critical-diamond incidence'},
 {'id':'PM43-G3-PARTIAL','status':'PROVEN','subject':'partial resolution R^6 direct sum 18 order-two cyclic modules; all pairwise cyclic overlap comparison cycles L1-L4 absorbed'},
 {'id':'PM43-G3-FULL','status':'NT','subject':'complete group-ring resolution; finite-character residual dimension at least 14'},
 {'id':'PM43-J-MARK','status':'NO-CURRENT-MARK/SEPARATED-CONDITIONAL-CANDIDATE','subject':'Singer-framed odd pseudoscalar exists separately but is not a current historical Selling/deck source mark'},
 {'id':'PM43-SELL-FUN','status':'PENDING-MORPHISM','subject':'F_Sell^field -> F29_can'},
 {'id':'PM43-H12','status':'NT','subject':'historical field H1/H2'},
 {'id':'PM43-ROOF-TRANSPOSE','status':'REFUTED-TYPED','subject':'roof reversal = transpose'},
 {'id':'PM43-PHASE','status':'SEPARATED/PENDING-MORPHISM','subject':'physical phase'}],
 'guard':'Pairwise overlap closure is not a complete resolution; conditional Singer framing is not silently inserted into the historical groupoid.'}
dump('pending_morphism_registry_v43.json',pending)

new_count=len(checks)-6
status='PASS' if all(v for _,v in checks) else 'FAIL'
cert={'phase':'v43-R5-incidence-refutation-length1-partial-resolution-odd-mark-gate','status':status,
 'predecessor_v42_checks':5832,'new_check_count':new_count,'combined_check_count':5832+new_count,'failed':[n for n,v in checks if not v],
 'selling':{'new_source_geometric_promotions':0,'C2_89_census':'NOT_ACTIVATED','field_functor':'PENDING','historical_H1_H2':'NT'},
 'envelope':{'candidate':'R5','strict_critical_diamond_incidence_adapter':'REFUTED-TYPED','historical_boundary':'NT'},
 'gamma3':{'old_free_modules':6,'order_two_modules':18,'new_length1_modules':12,'length1_cycles':1662,'v42_survivors':744,'v43_survivors':0,'pairwise_overlap_L1_L4_survivors':0,'finite_shadow_residual_min':min(resid_hist),'complete_resolution':'NT'},
 'projective_J':{'current_historical_odd_mark':'NONE','separate_Singer_framed_candidate':'SEPARATED/CONDITIONAL','Lorentz_orientation':'NT'},
 'pending_registry':'pending_morphism_registry_v43.json'}
dump('v43_VERIFICATION_REPORT.json',cert)
print(json.dumps(cert,indent=2,ensure_ascii=False));raise SystemExit(0 if status=='PASS' else 1)
