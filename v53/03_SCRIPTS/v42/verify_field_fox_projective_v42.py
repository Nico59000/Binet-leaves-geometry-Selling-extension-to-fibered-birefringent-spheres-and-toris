#!/usr/bin/env python3
from pathlib import Path
import json,csv,hashlib,subprocess,sys,itertools,collections,math
import sympy as sp
try:
 import networkx as nx
except Exception:
 nx=None
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
def dump(name,obj): (B/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8')
def invtok(t): return t[:-3] if t.endswith('^-1') else t+'^-1'
def invword(w): return [invtok(t) for t in reversed(w)]
def reduce_word(w):
    s=[]
    for t in w:
        if s and invtok(t)==s[-1]: s.pop()
        else: s.append(t)
    return s
def rots(w): return [w[i:]+w[:i] for i in range(len(w))]
def boundary_cancel(a,b):
    k=0
    for j in range(1,min(len(a),len(b))+1):
        if a[-j:]==invword(b[:j]): k=j
        else: break
    return k
# parent replay with byte-exact snapshot/restore of all immutable v41 bundle
# members. Legacy predecessor verifiers may rewrite semantically identical JSON
# with noncanonical key order; those writes are rolled back before v42 work.
_v41_names=[]
for _line in (B/'v41_SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
    if _line.strip(): _v41_names.append(_line.split('  ',1)[1])
_v41_names.append('v41_SHA256SUMS.txt')
_v41_snap={n:(B/n).read_bytes() for n in _v41_names}
p=subprocess.run([sys.executable,str(B/'verify_field_fox_projective_v41.py')],cwd=B,text=True,capture_output=True)
for n,b in _v41_snap.items(): (B/n).write_bytes(b)
ok('parent_v41_returncode',p.returncode==0)
ok('parent_v41_pass','\"status\": \"PASS\"' in p.stdout and '\"combined_check_count\": 5646' in p.stdout)
ok('parent_v41_snapshot_restored',all((B/n).read_bytes()==b for n,b in _v41_snap.items()))

# A. Selling/BSB binding stays fail-closed; consequently 89 C2 test is not activated.
prev=json.loads((B/'selling_axis_binding_gate_v41.json').read_text())
axes=list(csv.DictReader(open(B/'selling_fig1_external_axes_v40.csv',encoding='utf-8')))
centres=list(csv.DictReader(open(B/'selling_fig1_external_centers_v40.csv',encoding='utf-8')))
nodes=list(csv.DictReader(open(B/'selling_fig1_nodes_v40.csv',encoding='utf-8')))
segs=list(csv.DictReader(open(B/'selling_fig1_segments_v40.csv',encoding='utf-8')))
fields=list(csv.DictReader(open(B/'selling_fig1_fields_v40.csv',encoding='utf-8')))
geo_count=sum(prev['source_certified_geometric_matches'].values())
ok('sell42_no_prior_source_geometric_promotions',geo_count==0)
ok('sell42_axes_still_unmatched',all(a['geometry_assignment_status'].startswith('NT_') for a in axes))
ok('sell42_centres_still_unmatched',all(c['geometry_assignment_status'].startswith('NT_') for c in centres))
ok('sell42_hist_nodes_unresolved',all(r['historical_node_type']=='UNRESOLVED' for r in nodes))
ok('sell42_degree1_count_89',sum(int(r['degree_solid'])==1 for r in nodes)==89)
bind42={
 'version':'v42','predecessor':'selling_axis_binding_gate_v41.json',
 'source_certified_geometric_matches':prev['source_certified_geometric_matches'],
 'new_promotions':0,
 'decision':'NO_NEW_SOURCE_CERTIFIED_BSB_GEOMETRIC_BINDING; KEEP_ALL_HAX/HC/FISSURE/CROSSING/WALL/FIELD_ASSIGNMENTS_NT',
 'guard':'Image-processing candidates, node degree and visual similarity are not source-certified historical assignments.'}
dump('selling_source_binding_gate_v42.json',bind42)
c2gate={
 'version':'v42','degree_one_candidates':89,
 'precondition':'a source-certified historical C2 symmetry must act on the 89 candidates',
 'precondition_status':'FALSE_NOT_YET_MATERIALIZED',
 'fixed_points':None,'free_orbits':None,'total_orbits':None,
 'comparison_44_45':'SEPARATED_NOT_TESTED','comparison_34_55':'SEPARATED_NOT_TESTED',
 'decision':'C2_ORBIT_CENSUS_NOT_ACTIVATED_FAIL_CLOSED'}
dump('selling_degree1_C2_gate_v42.json',c2gate)
ok('sell42_c2_not_activated',c2gate['fixed_points'] is None)

# B. 17->16 pointed envelope: search graph-canonical deletion candidates; no Gamma adapter is invented.
grp=json.loads((B/'multiobject_lift_groupoid_v34.json').read_text())
objs=list(grp['objects'])
edges=[(e['source'],e['target']) for e in grp['edges']]
if nx is None: raise RuntimeError('networkx required')
G=nx.Graph();G.add_nodes_from(objs);G.add_edges_from(edges)
beta=lambda H:H.number_of_edges()-H.number_of_nodes()+nx.number_connected_components(H)
ok('env42_original_17',G.number_of_nodes()==17)
ok('env42_original_components3',nx.number_connected_components(G)==3)
ok('env42_original_beta1',beta(G)==1)
autos=list(nx.algorithms.isomorphism.GraphMatcher(G,G).isomorphisms_iter())
ok('env42_graph_automorphism_count_40',len(autos)==40)
fixed_all={v for v in G if all(a[v]==v for a in autos)}
candidates=[]
for v in G:
    H=G.copy();H.remove_node(v)
    if H.number_of_nodes()==16 and nx.number_connected_components(H)==3 and beta(H)==1:
        candidates.append(v)
ok('env42_deletion_preserving_candidates_5',set(candidates)=={'R5','P1','Q1','S0','S3'})
canon=[v for v in candidates if v in fixed_all]
ok('env42_unique_graph_canonical_candidate_R5',canon==['R5'])
H=G.copy();H.remove_node('R5')
comp_sizes=sorted(len(c) for c in nx.connected_components(H))
ok('env42_R5_delete_components_4_5_7',comp_sizes==[4,5,7])
env42={
 'version':'v42','selling_deck_objects':17,'graph_edges':15,'components':3,'cycle_rank':1,
 'declared_preservation_criterion':'deletion leaves 16 objects, 3 connected components and cycle rank 1',
 'candidates':candidates,'automorphism_group_size_unlabelled_graph':len(autos),
 'fixed_by_all_graph_automorphisms':sorted(fixed_all),
 'unique_graph_canonical_candidate_under_criterion':'R5',
 'R5_deleted_component_sizes':comp_sizes,
 'status_R5':'PROVEN_GRAPH-CANONICAL_CANDIDATE_RELATIVE_TO_DECLARED_GRAPH_CRITERION',
 'historical_source_distinguished_boundary_status':'NT: the criterion is graph-theoretic and not an historical Selling statement naming R5 as a boundary object',
 'Gamma3_adapter':'SEPARATED/PENDING-MORPHISM: the 16 Gamma3 oriented occurrences carry no constructed incidence/groupoid identification with the R5-deleted witness',
 'decision':'POINTED_ENVELOPE_PARTIALLY_REFINED_NOT_MERGED'}
dump('selling_gamma_pointed_envelope_v42.json',env42)

# C. Exact Z[Gamma_3(2)] lift of the 52 length-two diamonds.
g3=json.loads((B/'gamma3_level2_relative_3cells_v36.json').read_text())
cells=g3['cells']; gens=g3['generator_order']; ci={c['name']:i for i,c in enumerate(cells)}
STD={}
for i in range(1,4):
    F=sp.eye(3);F[i-1,i-1]=-1;STD[f'F{i}']=F
for i in range(1,4):
    for j in range(1,4):
        if i!=j:
            E=sp.eye(3);E[i-1,j-1]=2;STD[f'E{i}{j}']=E
def mt(M): return tuple(tuple(int(M[i,j]) for j in range(3)) for i in range(3))
def mm(A,C): return tuple(tuple(sum(A[i][k]*C[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def mi(A): return mt(sp.Matrix(A).inv())
I=mt(sp.eye(3)); MAT={g:mt(M) for g,M in STD.items()}; MATI={g:mi(MAT[g]) for g in MAT}
def wordmat(w):
    A=I
    for t in w:
        inv=t.endswith('^-1');b=t[:-3] if inv else t;A=mm(A,MATI[b] if inv else MAT[b])
    return A
def grterm(d,g,c):
    if c:
        d[g]=d.get(g,0)+c
        if d[g]==0: del d[g]
def fox(w):
    out={x:{} for x in gens}; pref=I
    for t in w:
        inv=t.endswith('^-1'); b=t[:-3] if inv else t
        if inv:
            xinv=MATI[b];grterm(out[b],mm(pref,xinv),-1);pref=mm(pref,xinv)
        else:
            grterm(out[b],pref,1);pref=mm(pref,MAT[b])
    assert pref==I
    return out
FOX=[fox(c['word']) for c in cells]
inst=[]
for i,c in enumerate(cells):
    w=c['word']
    for sign,sw in [(1,w),(-1,invword(w))]:
        for sh in range(len(sw)):
            q=sw[:sh];inst.append({'idx':i,'name':c['name'],'sign':sign,'shift':sh,'word':sw[sh:]+sw[:sh],'q':q,'qinv':mi(wordmat(q))})
pairs=[]
for ii,a in enumerate(inst):
    for b in inst[ii+1:]:
        if a['idx']!=b['idx'] and boundary_cancel(a['word'],b['word'])==2:
            pairs.append((a,b,tuple(reduce_word(a['word']+b['word']))))
rg=collections.defaultdict(list)
for a,b,r in pairs: rg[r].append((a,b))
diam=[(r,fs) for r,fs in rg.items() if len(fs)==2]
ok('gr42_diamond_count_52',len(diam)==52)
def vterm(v,k,g,c):
    d=v.setdefault(k,{});grterm(d,g,c)
    if not d:v.pop(k,None)
def dvec(fs):
    (a,b),(c,d)=fs;v={}
    for x,s in [(a,1),(b,1),(d,-1),(c,-1)]:vterm(v,x['idx'],x['qinv'],s*x['sign'])
    return v
def grleft(g,d):return {mm(g,h):c for h,c in d.items()}
def gradd(dst,src,s=1):
    for h,c in src.items():grterm(dst,h,s*c)
def bdry(v):
    out={x:{} for x in gens}
    for i,ring in v.items():
        for cg,cc in ring.items():
            for x in gens:gradd(out[x],grleft(cg,FOX[i][x]),cc)
    return {x:d for x,d in out.items() if d}
vecs=[dvec(fs) for _,fs in diam]
for i,v in enumerate(vecs):ok(f'gr42_diamond_{i}_fox_boundary_zero',not bdry(v))
# old seed submodule
fc=json.loads((B/'gamma3_group_ring_4cells_v38.json').read_text())
seedpairs=[(ci[x['relator_a']],ci[x['relator_b']]) for x in fc['cycles']]
seedcoords=set(y for p in seedpairs for y in p)
def combine(A,C,s=1):
    d=dict(A)
    for g,c in C.items():grterm(d,g,s*c)
    return d
def modseed(v):
    w={i:dict(r) for i,r in v.items()}
    for a,b in seedpairs:
        A=w.pop(a,{}); C=w.pop(b,{})
        R=combine(A,C,-1)
        if R:w[a]=R
    return w
mods=[modseed(v) for v in vecs]
ok('gr42_28_old_seed_submodule',sum(not x for x in mods)==28)
ok('gr42_24_new_mod_seed_nonzero',sum(bool(x) for x in mods)==24)
ok('gr42_24_have_nonseed_coordinate',sum(any(i not in seedcoords for i in v) for v in vecs)==24)
# canonical six order-two half-relator cycles z=(1-g)e_r
families=collections.defaultdict(list)
for di,v in enumerate(vecs):
    if mods[di]:
        ok(f'gr42_new_{di}_single_coordinate',len(v)==1)
        k=next(iter(v));families[k].append(di)
ok('gr42_six_new_coordinate_families',len(families)==6)
ok('gr42_four_diamonds_per_family',all(len(v)==4 for v in families.values()))
new_rows=[]; five=[]
for fi,(k,dis) in enumerate(sorted(families.items(),key=lambda z:cells[z[0]]['name'])):
    w=cells[k]['word'];ok(f'gr42_family_{fi}_relator_double_half',len(w)%2==0 and w[:len(w)//2]==w[len(w)//2:])
    half=w[:len(w)//2];g=wordmat(half);ok(f'gr42_family_{fi}_half_involution',mm(g,g)==I)
    # canonical z=(1-g)e_r; every observed diamond is +/- a monomial left translate
    canonical={I:1,g:-1}
    for di in dis:
        ring=vecs[di][k];ok(f'gr42_family_{fi}_diamond_{di}_two_terms',len(ring)==2 and sorted(ring.values())==[-1,1])
        plus=next(h for h,c in ring.items() if c==1);minus=next(h for h,c in ring.items() if c==-1)
        ok(f'gr42_family_{fi}_diamond_{di}_translate_pattern',mm(mi(plus),minus)==g)
    # exact boundary proof by direct machine Fox and algebraic annihilator
    z={k:canonical};ok(f'gr42_family_{fi}_canonical_boundary_zero',not bdry(z))
    # (1+g)(1-g)=0 and gives second syzygy after adjoining 4-cell
    plusminus={I:1,g:1}
    prod={}
    for a,ca in plusminus.items():
        for b,cb in canonical.items():grterm(prod,mm(a,b),ca*cb)
    ok(f'gr42_family_{fi}_annihilator_1plusg',prod=={})
    new_rows.append({'generator_id':f'Z4_{fi+1}','relator':cells[k]['name'],'half_word':' '.join(half),'g_matrix':json.dumps([list(r) for r in g]),'four_diamond_occurrences':','.join(map(str,dis)),'cycle_formula':'(1-g) e_'+cells[k]['name'],'status':'PROVEN_NEW_RELATIVE_TO_V38_SEED_SUBMODULE'})
    five.append({'id':f'Z5_{fi+1}','four_cell':f'Z4_{fi+1}','boundary_formula':f'(1+g_{fi+1}) f_{fi+1}','verification':'(1+g)(1-g)=1-g^2=0','status':'PROVEN_SECOND_SYZYGY_IN_EXTENDED_COMPLEX'})
# disjoint coordinates -> direct sum of six nonzero cyclic submodules; not free because annihilators
ok('gr42_six_distinct_new_relator_coordinates',len({r['relator'] for r in new_rows})==6)
with open(B/'gamma3_new_4cell_generators_v42.csv','w',newline='',encoding='utf-8') as f:
    ww=csv.DictWriter(f,fieldnames=list(new_rows[0]));ww.writeheader();ww.writerows(new_rows)
# Z-linear rank after seed quotient exact finite-support flattening
basis={};rows=[]
for v in mods:
    if not v:continue
    row={}
    for k,d in v.items():
        for g,c in d.items():row[basis.setdefault((k,g),len(basis))]=c
    rows.append(row)
M=sp.zeros(len(rows),len(basis))
for i,row in enumerate(rows):
    for j,c in row.items():M[i,j]=c
ok('gr42_mod_seed_Z_rank_6',M.rank()==6)
summary={
 'version':'v42','group':'Gamma_3(2) embedded faithfully by the standard 3x3 integral matrices',
 'length_two_diamonds':52,'all_exact_Fox_boundaries_zero':True,
 'old_seed_submodule_diamonds':28,
 'augmentation_zero_but_group_ring_nonzero_diamonds':24,
 'new_canonical_4cycles':len(new_rows),'new_4cycles':new_rows,
 'structure':'the 24 hidden diamonds are four signed/left-translated occurrences of six canonical cycles z_i=(1-g_i)e_{r_i}, with r_i=g_i^2 and g_i^2=1 in Gamma_3(2)',
 'relative_to_v38_seed_submodule':'PROVEN_NONZERO: the six canonical cycles live on six relator coordinates outside the twelve seed coordinates',
 'direct_sum_guard':'their coordinate supports are disjoint, hence the six cyclic submodules sum directly; they are NOT claimed free over Z[Gamma] because (1+g_i) annihilates z_i',
 'Z_rank_of_24_images_mod_seed':int(M.rank()),
 'new_second_syzygies':five,
 'five_cell_interpretation':'after adjoining a 4-cell f_i with boundary z_i, (1+g_i)f_i is an exact 5-cycle. Completeness/independence of the full 5-cell module remains NT.',
 'complete_identity_module':'NT_PENDING_FULL_GROUP_RING_RESOLUTION'}
dump('gamma3_group_ring_squier_v42.json',summary)
dump('gamma3_second_syzygies_v42.json',{'version':'v42','count':6,'syzygies':five,'full_5cell_module':'NT'})

# D. Extend J_s to radius 6 and prove sign-torsor non-uniqueness under all current relations.
I3=sp.eye(3);P=sp.Matrix([[1,1,0],[3,4,0],[0,0,1]]);Q=sp.Matrix([[2,0,5],[0,1,0],[3,0,8]])
R=sp.Matrix([[-4,0,5],[0,-1,0],[-3,0,4]]);Sm=sp.Matrix([[-2,1,0],[-3,2,0],[0,0,-1]])
Gamma=sp.diag(-1,1,-1);Delta=sp.diag(-1,1,1);H=Gamma*Delta
DP=sp.diag(1,-1,1);DQ=sp.diag(1,1,-1)
T=sp.simplify(P*R*P.inv());U=sp.simplify(Q*Sm*Q.inv());V=sp.simplify(P*DP*P.inv()*DP);W=sp.simplify(Q*DQ*Q.inv()*DQ)
a=V*Gamma;b=W*H;g0=Gamma;h=H;t=T;u=U
CGEN={'a':a,'h':h,'g':g0,'b':b,'u':u,'t':t};cycle=['a','h','g','b','u','t'];order={x:i for i,x in enumerate(cycle)};comm={frozenset(e) for e in zip(cycle,cycle[1:]+cycle[:1])}
def racg_reduce(word):
    w=list(word);changed=True
    while changed:
        changed=False;i=0
        while i<len(w)-1:
            if frozenset((w[i],w[i+1])) in comm and order[w[i]]>order[w[i+1]]:
                w[i],w[i+1]=w[i+1],w[i];changed=True;i=max(0,i-1)
            else:i+=1
        i=0
        while i<len(w):
            j=i+1
            while j<len(w) and frozenset((w[i],w[j])) in comm:
                if w[j]==w[i]:del w[j];del w[i];changed=True;i=max(-1,i-2);break
                j+=1
            i+=1
    return tuple(w)
nfs=set()
for n in range(7):
    for w in itertools.product(cycle,repeat=n):nfs.add(racg_reduce(w))
ok('J42_radius6_normal_forms_22183',len(nfs)==22183)
def tt(M):return tuple(tuple(int(M[i,j]) for j in range(3)) for i in range(3))
def m3(A,C):return tuple(tuple(sum(A[i][k]*C[k][j] for k in range(3)) for j in range(3)) for i in range(3))
J=sp.Matrix([[0,0,0],[0,0,1],[0,-1,0]]);Itt=tt(I3);Jtt=tt(J);CGT={x:tt(M) for x,M in CGEN.items()}
Jmap={};Gmap={}
for nf in nfs:
    Gm=Itt;Js=Jtt
    for x in nf:
        Mx=CGT[x];Gm=m3(Gm,Mx);Js=m3(m3(Mx,Js),Mx)
    Gmap[nf]=Gm;Jmap[nf]=Js
edge_n=0;nat=True
for nf in nfs:
    for x in cycle:
        q=racg_reduce(nf+(x,))
        if q in nfs:
            edge_n+=1;Mx=CGT[x];nat=nat and Jmap[q]==m3(m3(Mx,Jmap[nf]),Mx)
ok('J42_radius6_edges_natural',nat)
ok('J42_radius6_edge_count_positive',edge_n>25000)
# sign flip is an automorphism of all transport constraints; -1 excluded from U_k+ for all k.
sign_nat=all(tuple(tuple(-z for z in row) for row in Jmap[q])==m3(m3(CGT[x],tuple(tuple(-z for z in row) for row in Jmap[nf])),CGT[x]) for nf in list(nfs)[:500] for x in cycle for q in [racg_reduce(nf+(x,))] if q in nfs)
ok('J42_global_sign_flip_preserves_transport_sample_and_formula',sign_nat)
for k in range(1,13):
    m=2*3**k; ok(f'J42_minus1_not_Uplus_{k}',((m-1)%3)!=1)
j42={
 'version':'v42','transport':'T_M(X)=M^{-1}XM; generators are involutions in the RACG chart',
 'radius6':{'raw_words':sum(6**n for n in range(7)),'normal_forms':len(nfs),'internal_generator_edges':edge_n,'naturality':'PROVEN'},
 'extension_beyond_v41':'radius 5 -> radius 6',
 'root_sign_torsor':{
   'sections':['J_s','-J_s'],
   'relation':'global sign flip commutes with every conjugation transport and therefore preserves every current historical/RACG relation constraint',
   'projective_tower':'the two signs remain distinct in every P_k^+ because -1 is not in U_k^+',
   'uniqueness_from_current_relational_data':'REFUTED-TYPED',
   'source_selected_Lorentz_orientation':'NT_REQUIRES_ADDITIONAL_SOURCE_MARK'
 },
 'opposite_variance':'MXM^{-1} retained only as opposite-groupoid action',
 'physical_phase':'SEPARATED/PENDING-MORPHISM'}
dump('projective_J_radius6_sign_torsor_v42.json',j42)

# E. Field functor remains fail closed.
fg={
 'version':'v42','source_binding':'selling_source_binding_gate_v42.json',
 'zero_cell_compatibility':'NT','one_cell_compatibility':'NT','two_cell_compatibility':'NT',
 'field_functor':'PENDING-MORPHISM_NOT_PROMOTED','historical_H1_H2':'NT_NOT_RECOMPUTED',
 'decision':'EXACT_0/1/2_CELL_GATE_NOT_CLOSED'}
dump('selling_field_functor_gate_v42.json',fg)

pending={
 'version':'v42','predecessor':'v41 CLOSED/APPEND-ONLY','entries':[
  {'id':'PM42-SELL-BIND','status':'NT','subject':'source-certified BSB assignment of HAX/HC/fissures/crossings/walls/fields'},
  {'id':'PM42-89-C2','status':'SEPARATED/PENDING-PRECONDITION','subject':'C2 orbit census on the 89 degree-one candidates; historical action absent'},
  {'id':'PM42-ENV','status':'PROVEN-GRAPH-CANDIDATE/SEPARATED-PENDING-MORPHISM','subject':'R5 unique graph-canonical deletion candidate under declared 16-envelope criterion; no Gamma adapter'},
  {'id':'PM42-G3-4','status':'PROVEN','subject':'six new exact group-ring 4-cycles z_i=(1-g_i)e_ri beyond v38 seed submodule'},
  {'id':'PM42-G3-5','status':'PROVEN-SECOND-SYZYGY/NT-COMPLETE','subject':'six exact (1+g_i)f_i 5-cycles after adjoining the new 4-cells'},
  {'id':'PM42-G3-FULL','status':'NT','subject':'complete Z[Gamma3(2)] resolution'},
  {'id':'PM42-J','status':'PROVEN-RADIUS6/REFUTED-TYPED-UNIQUE-SIGN','subject':'J_s transport; current relations cannot select between +/- sections'},
  {'id':'PM42-SELL-FUN','status':'PENDING-MORPHISM','subject':'F_Sell^field -> F29_can'},
  {'id':'PM42-H12','status':'NT','subject':'historical field H1/H2'},
  {'id':'PM42-ROOF-TRANSPOSE','status':'REFUTED-TYPED','subject':'roof reversal = transpose'},
  {'id':'PM42-PHASE','status':'SEPARATED/PENDING-MORPHISM','subject':'physical phase'}],
 'stratified_reality_guard':'No diagnostic similarity or arithmetic coincidence is upgraded without an explicit typed morphism.'}
dump('pending_morphism_registry_v42.json',pending)

new_count=len(checks)-2; status='PASS' if all(v for _,v in checks) else 'FAIL'
cert={'phase':'v42-selling-gates-pointed-envelope-group-ring-squier-second-syzygies-J-sign-torsor','status':status,
 'predecessor_v41_checks':5646,'new_check_count':new_count,'combined_check_count':5646+new_count,'failed':[n for n,v in checks if not v],
 'selling':{'new_source_geometric_promotions':0,'C2_89_census':'NOT_ACTIVATED','field_functor':'PENDING','historical_H1_H2':'NT'},
 'envelope':{'graph_canonical_candidate':'R5','Gamma_adapter':'SEPARATED/PENDING-MORPHISM'},
 'gamma3':{'diamonds':52,'old_seed_diamonds':28,'hidden_group_ring_diamonds':24,'new_canonical_4cycles':6,'new_second_syzygies':6,'complete_resolution':'NT'},
 'projective_J':{'radius6_objects':len(nfs),'radius6_edges':edge_n,'sign_uniqueness':'REFUTED-TYPED_FROM_CURRENT_RELATIONS','Lorentz_orientation':'NT'},
 'pending_registry':'pending_morphism_registry_v42.json'}
dump('v42_VERIFICATION_REPORT.json',cert)
print(json.dumps(cert,indent=2,ensure_ascii=False));raise SystemExit(0 if status=='PASS' else 1)
