#!/usr/bin/env python3
import json,csv,itertools,hashlib
from pathlib import Path
import sympy as sp
ROOT=Path('/mnt/data/v54_gate2/r2pure')
audit=json.load(open('/mnt/data/v54_research/gamma3_r2pure_group_ring_audit_v52.json'))
target=json.load(open(ROOT/'multiobject_lift_groupoid_v34.json'))
# Source exact overlap 1-skeleton
verts=sorted(audit['fingerprints'])
edges=[]
for r in audit['pairwise_exact_overlap']:
    if r['exact_shared_coordinate_count']>0:
        edges.append((r['a'],r['b'],r['exact_shared_coordinate_count'],r['signed_dot']))
assert len(verts)==13 and len(edges)==54
# Target K component and exact v29 bigon.
E={e['label']:e for e in target['edges'] if e['source'].startswith('K') or e['target'].startswith('K')}
MM={k:sp.Matrix(e['matrix']) for k,e in E.items()}
chain_labels=['v29_S1','v29_S2','v29_S1','v29_S2']
chain=sp.eye(3)
for lab in chain_labels: chain=chain*MM[lab]
deck=MM['kappa_12_deck']
assert chain==deck
# Normalize an explicit non-injective quotient by source-native R2PURE index modulo five.
# This is deliberately declared normalization-dependent; it is an existence witness only.
objmap={v:f"K{(int(v.split('_')[1])-1)%5}" for v in verts}
classes={f'K{i}':[] for i in range(5)}
for v,k in objmap.items():classes[k].append(v)
assert sorted(map(len,classes.values()))==[2,2,3,3,3]
# Potentials from K0. K4 uses the direct deck arrow, so the exact bigon is actively available.
pot={
 'K0':[],
 'K1':['v29_S1'],
 'K2':['v29_S1','v29_S2'],
 'K3':['v29_S1','v29_S2','v29_S1'],
 'K4':['kappa_12_deck'],
}
def invtok(t): return t[:-3] if t.endswith('^-1') else t+'^-1'
def reduce_word(word):
    st=[]
    for t in word:
        if st and invtok(t)==st[-1]: st.pop()
        else: st.append(t)
    return st
def path_word(a,b): return reduce_word([invtok(t) for t in reversed(pot[a])] + pot[b])
def mat_tok(t):
    return MM[t[:-3]].inv() if t.endswith('^-1') else MM[t]
def path_matrix(word):
    A=sp.eye(3)
    for t in word:A=A*mat_tok(t)
    return sp.simplify(A)
# verify potentials hit stored target charts and every image is a legitimate congruence.
Q={k:sp.Matrix(v['Q']) for k,v in target['objects'].items() if k.startswith('K')}
charts={k:sp.Matrix(v['chart']) for k,v in target['objects'].items() if k.startswith('K')}
Pmat={k:path_matrix(v) for k,v in pot.items()}
for k in pot:
    assert charts['K0']*Pmat[k]==charts[k]
images=[]; lens={}; identity=0; deck_used=0
for a,b,shared,dot in edges:
    ka,kb=objmap[a],objmap[b]; w=path_word(ka,kb); M=path_matrix(w)
    assert charts[ka]*M==charts[kb]
    assert sp.simplify(M.T*Q[ka]*M-Q[kb])==sp.zeros(3)
    if not w:identity+=1
    if any('kappa_12_deck' in t for t in w):deck_used+=1
    lens[len(w)]=lens.get(len(w),0)+1
    images.append({'source':[a,b],'source_shared':shared,'source_signed_dot':dot,'target_objects':[ka,kb],'target_path':w,'target_path_length':len(w)})
# Diagnostic only: every clique triangle of overlap graph has potential-telescoping identity.
adj={v:set() for v in verts}
for a,b,_,_ in edges:adj[a].add(b);adj[b].add(a)
tris=[]
for t in itertools.combinations(verts,3):
    if t[1] in adj[t[0]] and t[2] in adj[t[0]] and t[2] in adj[t[1]]: tris.append(t)
assert len(tris)==107
for a,b,c in tris:
    M=path_matrix(path_word(objmap[a],objmap[b]))*path_matrix(path_word(objmap[b],objmap[c]))*path_matrix(path_word(objmap[c],objmap[a]))
    assert M==sp.eye(3)
# Whiskering-forgetful mod-2 shadow: augment every group word to 1 and reduce coefficients mod 2.
rows=list(csv.DictReader(open(ROOT/'gamma3_radius2_pure_lifts_v47.csv')))
agg={}
for r in rows:
    key=(r['lift_id'],r['relator']);agg[key]=agg.get(key,0)+int(r['coefficient'])
odd=[{'lift':k[0],'relator':k[1],'augmentation':v} for k,v in agg.items() if v%2]
assert not odd
# There is no source-declared R2PURE 2-cell presentation in the frozen inputs.
cert={
 'version':'v54',
 'status':'PASS-EXPLICIT-NONINJECTIVE-EDGE-TO-PATH-1-SKELETON / 2CELL-WHISKERING-GATE-PENDING',
 'source':{'objects':13,'overlap_edges':54,'definition':'exact shared group-ring coordinate > 0','declared_2cell_presentation':'NOT-MATERIALIZED-IN-FROZEN-R2PURE-INPUTS'},
 'target':{'component':['K0','K1','K2','K3','K4'],'edges':5,'v29_bigon':{'chain':chain_labels,'direct':'kappa_12_deck','matrix_equality':True,'status':'PROVEN-EXACT-TARGET-2CELL'}},
 'quotient':{
   'normalization':'R2PURE source index modulo 5; existence witness, not claimed canonical',
   'object_map':objmap,'fibres':classes,'noninjective':True,
   'edge_rule':'F(a->b)=p(Fa)^-1 p(Fb) with p(K4)=direct kappa_12_deck; all paths lie in stored K component',
   'edge_images':images,'path_length_histogram':lens,'identity_edge_images':identity,'edge_images_using_deck_arrow':deck_used,
   'all_54_target_congruences':'PASS'
 },
 'diagnostic_clique_triangles':{
   'count':len(tris),'all_images_matrix_identity':True,
   'status':'PASS-DIAGNOSTIC-ONLY/NOT-PROMOTED-AS-SOURCE-2CELLS'
 },
 'whiskering_shadow':{
   'operation':'group words -> identity, coefficients -> mod 2, retain relator label',
   'rows':len(rows),'odd_relator_augmentation_residuals':len(odd),
   'status':'ZERO-ALL-13-LIFTS',
   'interpretation':'the simplest exact whiskering-forgetful C2 shadow annihilates all thirteen lifts; it supplies no nonzero Selling/deck invariant.'
 },
 'decision':{
   'explicit_noninjective_edge_to_path_functor_on_1_skeleton':'PROVEN-EXISTS/NORMALIZATION-DEPENDENT',
   'target_v29_bigon_compatibility':'PROVEN-AVAILABLE',
   'source_2cell_compatibility':'NT/BLOCKED-BY-MISSING-SOURCE-R2PURE-2CELL-PRESENTATION',
   'source_generator_whiskering_equivariance':'NT; prior literal elementary whiskering permutation action on the 13 objects is REFUTED-TYPED',
   'R2PURE13_to_SellingDeck17_adapter':'SEPARATED/PENDING-MORPHISM'
 },
 'guard':'This construction proves existence of a noninjective path-valued quotient of the exact overlap 1-skeleton. It is not promoted to the requested incidence/groupoid adapter because the frozen R2PURE carrier does not supply a 2-cell presentation and the source group-ring whiskering data is not transported nontrivially.'
}
out=ROOT/'r2pure13_to_sellingdeck17_noninjective_path_quotient_v54.json'
out.write_text(json.dumps(cert,indent=2,sort_keys=True)+"\n")
print(json.dumps({'status':cert['status'],'fibres':{k:len(v) for k,v in classes.items()},'path_lengths':lens,'identity':identity,'deck_used':deck_used,'clique_triangles_diag':len(tris),'mod2_odd':len(odd)},indent=2))
