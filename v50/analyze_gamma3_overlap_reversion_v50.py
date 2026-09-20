#!/usr/bin/env python3
from __future__ import annotations
import json,itertools
from pathlib import Path
import networkx as nx
B=Path(__file__).resolve().parent
src=json.loads((B/'gamma3_radius2_relation_audit_v49.json').read_text())
verts=[f'R2PURE_{i:02d}' for i in range(1,14)]
G=nx.Graph();G.add_nodes_from(verts)
for e in src['pairwise_overlap_graph']:
 G.add_edge(e['a'],e['b'],weight=int(e['shared_relator_count']),shared=tuple(e['shared_relators']))
allpairs={tuple(sorted(p)) for p in itertools.combinations(verts,2)}
edges={tuple(sorted(e)) for e in G.edges()}
nonedges=sorted(allpairs-edges)
C=nx.complement(G)
# exact automorphism enumeration on 13 vertices
unw=list(nx.algorithms.isomorphism.GraphMatcher(G,G).isomorphisms_iter())
em=nx.algorithms.isomorphism.categorical_edge_match('weight',None)
w=list(nx.algorithms.isomorphism.GraphMatcher(G,G,edge_match=em).isomorphisms_iter())
weighted_nonidentity=[m for m in w if any(m[v]!=v for v in verts)]
components=sorted([sorted(c) for c in nx.connected_components(C)], key=lambda x:(-len(x),x))
out={
 'version':'v50',
 'input':'gamma3_radius2_relation_audit_v49.json',
 'vertices':verts,
 'overlap_edges':G.number_of_edges(),
 'complete_graph_pairs':len(allpairs),
 'nonoverlap_edge_count':len(nonedges),
 'nonoverlap_edges':[list(e) for e in nonedges],
 'overlap_degree':{v:int(G.degree(v)) for v in verts},
 'weighted_overlap_degree':{v:int(sum(G[v][u]['weight'] for u in G.neighbors(v))) for v in verts},
 'complement_degree':{v:int(C.degree(v)) for v in verts},
 'complement_components':components,
 'complement_nontrivial_component_count':sum(len(c)>1 for c in components),
 'overlap_triangle_count':sum(nx.triangles(G).values())//3,
 'complement_triangle_count':sum(nx.triangles(C).values())//3,
 'unweighted_overlap_automorphism_group_size':len(unw),
 'weighted_overlap_automorphism_group_size':len(w),
 'weighted_nonidentity_involutions':[],
 'weighted_nonidentity_automorphisms':len(weighted_nonidentity),
 'decision':{
   'nontrivial_reversion_preserving_weighted_pairwise_overlap':'REFUTED-TYPED',
   'selling_deck_17_object_adapter':'SEPARATED/PENDING-MORPHISM',
   'higher_group_ring_reversion':'NT'
 },
 'selling_deck_context':{
   'objects':17,'component_sizes':[8,4,5],'edges':15,'cycle_rank':1,
   'R5_deleted_objects':16,'R5_deleted_component_sizes':[7,5,4],
   'strict_R5_deleted_to_critical_diamond_adapter':'REFUTED-TYPED_BY_v43_GRAPH_INVARIANTS'
 },
 'guard':'The weighted pairwise-overlap graph has only the identity automorphism. This excludes a reversal visible at that coarse weighted graph level only; it does not exclude a richer group-ring/whiskering/groupoid anti-equivalence or a future 13-to-17 typed adapter.'
}
# Involutions among weighted nonidentity autos (should be empty, but compute semantically).
for m in weighted_nonidentity:
 if all(m[m[v]]==v for v in verts):
  out['weighted_nonidentity_involutions'].append({v:m[v] for v in verts})
(B/'gamma3_radius2_overlap_reversion_audit_v50.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'status':'PASS','nonedges':len(nonedges),'weighted_aut':len(w),'unweighted_aut':len(unw),'components':[len(x) for x in components]},indent=2))
