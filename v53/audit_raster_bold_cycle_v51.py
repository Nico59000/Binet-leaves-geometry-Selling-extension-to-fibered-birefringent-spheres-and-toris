#!/usr/bin/env python3
import csv,json,math,hashlib,itertools
from pathlib import Path
import networkx as nx
B=Path('/mnt/data/v51_work');V50=Path('/mnt/data/v50_work')
nodes={}
with open(V50/'selling_fig1_nodes_v40.csv') as f:
 for r in csv.DictReader(f):nodes[r['node_id']]=(float(r['x_px300']),float(r['y_px300']))
G=nx.Graph(); seg={}
with open(V50/'selling_fig1_segments_v40.csv') as f:
 for r in csv.DictReader(f):
  if r['layer']!='reinforced':continue
  x1,y1=nodes[r['node_u']];x2,y2=nodes[r['node_v']]
  # central bold carrier bounding box
  if min(x1,x2)>=630 and max(x1,x2)<=1300 and min(y1,y2)>=1250 and max(y1,y2)<=2070:
   G.add_edge(r['node_u'],r['node_v'],id=r['segment_id'],length=float(r['length_px300']));seg[r['segment_id']]=r
comps=sorted(nx.connected_components(G),key=lambda c:-len(c))
assert [len(c) for c in comps]==[11,5], [len(c) for c in comps]
ends=[]
for ci,c in enumerate(comps):
 S=G.subgraph(c); ee=[v for v,d in S.degree() if d==1]; assert len(ee)==2
 for v in ee:ends.append((ci,v))
# minimal cross-component perfect matching between endpoints: two short gaps
a=[v for ci,v in ends if ci==0];b=[v for ci,v in ends if ci==1]
def dist(u,v):
 x,y=nodes[u];X,Y=nodes[v];return math.hypot(x-X,y-Y)
matchings=[[(a[0],b[0]),(a[1],b[1])],[(a[0],b[1]),(a[1],b[0])]]
ms=sorted([(sum(dist(u,v) for u,v in m),m) for m in matchings],key=lambda z:z[0])
best=ms[0][1];other=ms[1][0]
for u,v in best:G.add_edge(u,v,id='BRIDGE',length=dist(u,v),bridge=True)
assert nx.is_connected(G) and all(d==2 for _,d in G.degree()) and G.number_of_edges()==G.number_of_nodes()==16
cycle=nx.cycle_basis(G);assert len(cycle)==1 and len(cycle[0])==16
# canonical ordered cycle starting at lexicographically smallest node, both orientations choose lexicographically smaller sequence
cyc=cycle[0]
def canon(seq):
 mi=min(range(len(seq)),key=lambda i:seq[i]);seq=seq[mi:]+seq[:mi];return seq
c1=canon(cyc);c2=canon([cyc[0]]+list(reversed(cyc[1:])))
order=min(c1,c2)
# list edges in order
edges=[]
for i,u in enumerate(order):
 v=order[(i+1)%len(order)];d=G[u][v]
 edges.append({'u':u,'v':v,'id':d['id'],'length_px300':d['length'],'bridge':bool(d.get('bridge',False))})
out={
 'version':'v51','status':'PROVEN-RASTER-CENTRAL-REINFORCED-16-CYCLE; HC/HAX-SPLIT-NT',
 'source_tables':['selling_fig1_nodes_v40.csv','selling_fig1_segments_v40.csv'],
 'central_bbox_px300':[630,1250,1300,2070],
 'reinforced_components_before_gap_bridge':[11,5],
 'bridge_pairs':[{'nodes':[u,v],'distance_px300':dist(u,v)} for u,v in best],
 'alternative_matching_total_gap_px300':other,
 'cycle_node_count':16,'cycle_edge_count':16,'cycle_order':order,'cycle_edges':edges,
 'source_axis_centre_required_count':6,
 'ordered_six_cut_subsets_without_extra_semantics':math.comb(16,6),
 'decision':{
   'central_reinforced_raster_cycle':'PROVEN-DIAGNOSTIC-GEOMETRIC-CARRIER',
   'unique_binding_of_six_HC_to_cycle_nodes':'NT/NONUNIQUE-FROM-V40-VECTOR-TABLES',
   'source_to_raster_62_face_realization':'PENDING-DIRECT-HC-MARK/AXIS-BRANCH-ANCHORS'
 },
 'guard':'The two short bridges only repair segmentation gaps in the already visible reinforced stroke. They do not identify historical HC_i. Choosing six of the 16 cycle nodes without direct right-angle/axis evidence would be arbitrary.'
}
p=B/'selling_raster_central_bold_cycle_gate_v51.json';p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':out['status'],'bridges':out['bridge_pairs'],'order':order,'six_cut_subsets':out['ordered_six_cut_subsets_without_extra_semantics'],'sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
