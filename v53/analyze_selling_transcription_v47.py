#!/usr/bin/env python3
from pathlib import Path
import json,csv,itertools,collections
B=Path(__file__).resolve().parent
v45=json.loads((B/'selling_bsb_axis_binding_attempt_v45.json').read_text(encoding='utf-8'))
traces=v45['trace_ids']
edges=[]
for e in v45['segment_extent_intersections']:
    edges.append((e['a'],e['b'],e['point_px300']))
adj={t:set() for t in traces}
for a,b,p in edges:
    adj[a].add(b);adj[b].add(a)
# connected components
seen=set(); comps=[]
for t in traces:
    if t in seen: continue
    q=[t];seen.add(t);c=[]
    while q:
        u=q.pop();c.append(u)
        for v in sorted(adj[u]):
            if v not in seen:seen.add(v);q.append(v)
    comps.append(sorted(c))
# enumerate undirected simple cycles (canonicalized), length>=3
def canon_cycle(c):
    # c list without repeated endpoint
    rots=[]
    n=len(c)
    for s in (c,list(reversed(c))):
        for i in range(n):rots.append(tuple(s[i:]+s[:i]))
    return min(rots)
cycles=set()
for start in traces:
    stack=[(start,[start])]
    while stack:
        u,path=stack.pop()
        for v in adj[u]:
            if v==start and len(path)>=3:
                cycles.add(canon_cycle(path))
            elif v not in path and len(path)<7:
                stack.append((v,path+[v]))
cycles=sorted(cycles,key=lambda c:(len(c),c))
cycle_lengths=collections.Counter(map(len,cycles))
# Explicit abstract source carrier: 6-cycle of historical boundary axes / fourfold centres.
abstract=[{'axis':f'HAX{i}','boundary':[f'HC{i}',f'HC{(i+1)%6}']} for i in range(6)]
# diagnostic incidence table
with (B/'selling_axis_stitching_incidence_v47.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=['trace_a','trace_b','x_px300','y_px300','status']);w.writeheader()
    for a,b,p in edges:
        w.writerow({'trace_a':a,'trace_b':b,'x_px300':p[0],'y_px300':p[1],
                    'status':'DIAGNOSTIC_DAX_EXTENT_INTERSECTION_NOT_HISTORICAL_HC'})
trans={
 'version':'v47',
 'primary_source':'E. Selling, Des formes quadratiques binaires et ternaires, JMPA 3e serie 3 (1877), pp.153-206; Fig.1 / Tafel III',
 'source_rules':[
  {'printed_page':'175','statement':'The six external boundary lines of Fig.1 are bipartite symmetry axes and their six intersections are fourfold symmetry centres.','status':'PROVEN_SOURCE_TEXT'},
  {'printed_pages':'182-184','statement':'Dotted lines denote symmetry axes; the six indicated axis intersections are fourfold centres; the drawing may be bent when imagining the corresponding perpendicular axes as straight.','status':'PROVEN_SOURCE_DRAWING_AND_TOPOLOGICAL_RULE'},
  {'printed_pages':'189-191','statement':'External-boundary routes P,Q,R,S generate the repeated-boundary transformation grammar; P is attached to the lower-right boundary and Q to the extreme-right boundary; R,S use doubled leftward routes across symmetry axes.','status':'PROVEN_SOURCE_ROUTE_SEMANTICS'}],
 'historical_abstract_axis_centre_carrier':{
   'vertices':[f'HC{i}' for i in range(6)],'edges':abstract,'topology':'C6','status':'PROVEN_SOURCE_SEMANTIC_INCIDENCE'},
 'plate_transcription_model':'HAX_i are source boundary/symmetry curves. DAX_j are local straight diagnostic traces only; Euclidean straightness/angle of a DAX is not a source invariant.',
 'diagnostic_DAX_graph':{
   'vertices':traces,'edges':[{'a':a,'b':b,'point_px300':p} for a,b,p in edges],
   'degrees':{t:len(adj[t]) for t in traces},'connected_components':comps,
   'simple_cycle_lengths':dict(sorted(cycle_lengths.items())),
   'simple_cycles':[list(c) for c in cycles],
   'has_six_cycle':any(len(c)==6 for c in cycles),
   'isolated_traces':[t for t in traces if not adj[t]]},
 'raster_HC_coordinates':None,
 'raster_HAX_curves':None,
 'boundary_letter_attachment':None,
 'PQRS_raster_routes':None,
 'decision':'SOURCE ABSTRACT C6 CARRIER RETAINED, BUT CURRENT SEVEN DAX FRAGMENTS DO NOT CONTAIN A SIX-CYCLE OF PAIRWISE EXTENT INTERSECTIONS AND DAX00 IS ISOLATED; THEREFORE DAX-ONLY STITCHING CANNOT SOURCE-CERTIFY THE SIX HC/HAX. DIRECT PLATE TRANSCRIPTION OF THE SIX INDICATED CENTRES AND THEIR LOCAL DOTTED BRANCHES IS STILL REQUIRED.',
 'status':'NT/PENDING-SOURCE-PLATE-TRANSCRIPTION',
 'guard':'Absence of a C6 in the diagnostic DAX graph is a no-go only for stitching the currently detected straight trace extents. It does not refute the historical curved axes visible in the source.'}
(B/'selling_axis_curve_transcription_v47.json').write_text(json.dumps(trans,indent=2,ensure_ascii=False),encoding='utf-8')
routes={
 'version':'v47','status':'PROVEN_SOURCE_ROUTE_SEMANTICS / RASTER_BINDING_NT',
 'routes':{
  'P':{'source_description':'lower-right external boundary route','source_status':'PROVEN','raster_polyline':None},
  'Q':{'source_description':'extreme-right external boundary route, traversed from lower to upper endpoint','source_status':'PROVEN','raster_polyline':None},
  'R':{'source_description':'leftward external boundary route interpreted after crossing/doubling across a symmetry axis','source_status':'PROVEN','raster_polyline':None},
  'S':{'source_description':'analogous doubled leftward external boundary route','source_status':'PROVEN','raster_polyline':None}},
 'relations':{'T':'P R P^-1','U':'Q S Q^-1','TU_equals_UT':'PROVEN_SOURCE_RELATION'},
 'decision':'SEMANTIC ROUTE ANCHORS ARE SOURCE-CERTIFIED, BUT NO P/Q/R/S LABEL IS YET ATTACHED TO A MACHINE-READABLE FIG.1 ARC.'}
(B/'selling_PQRS_route_semantics_v47.json').write_text(json.dumps(routes,indent=2,ensure_ascii=False),encoding='utf-8')
c2={
 'version':'v47','degree_one_candidates':89,
 'precondition':'six source-certified raster HAX curves and reflection action preserving the 136-node/105-edge carrier',
 'precondition_status':'FALSE_HC_HAX_RASTER_TRANSCRIPTION_NOT_TOTAL',
 'activated':False,'fixed_points':None,'free_orbits':None,'total_orbits':None,
 'comparison_44_45':'SEPARATED_NOT_TESTED','comparison_34_55':'SEPARATED_NOT_TESTED',
 'decision':'C2_ORBIT_CENSUS_NOT_ACTIVATED_FAIL_CLOSED'}
(B/'selling_degree1_C2_gate_v47.json').write_text(json.dumps(c2,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'trace_count':len(traces),'edge_count':len(edges),'components':comps,'cycle_lengths':dict(cycle_lengths),'has_six_cycle':trans['diagnostic_DAX_graph']['has_six_cycle'],'isolated':trans['diagnostic_DAX_graph']['isolated_traces'],'status':trans['status']},indent=2))
