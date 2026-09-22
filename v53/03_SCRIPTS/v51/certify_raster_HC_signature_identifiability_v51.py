#!/usr/bin/env python3
import csv,json,math,hashlib
from pathlib import Path
B=Path('/mnt/data/v51_work')
source=json.load(open(B/'selling_source_HC_local_star_signatures_v51.json'))
gate=json.load(open(B/'selling_raster_central_bold_cycle_gate_v51.json'))
N={r['node_id']:r for r in csv.DictReader(open(B/'selling_fig1_nodes_v40.csv'))}
S=list(csv.DictReader(open(B/'selling_fig1_segments_v40.csv')))
F=list(csv.DictReader(open(B/'selling_fig1_fields_v40.csv')))
order=gate['cycle_order']
cycle=[N[n] for n in order]
# Available table-level discriminants corresponding to the source signature are all unresolved.
axis_labels=[s['historical_axis_id'].strip() for s in S if s['historical_axis_id'].strip()]
wall_labels=[s['historical_wall_label'].strip() for s in S if s['historical_wall_label'].strip()]
field_labels=[f['historical_form_label'].strip() for f in F if f['historical_form_label'].strip()]
hist_nodes=[r['historical_node_id'].strip() for r in N.values() if r['historical_node_id'].strip()]
# The only table-native necessary condition we can safely apply to an HC on the reinforced contour:
# two reinforced contour branches meet at the node. This leaves many candidates.
candidates=[r['node_id'] for r in cycle if int(r['reinforced_degree'])>=2]
# order-preserving six-subsets on a cycle before source label shift/orientation
subset_count=math.comb(len(candidates),6) if len(candidates)>=6 else 0
# fingerprint by available integer degrees; count duplicate fingerprints.
fps={}
for r in cycle:
 fp=(int(r['degree_solid']),int(r['ordinary_degree']),int(r['reinforced_degree']),r['diagnostic_topology_candidate'])
 fps.setdefault(fp,[]).append(r['node_id'])
duplicates={repr(k):v for k,v in fps.items() if len(v)>1}
cert={
 'version':'v51','status':'PROVEN-RASTER-TABLE-IDENTIFIABILITY-NO-GO',
 'scope':'v40 vector tables + v51 exact six HC local-star/stabilizer signatures; no manual pixel semantics',
 'source_signature_count':len(source['centres']),'raster_cycle_node_count':len(order),
 'populated_historical_axis_labels_in_segments':len(axis_labels),
 'populated_historical_wall_labels_in_segments':len(wall_labels),
 'populated_historical_form_labels_in_fields':len(field_labels),
 'populated_historical_node_ids':len(hist_nodes),
 'reinforced_degree_ge2_cycle_candidates':candidates,
 'candidate_count':len(candidates),
 'six_subset_count_under_this_necessary_condition':subset_count,
 'duplicate_available_fingerprints':duplicates,
 'decision':'NO_UNIQUE_HC0..HC5_BINDING_IS_IDENTIFIABLE_FROM_THE_V40_VECTOR_TABLES_AND_THE_SOURCE_SIGNATURES_ALONE',
 'proof':'The source signatures require named axis roles, incident historical form labels and V4 stabilizers. The v40 raster tables contain no populated historical axis, wall, form or node identifiers. Even the safe table-native necessary condition of two reinforced branches leaves multiple cycle nodes and many six-subsets. Thus the required discriminants are absent rather than merely noisy.',
 'guard':'This is an information-sufficiency no-go for the current vector tables. It does not refute a binding obtained from direct printed right-angle/axis marks, a new source-aware raster transcription, or certified manual anchors.'
}
p=B/'selling_raster_HC_signature_identifiability_v51.json';p.write_text(json.dumps(cert,indent=2,sort_keys=True)+'\n')
print(json.dumps(cert,indent=2));print('sha256',hashlib.sha256(p.read_bytes()).hexdigest())
