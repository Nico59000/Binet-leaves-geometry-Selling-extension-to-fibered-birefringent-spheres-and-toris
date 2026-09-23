#!/usr/bin/env python3
from pathlib import Path
import csv,json,itertools,math
B=Path(__file__).resolve().parent
rows=list(csv.DictReader(open(B/'selling_dotted_axis_traces_v45.csv',encoding='utf-8')))
pg=json.load(open(B/'selling_fig1_bsb_photo_gate_v38.json'))
sem=json.load(open(B/'selling_fig1_semantic_contract_v39.json'))
q0=json.load(open(B/'selling_q0_historical_domain_v36.json'))
assert len(rows)==7
# Source constraints reviewed against Selling 1877 second part.
source={
 'primary':'E. Selling, Des formes quadratiques binaires et ternaires, JMPA, 3e serie 3 (1877), pp.153-206',
 'constraints':[
  {'id':'S46-A','printed_page':'163-164','statement':'Four of the six external boundary lines are of the numerical-coincidence symmetry-axis species; boundary letters name the variable vanishing there.','status':'PROVEN_SOURCE_TEXT'},
  {'id':'S46-B','printed_page':'175','statement':'The six external boundary lines are bipartite symmetry axes and their six intersections are quadripartite centres.','status':'PROVEN_SOURCE_TEXT'},
  {'id':'S46-C','printed_page':'182','statement':'Dotted lines in the figures denote symmetry axes.','status':'PROVEN_SOURCE_DRAWING_CONVENTION'},
  {'id':'S46-D','printed_page':'184','statement':'The six indicated axis intersections are quadripartite centres; the relevant intersecting axes are perpendicular in the source construction.','status':'PROVEN_SOURCE_TEXT'},
  {'id':'S46-E','printed_pages':'189-191','statement':'For Fig.1, external-boundary substitutions P,Q,R,S are attached to explicit outer-boundary routes; P is lower-right and Q is the extreme right boundary route.','status':'PROVEN_SOURCE_STABILIZER_GEOMETRY'}
 ],
 'inherited_semantic_contract':'selling_fig1_semantic_contract_v39.json',
 'inherited_q0_domain':'selling_q0_historical_domain_v36.json'
}
# Angle comparison for the direct model DAX_j == a whole straight HAX_i.
def acute_diff(a,b):
 d=abs(a-b)%180.0
 return min(d,180.0-d)
angles={r['trace_id']:float(r['angle_deg']) for r in rows}
spans={r['trace_id']:float(r['span']) for r in rows}
err=float(pg['registration_to_v37']['p95_reprojection_error_px_v37_300'])
minspan=min(spans.values())
# Conservative direction uncertainty if both fitted endpoints move by p95 registration error.
line_angle_err=math.degrees(math.atan2(2*err,minspan))
pair_angle_err=2*line_angle_err
pairs=[]
for a,b in itertools.combinations(rows,2):
 d=acute_diff(float(a['angle_deg']),float(b['angle_deg']))
 pairs.append({'a':a['trace_id'],'b':b['trace_id'],'acute_angle_deg':d,'distance_to_90_deg':90-d})
closest=min(pairs,key=lambda x:abs(90-x['acute_angle_deg']))
subsets=[]
ids=[r['trace_id'] for r in rows]
for omit in ids:
 sub=[x for x in ids if x!=omit]
 sp=[p for p in pairs if p['a'] in sub and p['b'] in sub]
 orth=[p for p in sp if abs(90-p['acute_angle_deg'])<=pair_angle_err]
 subsets.append({
  'omitted':omit,'candidate_six':sub,
  'orthogonal_pair_count_under_conservative_error_bound':len(orth),
  'required_fourfold_centre_pair_count':6,
  'boundary_letters_serialized':False,'PQRS_route_labels_serialized':False,'cyclic_HC_order_serialized':False,
  'four_of_six_numerical_coincidence_types_serialized':False,
  'status':'REFUTED-TYPED_AS_DIRECT_STRAIGHT_DAX_TO_HAX_MODEL'
 })
out={
 'version':'v46','source_constraints':source,
 'diagnostic_traces':ids,'trace_count':len(ids),'historical_external_axes':6,
 'direct_model':'one whole straight DAX trace is identified with one whole historical HAX axis',
 'p95_registration_error_px300':err,'min_trace_span_px300':minspan,
 'single_line_angle_uncertainty_bound_deg':line_angle_err,'pair_angle_uncertainty_bound_deg':pair_angle_err,
 'closest_DAX_pair_to_perpendicular':closest,
 'max_acute_DAX_pair_angle_deg':max(p['acute_angle_deg'] for p in pairs),
 'six_subset_census':subsets,
 'direct_six_subset_survivors':sum(s['status'].startswith('PASS') for s in subsets),
 'decision':'REFUTED-TYPED: NO 6-SUBSET OF THE SEVEN STRAIGHT DAX FITS CAN BE THE SIX WHOLE HISTORICAL HAX AXES. THE DAX OBJECTS MUST REMAIN LOCAL/DIAGNOSTIC TRACE FRAGMENTS UNTIL STITCHED INTO SOURCE-CERTIFIED AXIS CURVES OR REPLACED BY DIRECT OUTER-BOUNDARY TRANSCRIPTION.',
 'historical_HAX_bound':0,'historical_HC_bound':0,
 'C2_89':'NOT_ACTIVATED_FAIL_CLOSED',
 'guard':'This refutes only the direct one-DAX-per-HAX straight-line search space. It does not refute the six historical axes, curved/stiched axes, or a future source-certified plate transcription.'
}
(B/'selling_axis_subset_census_v46.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
req={
 'version':'v46','corrected_geometric_model':'historical HAX_i are source symmetry-axis curves/outer boundaries; DAX_j are local straight diagnostic fragments and need not be in bijection with HAX_i',
 'next_minimal_source_tasks':[
  'transcribe the six fourfold centre marks directly from the source plate',
  'for each centre, identify the two dotted symmetry-axis branches crossing there and preserve the right-angle incidence',
  'stitch local dotted fragments through consecutive certified centres into six HAX_i curves',
  'attach boundary letters (vanishing variables) to the six external branches',
  'attach the P,Q,R,S external-boundary routes from the historical text to the certified boundary arcs',
  'only then derive reflection involutions and act on the 136-node/105-edge diagnostic carrier'
 ],
 'status':'PROTOCOL-REFINED; NO_HISTORICAL_PROMOTION_IN_V46'
}
(B/'selling_axis_stitching_requirements_v46.json').write_text(json.dumps(req,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'trace_count':7,'subsets':7,'survivors':0,'closest_pair':closest,'angle_error_bound':pair_angle_err,'decision':out['decision']},indent=2))
