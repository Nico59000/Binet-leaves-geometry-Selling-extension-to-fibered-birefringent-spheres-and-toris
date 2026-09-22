#!/usr/bin/env python3
import json,pathlib,collections,sympy as sp,hashlib
B=pathlib.Path(__file__).resolve().parent
reg=json.load(open(B/'selling_regular_TU_quotient_v50.json'));carrier=json.load(open(B/'selling_source_axis_centre_carrier_v51.json'));bfs=json.load(open('/mnt/data/v50_roundtrip/exact_TU_field_bfs_v50.json'))['records'];Q=sp.Matrix([[7,0,2],[0,-1,0],[2,0,-8]])
edge_by_id={int(e['id']):e for e in reg['edges']};v_edges=collections.defaultdict(list);v_faces=collections.defaultdict(set)
for e in reg['edges']:
 a,b=e['endpoints'];v_edges[a].append(e);v_edges[b].append(e)
for fs,cyc in reg['face_cycles'].items():
 f=int(fs)
 for item in cyc:
  e=edge_by_id[int(item['edge'])]
  for v in e['endpoints']:v_faces[v].add(f)
def form(f):
 A=sp.Matrix([[int(x) for x in row] for row in bfs[f]['path_matrix']]);M=A.T*Q*A;return [[int(M[i,j]) for j in range(3)] for i in range(3)]
out=[]
for c in carrier['centres']:
 v=int(c['vertex']);edges=sorted(v_edges[v],key=lambda e:int(e['id']));faces=sorted(v_faces[v]);out.append({'centre':c['id'],'vertex':v,'q':c['q'],'r':c['r'],'incident_axes':c['incident_axes'],'incident_edge_count':len(edges),'incident_edges':[{'id':int(e['id']),'kind':e['kind'],'endpoints':e['endpoints']} for e in edges],'incident_faces':faces,'incident_face_count':len(faces),'incident_form_matrices':{str(f):form(f) for f in faces}})
cert={'version':'v51','status':'PROVEN-SOURCE-HC-LOCAL-STAR-SIGNATURES','source':'v50 regular Selling complex + v51 exact source C6 axis/centre carrier','centres':out,'decision':'THESE_SIX_EXACT_LINK_SIGNATURES_ARE_MANDATORY_TARGETS_FOR_ANY_RASTER_HC_BINDING','guard':'No raster node is promoted here. Raster assignment requires matching the local incidence/axis roles and stabilizers, not visual proximity alone.'}
p=B/'selling_source_HC_local_star_signatures_v51.json';p.write_text(json.dumps(cert,indent=2,sort_keys=True)+'\n');print(json.dumps(cert,indent=2));print('sha256',hashlib.sha256(p.read_bytes()).hexdigest())
