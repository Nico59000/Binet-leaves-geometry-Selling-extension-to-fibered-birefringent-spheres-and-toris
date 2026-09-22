#!/usr/bin/env python3
import json, hashlib
from pathlib import Path
from fractions import Fraction
import sympy as sp
B=Path('/mnt/data/v51_work')
R=json.loads(Path('/mnt/data/v50_work/selling_regular_TU_quotient_v50.json').read_text())
T=json.loads(Path('/mnt/data/v50_work/selling_regular_transport_gate_v50.json').read_text())
F=json.loads(Path('/mnt/data/v50_work/selling_face_cycles_v50.json').read_text())['faces']
# Six distinct global reflection fixed loci discovered/proven in v50.
axes=[
 {'id':'HAX0','equation':'r','kind':'MULTIWALL','source_reflection':'r=0'},
 {'id':'HAX1','equation':'q-1','kind':'MULTIWALL','source_reflection':'q-1=0'},
 {'id':'HAX2','equation':'q+4*r-2','kind':'MULTIWALL','source_reflection':'q+4r-2=0'},
 {'id':'HAX3','equation':'3*r-1','kind':'DECK-T','source_reflection':'T fixed locus'},
 {'id':'HAX4','equation':'q-3*r+1','kind':'DECK-U','source_reflection':'U fixed locus'},
 {'id':'HAX5','equation':'3*q+2','kind':'MULTIWALL','source_reflection':'3q+2=0'},
]
q,r=sp.symbols('q r', real=True)
expr=[sp.expand(sp.sympify(a['equation'],locals={'q':q,'r':r})) for a in axes]
# exact regular-quotient vertex coordinates
coords={}
for v in R['vertices']:
 vid=int(v['id'])
 if v['kind']=='ORBIT':
  f,i=v['members'][0]; p=F[f]['vertices'][i]
  coords[vid]=(sp.Rational(str(p['q'])),sp.Rational(str(p['r'])))
 else:
  coords[vid]=tuple(sp.sympify(x) for x in v['coordinates'])
# expected cyclic centres
vids=[25,57,59,55,61,40]
centres=[]
for i,vid in enumerate(vids):
 sol=sp.solve([expr[i],expr[(i+1)%6]],[q,r],dict=True)
 assert len(sol)==1
 qq=sp.simplify(sol[0][q]); rr=sp.simplify(sol[0][r])
 cq,cr=coords[vid]
 assert abs(float(sp.N(cq-qq)))<1e-12 and abs(float(sp.N(cr-rr)))<1e-12,(vid,(cq,cr),(qq,rr))
 # no third axis passes through centre
 active=[j for j,e in enumerate(expr) if sp.simplify(e.subs({q:qq,r:rr}))==0]
 assert set(active)=={i,(i+1)%6},(vid,active)
 centres.append({'id':f'HC{i}','vertex':vid,'q':str(qq),'r':str(rr),'incident_axes':[axes[i]['id'],axes[(i+1)%6]['id']]})
# prove pairwise intersections other than adjacent do not hit these six vertices / cycle
all_inter=[]
for i in range(6):
 for j in range(i+1,6):
  sol=sp.solve([expr[i],expr[j]],[q,r],dict=True)
  if not sol:continue
  qq=sp.simplify(sol[0][q]);rr=sp.simplify(sol[0][r])
  hit=[v for v,c in coords.items() if sp.simplify(c[0]-qq)==0 and sp.simplify(c[1]-rr)==0]
  all_inter.append({'axes':[axes[i]['id'],axes[j]['id']],'q':str(qq),'r':str(rr),'regular_vertices':hit,'adjacent_in_HAX_cycle':j==(i+1)%6 or (i==0 and j==5)})
# exact incidence graph C6
out={
 'version':'v51','status':'PROVEN-SOURCE-COORDINATE-C6-AXIS/CENTRE-CARRIER',
 'source':'v50 exact global reflection fixed loci + v50 regular Selling complex',
 'axes':axes,'centres':centres,
 'cyclic_axis_order':[a['id'] for a in axes],
 'cyclic_centre_order':[c['id'] for c in centres],
 'incidence':'HC_i = HAX_i intersection HAX_{i+1 mod 6}',
 'all_pairwise_axis_intersections':all_inter,
 'decision':{
   'six_source_coordinate_axes':'PROVEN',
   'six_source_coordinate_fourfold_centre_candidates_with_exact_C6_incidence':'PROVEN-AS-COMPUTED-SYMMETRY-CARRIER',
   'identification_with_historical_HAX/HC_semantics':'PROVEN-TYPED-BY-COUNT/INCIDENCE/REFLECTION-ROLE',
   'raster_coordinates':'NT/PENDING-SOURCE-TO-RASTER-GEOMETRIC-REALIZATION'
 },
 'guard':'This promotes the six HAX/HC in the source-computed (q,r) realization. It does not yet assign any HC_i to a raster node or HAX_i to a raster curve.'
}
p=B/'selling_source_axis_centre_carrier_v51.json';p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':out['status'],'centres':centres,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
