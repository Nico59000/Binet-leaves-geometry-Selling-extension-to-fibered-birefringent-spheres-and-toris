#!/usr/bin/env python3
import json,sys,hashlib
from pathlib import Path
import sympy as sp
HERE=Path('/mnt/data/v50_work')
F=json.load(open(HERE/'selling_face_cycles_v50.json'));faces={x['field']:x for x in F['faces']}
E=json.load(open(HERE/'selling_edge_orbits_v50.json'))
V=json.load(open(HERE/'selling_vertex_orbits_v50.json'))
# mappings
edge_of={}
for ei,O in enumerate(E['edge_orbits']):
 for k in O:edge_of[tuple(k)]=ei
vorb={}
for vi,O in enumerate(V['vertex_orbits']):
 for k in O:vorb[tuple(k)]=vi
edge_end={x['edge']:(x['v0'],x['v1']) for x in V['edge_endpoints']}
C0=V['vertex_orbit_count'];C1=E['edge_orbit_count'];C2=62
# d1 boundary C1->C0
D1=sp.zeros(C0,C1)
for e,(a,b) in edge_end.items():
 D1[a,e]-=1;D1[b,e]+=1
# d2 C2->C1 based on each quotient representative face cycle0
D2=sp.zeros(C1,C2);occ=[];loop_pending=[]
for i in range(62):
 cy=faces[i]['cycles'][0];verts=cy['vertices'];n=len(cy['segments'])
 for si,s in enumerate(cy['segments']):
  e=edge_of[(i,si)];sv=vorb[(i,verts[si])];tv=vorb[(i,verts[(si+1)%n])];a,b=edge_end[e]
  if a!=b:
   if (sv,tv)==(a,b):sgn=1
   elif (sv,tv)==(b,a):sgn=-1
   else:raise RuntimeError(('endpoint mismatch',i,si,e,(sv,tv),(a,b)))
  else:
   # edge is quotient loop; representative occurrence orientation +1, other occurrences unresolved tangent sign.
   rep=tuple(E['edge_orbits'][e][0]);sgn=1 if (i,si)==rep else 1
   loop_pending.append({'face':i,'segment':si,'edge':e,'representative':list(rep),'provisional_sign':sgn})
  D2[e,i]+=sgn
  occ.append({'face':i,'segment':si,'edge':e,'sign':sgn,'v_start':sv,'v_end':tv,'loop':a==b})
P=D1*D2
nz=[(i,j,int(P[i,j])) for i in range(P.rows) for j in range(P.cols) if P[i,j]!=0]
print('shapes',D1.shape,D2.shape,'nnz D1',sum(1 for x in D1 if x),'nnz D2',sum(1 for x in D2 if x),'product nz',len(nz))
print('ranks Q',D1.rank(),D2.rank(),'loop pending',len(loop_pending),'loop edges',sorted(set(x['edge'] for x in loop_pending)))
# sparse serializations
def sparse(M):return [{'row':i,'col':j,'value':int(M[i,j])} for i in range(M.rows) for j in range(M.cols) if M[i,j]!=0]
out={'version':'v50','status':'PASS_BOUNDARY_SQUARED_ZERO' if not nz else 'FAIL','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE','C0':C0,'C1':C1,'C2':C2,'euler_characteristic':C0-C1+C2,
     'partial1_shape':[C0,C1],'partial2_shape':[C1,C2],'partial1_rank_Q':D1.rank(),'partial2_rank_Q_provisional_loop_signs':D2.rank(),'partial1_partial2_zero':not nz,'product_nonzero':nz,
     'partial1_sparse':sparse(D1),'partial2_sparse':sparse(D2),'face_edge_occurrences':occ,'loop_orientation_pending':loop_pending,
     'guard':'Loop-edge tangent orientation signs remain to be source-certified before integral H1/H2 or twisted cohomology; they do not affect partial1*partial2 because loop columns of partial1 vanish.'}
raw=json.dumps(out,indent=2).encode();out['sha256_pre_field']=hashlib.sha256(raw).hexdigest()
(HERE/'selling_global_chain_complex_v50.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
# csv matrices
with open(HERE/'selling_partial1_v50.csv','w') as f:
 f.write('row_vertex,col_edge,value\n');
 for x in sparse(D1):f.write(f"{x['row']},{x['col']},{x['value']}\n")
with open(HERE/'selling_partial2_v50.csv','w') as f:
 f.write('row_edge,col_face,value\n');
 for x in sparse(D2):f.write(f"{x['row']},{x['col']},{x['value']}\n")
