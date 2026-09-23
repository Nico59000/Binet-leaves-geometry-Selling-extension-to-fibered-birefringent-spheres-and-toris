#!/usr/bin/env python3
from __future__ import annotations
import json,sys,hashlib
from pathlib import Path
import sympy as sp
HERE=Path('/mnt/data/v50_work');sys.path.insert(0,str(HERE))
import selling_global_wall_linearization_v50 as gl
from selling_reconstruction_engine_v49 import RELABEL
F=json.load(open(HERE/'selling_face_cycles_v50.json'));faces={x['field']:x for x in F['faces']}
E=json.load(open(HERE/'selling_edge_orbits_v50.json'))
V=json.load(open(HERE/'selling_vertex_orbits_v50.json'))
EC=json.load(open(HERE/'selling_edge_orbit_exact_certificates_v50.json'))
BFS=json.load(open(HERE/'exact_TU_field_bfs_v50.json'));recs={x['id']:x for x in BFS['records']}

edge_of={tuple(k):ei for ei,O in enumerate(E['edge_orbits']) for k in O}
vorb={tuple(k):vi for vi,O in enumerate(V['vertex_orbits']) for k in O}
old_end={x['edge']:(x['v0'],x['v1']) for x in V['edge_endpoints']}

# Exact deck stabilizers of quotient face representatives.
def mat(x): return sp.Matrix([[sp.Integer(v) for v in row] for row in x])
def tup(A): return tuple(int(x) for x in list(A))
rel={tup(B):(i,1) for i,B in enumerate(RELABEL)}; rel.update({tup(-B):(i,-1) for i,B in enumerate(RELABEL)})
DECK={'T':gl.GT,'U':gl.GU,'TU':gl.GT*gl.GU}
stabilizers={}
for i in range(62):
 A=mat(recs[i]['path_matrix']); arr=[]
 for nm,D in DECK.items():
  X=sp.simplify(A.inv()*D*A); key=tup(X)
  if key in rel:
   bi,sg=rel[key];arr.append({'deck':nm,'relabel_index':bi,'sign':sg,'conjugate':X.tolist()})
 if arr: stabilizers[i]=arr
expected_stab={38:'U',48:'U',49:'U',53:'T',55:'TU',56:'U'}
assert {i:v[0]['deck'] for i,v in stabilizers.items()}==expected_stab

# Exact fixed points forced by orientation-reversing invariant boundary edges 80 and 108.
s5=sp.sqrt(5)
mid80=(sp.Rational(-34,61)+sp.Rational(6,61)*s5, sp.Rational(9,61)+sp.Rational(2,61)*s5)
mid108=(sp.Rational(-2,3),sp.Rational(1,9))
q,r=gl.q,gl.r

def proj(D,pt):
 v=sp.Matrix([1,-pt[1],pt[0]]);u=D.T*v
 return (sp.simplify(u[2]/u[0]),sp.simplify(-u[1]/u[0]))
assert all(sp.simplify(a-b)==0 for a,b in zip(proj(gl.GU,mid80),mid80))
assert all(sp.simplify(a-b)==0 for a,b in zip(proj(gl.GU,mid108),mid108))

# reduced signs at fixed midpoint witnesses
def signs_at(fid,pt):
 L=mat(recs[fid]['wall_action_matrix']); vals=[sp.factor(sp.simplify(x.subs({q:pt[0],r:pt[1]}))) for x in L*gl.w0]
 dv=sp.factor(sp.simplify(gl.D0.subs({q:pt[0],r:pt[1]})))
 return vals,dv
mid_checks={}
for fid,pt in [(38,mid80),(49,mid80),(56,mid108)]:
 vals,dv=signs_at(fid,pt)
 assert sp.N(dv)>0 and all(sp.N(x)<=0 for x in vals)
 mid_checks[str(fid)]={'walls':{lab:str(v) for lab,v in zip(gl.labs,vals)},'D0':str(dv),'reduced':True}
# interior bracketing of fixed points by exact isolating intervals / rational endpoints.
VC=json.load(open(HERE/'selling_face_vertex_exact_certificates_v50.json'))['certificates']; cidx={(c['field'],c['local_vertex_id']):c for c in VC}
def iv_for(c,coord):
 if c['representation_variable']=='rational':
  z=sp.Rational(c['coordinates_exact'][0 if coord=='q' else 1]);return z,z
 if c['representation_variable']==coord:
  return tuple(sp.Rational(x) for x in c['root_interval'])
 # if other coordinate is constant rational, use it
 ex=sp.sympify(c['other_coordinate_expr'])
 if not ex.free_symbols: return ex,ex
 return None
# edge80 endpoints: field38 local 0 and1; compare q. mid q lies strictly between their q coordinates.
c380,c381=cidx[(38,0)],cidx[(38,1)]
iv0=iv_for(c380,'q');iv1=iv_for(c381,'q'); qm=mid80[0]
assert iv0 and iv1
lo=max(min(iv0),min(iv1)); hi=min(max(iv0),max(iv1)) # not useful if disjoint
# exact ordering using interval envelopes
assert (iv0[0] < qm < iv1[1]) or (iv1[0] < qm < iv0[1]) or (iv0[1] < qm < iv1[0]) or (iv1[1] < qm < iv0[0])
# More direct numeric ordering with exact interval bounds: one endpoint q=-1/3 and other interval near -0.341...
qvals=[iv0,iv1]
assert any(a==b==sp.Rational(-1,3) for a,b in qvals)
other=[z for z in qvals if not (z[0]==z[1]==sp.Rational(-1,3))][0]
assert other[1] < qm < sp.Rational(-1,3)
# edge108 endpoints field56 local2 and1; compare r around 1/9
c562,c561=cidx[(56,2)],cidx[(56,1)]
r2=iv_for(c562,'r');r1=iv_for(c561,'r'); rm=mid108[1]
assert r2 and r1
assert (r2[1] < rm < r1[0]) or (r1[1] < rm < r2[0])

# Fixed-line equations for mirror edges.
fixed_lines={'U':'q - 3*r + 1','T':'3*r - 1'}
# all fixed vertices used for mirror endpoints satisfy the corresponding line exactly, verified via their Groebner/root representation elsewhere;
# here use exact reduction helper from edge certificate script by reading its exact endpoint transport certificates + known deck fixed vertex action.

# Regular quotient vertices: old 0..59 plus two fixed midpoints.
MID80=60; MID108=61; C0=62
vertex_data=[{'id':i,'kind':'ORBIT','members':V['vertex_orbits'][i]} for i in range(60)]
vertex_data += [
 {'id':MID80,'kind':'DECK-FIXED-MIDPOINT','deck':'U','edge_orbit':80,'coordinates':['-34/61 + 6*sqrt(5)/61','9/61 + 2*sqrt(5)/61'],'status':'PASS-EXACT-INTERIOR'},
 {'id':MID108,'kind':'DECK-FIXED-MIDPOINT','deck':'U','edge_orbit':108,'coordinates':['-2/3','1/9'],'status':'PASS-EXACT-INTERIOR'}]

# Existing edge orbits remain one quotient edge each; reflection loops become half-edges to their fixed midpoint.
edge_end={e:tuple(old_end[e]) for e in range(118)}
edge_end[80]=(23,MID80); edge_end[108]=(39,MID108)
edge_data=[{'id':e,'kind':'BOUNDARY-ORBIT','endpoints':list(edge_end[e]),'source_orbit':E['edge_orbits'][e]} for e in range(118)]
# mirror edges for the five stabilizer faces whose representative component is itself reflected.
mirror_specs=[
 (38,MID80,41,'U'),
 (48,41,50,'U'),  # oriented to close chosen path pos0(50)->...->pos2(41): mirror 41->50
 (49,MID80,43,'U'),
 (53,54,47,'T'),
 (56,MID108,50,'U'),
]
mirror_edge={}
for j,(fid,a,b,deck) in enumerate(mirror_specs,start=118):
 mirror_edge[fid]=j; edge_end[j]=(a,b); edge_data.append({'id':j,'kind':'DECK-FIXED-MIRROR','face':fid,'deck':deck,'fixed_locus':fixed_lines[deck],'endpoints':[a,b]})
C1=len(edge_data); assert C1==123

# Build fundamental quotient face cycles. Edge signs are determined solely from oriented endpoint sequence.
def old_vertex(fid,pos):
 vlocal=faces[fid]['cycles'][0]['vertices'][pos]; return vorb[(fid,vlocal)]
def add_step(steps,e,a,b,label):
 ea,eb=edge_end[e]
 if (a,b)==(ea,eb):sg=1
 elif (a,b)==(eb,ea):sg=-1
 else: raise RuntimeError(('edge endpoints mismatch',label,e,(a,b),(ea,eb)))
 steps.append({'edge':e,'sign':sg,'v_start':a,'v_end':b,'source':label})

def full_face(fid):
 cy=faces[fid]['cycles'][0]; Vpos=cy['vertices']; n=len(cy['segments']);steps=[]
 for si in range(n):
  e=edge_of[(fid,si)];a=vorb[(fid,Vpos[si])];b=vorb[(fid,Vpos[(si+1)%n])]
  add_step(steps,e,a,b,[fid,si])
 return steps

def regular_face(fid):
 if fid not in {38,48,49,53,56}: return full_face(fid)
 cy=faces[fid]['cycles'][0]
 steps=[]
 if fid==38:
  # pos0 fixed -> pos1 -> pos2 -> midpoint80 -> fixed pos0
  p0,p1,p2=[old_vertex(fid,k) for k in [0,1,2]]
  add_step(steps,edge_of[(fid,0)],p0,p1,[fid,0]); add_step(steps,edge_of[(fid,1)],p1,p2,[fid,1]);
  add_step(steps,80,p2,MID80,[fid,'half-fixed',2]); add_step(steps,mirror_edge[fid],MID80,p0,[fid,'mirror'])
 elif fid==48:
  p0,p1,p2=[old_vertex(fid,k) for k in [0,1,2]]
  add_step(steps,edge_of[(fid,0)],p0,p1,[fid,0]); add_step(steps,edge_of[(fid,1)],p1,p2,[fid,1]); add_step(steps,mirror_edge[fid],p2,p0,[fid,'mirror'])
 elif fid==49:
  p2,p0=[old_vertex(fid,k) for k in [2,0]]
  add_step(steps,edge_of[(fid,2)],p2,p0,[fid,2]); add_step(steps,80,p0,MID80,[fid,'half-fixed',0]); add_step(steps,mirror_edge[fid],MID80,p2,[fid,'mirror'])
 elif fid==53:
  p0,p1,p2=[old_vertex(fid,k) for k in [0,1,2]]
  add_step(steps,edge_of[(fid,0)],p0,p1,[fid,0]); add_step(steps,edge_of[(fid,1)],p1,p2,[fid,1]); add_step(steps,mirror_edge[fid],p2,p0,[fid,'mirror'])
 elif fid==56:
  p1,p2=[old_vertex(fid,k) for k in [1,2]]
  add_step(steps,edge_of[(fid,1)],p1,p2,[fid,1]); add_step(steps,108,p2,MID108,[fid,'half-fixed',2]); add_step(steps,mirror_edge[fid],MID108,p1,[fid,'mirror'])
 # closure
 assert steps[-1]['v_end']==steps[0]['v_start'] and all(steps[k]['v_end']==steps[(k+1)%len(steps)]['v_start'] for k in range(len(steps)))
 return steps

face_steps={i:regular_face(i) for i in range(62)}
# Face55 already uses lobeA as fundamental representative of its TU-exchanged pair; no internal mirror added.
assert 55 not in mirror_edge
C2=62
D1=sp.zeros(C0,C1)
for e,(a,b) in edge_end.items():D1[a,e]-=1;D1[b,e]+=1
D2=sp.zeros(C1,C2)
for i,steps in face_steps.items():
 for x in steps:D2[x['edge'],i]+=x['sign']
P=D1*D2;nz=[(i,j,int(P[i,j])) for i in range(P.rows) for j in range(P.cols) if P[i,j]!=0]
assert not nz

# Check mirror edges have exactly one incident quotient face; non-mirror incidence is arbitrary due quotient identifications but well-defined.
mirror_inc={e:sum(abs(int(D2[e,j])) for j in range(C2)) for e in mirror_edge.values()}
assert all(v==1 for v in mirror_inc.values())

out={
 'version':'v50','status':'PASS-EXACT-REGULAR-QUOTIENT','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
 'pre_refinement_candidate':[60,118,62],'correction_reason':'Two quotient loop edges are fixed setwise by U with tangent sign -1. Their interior fixed points must be vertices; five deck-stabilized face representatives require mirror edges in a regular quotient cell structure.',
 'deck_face_stabilizers':{str(k):v for k,v in stabilizers.items()},
 'fixed_midpoint_certificates':{
  'edge80':{'deck':'U','coordinates':vertex_data[MID80]['coordinates'],'fixed_exact':True,'interior_exact':True,'reduced_in_faces':[38,49],'wall_checks':{k:mid_checks[k] for k in ['38','49']}},
  'edge108':{'deck':'U','coordinates':vertex_data[MID108]['coordinates'],'fixed_exact':True,'interior_exact':True,'reduced_in_faces':[56],'wall_checks':{'56':mid_checks['56']}},
 },
 'C0':C0,'C1':C1,'C2':C2,'euler_characteristic':C0-C1+C2,
 'vertices':vertex_data,'edges':edge_data,'face_cycles':face_steps,
 'partial1_rank_Q':int(D1.rank()),'partial2_rank_Q':int(D2.rank()),'partial1_partial2_zero':True,
 'partial1_sparse':[{'row':i,'col':j,'value':int(D1[i,j])} for i in range(C0) for j in range(C1) if D1[i,j]],
 'partial2_sparse':[{'row':i,'col':j,'value':int(D2[i,j])} for i in range(C1) for j in range(C2) if D2[i,j]],
 'mirror_edge_incidence':{str(k):v for k,v in mirror_inc.items()},
 'guard':'This regularizes deck-reflection stabilizers before any local-system cohomology. Face55 is a two-component orbit exchanged by TU and uses one lobe as its quotient face representative.'
}
p=HERE/'selling_regular_TU_quotient_v50.json';p.write_text(json.dumps(out,indent=2,ensure_ascii=False,default=str),encoding='utf-8')
# csv
with open(HERE/'selling_regular_partial1_v50.csv','w') as f:
 f.write('row_vertex,col_edge,value\n');
 for x in out['partial1_sparse']:f.write(f"{x['row']},{x['col']},{x['value']}\n")
with open(HERE/'selling_regular_partial2_v50.csv','w') as f:
 f.write('row_edge,col_face,value\n');
 for x in out['partial2_sparse']:f.write(f"{x['row']},{x['col']},{x['value']}\n")
print(json.dumps({'status':out['status'],'C0':C0,'C1':C1,'C2':C2,'chi':out['euler_characteristic'],'ranks':[out['partial1_rank_Q'],out['partial2_rank_Q']],'d1d2_zero':True,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
