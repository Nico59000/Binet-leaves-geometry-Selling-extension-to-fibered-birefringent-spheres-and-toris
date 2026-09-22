#!/usr/bin/env python3
import json, sympy as sp
from collections import defaultdict, Counter
from pathlib import Path
B=Path('/mnt/data/v50_work')
reg=json.load(open(B/'selling_regular_TU_quotient_v50.json'))
fc=json.load(open(B/'selling_face_cycles_v50.json'))
bfs=json.load(open(B/'exact_TU_field_bfs_v50.json'))
glob=json.load(open(B/'selling_global_wall_linearization_v50.json'))
tg=json.load(open(B/'selling_regular_transport_gate_v50.json'))
records={r['id']:r for r in bfs['records']}
# metadata for source face cycle segments
meta={}
for fo in fc['faces']:
    cyc=fo['cycles'][0]
    for si,seg in enumerate(cyc['segments']):
        comp=fo['components'][seg['component']]
        meta[(fo['field'],si)]={'kind':comp['kind'],'walls':comp['walls'],'factor':comp['factor'],'component':seg['component']}
# incidence of quotient edges in quotient faces
inc=defaultdict(list)
for fs,cyc in reg['face_cycles'].items():
    f=int(fs)
    for z in cyc: inc[z['edge']].append({'face':f,'sign':z['sign'],'source':z['source']})
# base Lie algebra basis for Q2
M0=sp.Matrix([[7,0,2],[0,-1,0],[2,0,-8]])
E=[sp.Matrix([[0,1,0],[7,0,2],[0,0,0]]),sp.Matrix([[0,4,0],[30,0,0],[0,1,0]]),sp.Matrix([[-2,0,8],[0,0,0],[7,0,2]])]
Bmat=sp.Matrix.hstack(*[x.reshape(9,1) for x in E])
def mat(x): return sp.Matrix([[sp.Rational(v) for v in row] for row in x])
def ad_matrix(A):
    cols=[]
    for X in E:
        Y=sp.simplify(A.inv()*X*A)
        sol=sp.linsolve((Bmat,Y.reshape(9,1)))
        vv=sp.Matrix(list(next(iter(sol))))
        assert Bmat*vv==Y.reshape(9,1)
        cols.append(vv)
    return sp.Matrix.hstack(*cols)
GT=mat(glob['deck_in_Q2_chart']['G_T']); GU=mat(glob['deck_in_Q2_chart']['G_U'])
I3=sp.eye(3)
decks={(0,0):I3,(1,0):GT,(0,1):GU,(1,1):GT*GU}
deck_name={(0,0):'I',(1,0):'T',(0,1):'U',(1,1):'TU'}
# exact assertions
assert GT*GT==I3 and GU*GU==I3 and GT*GU==GU*GT
ADdeck={k:ad_matrix(A) for k,A in decks.items()}
# cert maps
mw={x['edge_orbit']:x for x in tg['multiwall_isotropy_certificates']}
mir={x['edge_id']:x for x in tg['mirror_isotropy_certificates']}
# ordinary edge metadata helper
arrows=[]
fail=[]
for e in reg['edges']:
    eid=e['id']; ii=inc[eid]
    rec={'id':eid,'primal_endpoints':e['endpoints'],'primal_kind':e['kind']}
    if eid in mir:
        c=mir[eid]; f=c['face']; d=(1,0) if c['deck']=='T' else (0,1) if c['deck']=='U' else (1,1)
        rec.update(kind='DECK-MIRROR-ISOTROPY',source=f,target=f,wall=None,deck=c['deck'],deck_exponents=list(d),projective_sign=1,
                   binet_transport=int(c['binet_time_cone_sign']),ad_transport=ADdeck[d],field_matrix=mat(c['field_stabilizer_K']))
    elif eid in mw:
        c=mw[eid]; f=c['field']; J=mat(c['global_congruence_J'])
        rec.update(kind='MULTIWALL-ISOTROPY',source=f,target=f,wall='+'.join(c['walls']),deck='I',deck_exponents=[0,0],projective_sign=1,
                   binet_transport=int(c['binet_time_cone_sign']),ad_transport=ad_matrix(J),field_matrix=mat(c['field_stabilizer_K']),factor=c['factor'])
    else:
        # ordinary crossing: choose first source occurrence whose metadata is ordinary
        candidates=[]
        for src in e['source_orbit']:
            if len(src)>=2 and isinstance(src[1],int):
                mm=meta.get((src[0],src[1]))
                if mm and mm['kind']=='ORDINARY': candidates.append((tuple(src),mm))
        if not candidates:
            fail.append((eid,'no ordinary source')); continue
        (f,si),mm=candidates[0]; wall=mm['walls'][0]
        n=records[f]['neighbors'][wall]
        t=int(n['target'])
        # Verify incidence target if two-face edge, or self if single ordinary
        faces=sorted({x['face'] for x in ii})
        if len(faces)==2 and t not in faces:
            fail.append((eid,'target mismatch',f,wall,t,faces)); continue
        if len(faces)==1 and t!=f:
            fail.append((eid,'single target mismatch',f,wall,t,faces)); continue
        if n.get('new',False): d=(0,0); sig=1; ridx=0
        else:
            d=tuple(int(x)%2 for x in n.get('deck_exponents',[0,0])); sig=int(n.get('sign',1)); ridx=int(n.get('relabel_index',0))
        eps=sig*((-1)**(d[0]+d[1]))
        rec.update(kind='ORDINARY-CROSSING',source=f,target=t,wall=wall,deck=deck_name[d],deck_exponents=list(d),projective_sign=sig,
                   relabel_index=ridx,binet_transport=int(eps),ad_transport=ADdeck[d],source_segment=[f,si],factor=mm['factor'])
    # stringify matrices later
    arrows.append(rec)
assert not fail, fail
assert len(arrows)==123
# validate transports and counts
assert Counter(a['kind'] for a in arrows)==Counter({'ORDINARY-CROSSING':102,'MULTIWALL-ISOTROPY':16,'DECK-MIRROR-ISOTROPY':5})
# d0 unrestricted (for diagnostic/local link relation dimensions)
D0B=sp.zeros(123,62)
D0A=sp.zeros(123*3,62*3)
for a in arrows:
    e=a['id']; s=a['source']; t=a['target']; eps=sp.Integer(a['binet_transport']); R=a['ad_transport']
    D0B[e,t]+=1; D0B[e,s]-=eps
    D0A[3*e:3*e+3,3*t:3*t+3]+=sp.eye(3)
    D0A[3*e:3*e+3,3*s:3*s+3]-=R
# local left kernels at each primal vertex
local=[]
for v in range(reg['C0']):
    es=[e['id'] for e in reg['edges'] if v in e['endpoints']]
    fs=sorted(set(sum(([a['source'],a['target']] for a in arrows if a['id'] in es),[])))
    bmat=D0B.extract(es,fs)
    bdim=len(bmat.T.nullspace())
    erows=[]
    for e in es: erows += [3*e+i for i in range(3)]
    fcols=[]
    for f in fs: fcols += [3*f+i for i in range(3)]
    amat=D0A.extract(erows,fcols)
    adim=len(amat.T.nullspace())
    local.append({'vertex':v,'edges':es,'faces':fs,'binet_local_relation_dim_unrestricted':bdim,'ad_local_relation_dim_unrestricted':adim})
# serialize
for a in arrows:
    if isinstance(a.get('ad_transport'),sp.MatrixBase): a['ad_transport']=[[str(a['ad_transport'][i,j]) for j in range(3)] for i in range(3)]
    if isinstance(a.get('field_matrix'),sp.MatrixBase): a['field_matrix']=[[str(a['field_matrix'][i,j]) for j in range(3)] for i in range(3)]
out={'version':'v50','status':'PASS-EXACT-DUAL-ARROWS/LOCAL-RELATION-DIAGNOSTIC','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
     'objects':62,'arrows':arrows,'arrow_counts':dict(Counter(a['kind'] for a in arrows)),
     'base_soQ_basis':[[[str(x) for x in row] for row in E0.tolist()] for E0 in E],
     'deck_ad':{deck_name[k]:[[str(ADdeck[k][i,j]) for j in range(3)] for i in range(3)] for k in decks},
     'd0_unrestricted_ranks':{'Binet':D0B.rank(),'Ad':D0A.rank()},'local_relation_diagnostics':local,
     'guards':['unrestricted d0 is diagnostic only until stabilizer/orientation modules are imposed','H1/H2 forbidden']}
json.dump(out,open(B/'selling_field_orbigroupoid_v50.json','w'),indent=2)
print('arrow counts',out['arrow_counts'])
print('d0 ranks',out['d0_unrestricted_ranks'])
print('Binet local dims',Counter(x['binet_local_relation_dim_unrestricted'] for x in local))
print('Ad local dims',Counter(x['ad_local_relation_dim_unrestricted'] for x in local))
