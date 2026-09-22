#!/usr/bin/env python3
import json, sympy as sp
from collections import defaultdict,Counter
from pathlib import Path
B=Path('/mnt/data/v50_work')
reg=json.load(open(B/'selling_regular_TU_quotient_v50.json'))
g=json.load(open(B/'selling_field_orbigroupoid_v50.json'))
ar={a['id']:a for a in g['arrows']}
def M(x): return sp.Matrix([[sp.Rational(v) for v in row] for row in x])
# occurrence lookup edge
edgeocc=defaultdict(list)
for fs,cyc in reg['face_cycles'].items():
    f=int(fs)
    for i,z in enumerate(cyc): edgeocc[z['edge']].append((f,i,z))
# transition: cross current edge then turn across target face corner
trans={}; crossed={}
for fstr,cyc in reg['face_cycles'].items():
    f=int(fstr); n=len(cyc)
    for i,z in enumerate(cyc):
        e=z['edge']
        for v in [z['v_start'],z['v_end']]:
            occs=edgeocc[e]
            if len(occs)==2:
                target=[o for o in occs if not (o[0]==f and o[1]==i)][0]
            elif len(occs)==1: target=occs[0]
            else: raise RuntimeError((e,len(occs)))
            tf,ti,tz=target
            assert v in [tz['v_start'],tz['v_end']]
            tcyc=reg['face_cycles'][str(tf)]; tn=len(tcyc)
            oi=(ti-1)%tn if v==tz['v_start'] else (ti+1)%tn
            oz=tcyc[oi]; assert v in [oz['v_start'],oz['v_end']]
            trans[(f,i,v)]=(tf,oi,v)
            # direction crossing current edge relative to chosen arrow
            aa=ar[e]
            if aa['source']==f and aa['target']==tf: sgn=1
            elif aa['source']==tf and aa['target']==f: sgn=-1
            elif aa['source']==aa['target']==f==tf: sgn=1
            else: raise RuntimeError(('arrow mismatch',e,f,tf,aa['source'],aa['target']))
            crossed[(f,i,v)]={'edge':e,'from_face':f,'to_face':tf,'direction':sgn}
# extract cycles; interior gives 2 inverse cycles, mirror/orbifold 1 doubled cycle. choose lexicographic canonical among cycles by edge+face sequence
links=[]
for v in range(reg['C0']):
    flags=[k for k in trans if k[2]==v]; seen=set(); cycles=[]
    for x in sorted(flags):
        if x in seen: continue
        cur=x; arr=[]
        while cur not in seen:
            seen.add(cur); arr.append(cur); cur=trans[cur]
        assert cur==x
        cycles.append(arr)
    words=[]
    for arr in cycles:
        word=[crossed[x] for x in arr]
        # holonomy
        hb=sp.Integer(1); ha=sp.eye(3)
        for z in word:
            aa=ar[z['edge']]; eps=sp.Integer(aa['binet_transport']); R=M(aa['ad_transport'])
            if z['direction']==1: hb*=eps; ha=sp.simplify(R*ha)
            else: hb*=1/eps; ha=sp.simplify(R.inv()*ha)
        words.append({'flags':[list(x) for x in arr],'word':word,'binet_holonomy':str(sp.simplify(hb)),
                      'ad_holonomy':[[str(ha[i,j]) for j in range(3)] for i in range(3)],'ad_identity':ha==sp.eye(3)})
    # choose canonical: shortest lex key among cycles; if two, they should be reverses and holonomy inverses
    def key(w): return tuple((z['edge'],z['from_face'],z['to_face'],z['direction']) for z in w['word'])
    chosen=min(words,key=key)
    links.append({'vertex':v,'cycle_count':len(cycles),'cycle_lengths':[len(x) for x in cycles], 'chosen':chosen,'all_cycles':words})
# output summary
hout=Counter((x['chosen']['binet_holonomy'],x['chosen']['ad_identity']) for x in links)
print('cycle count hist',Counter(x['cycle_count'] for x in links),'length patterns',Counter(tuple(x['cycle_lengths']) for x in links))
print('holonomy summary',hout)
for x in links:
    if x['chosen']['binet_holonomy']!='1' or not x['chosen']['ad_identity']:
        print('nontrivial',x['vertex'],'len',x['cycle_lengths'],'B',x['chosen']['binet_holonomy'],'A',x['chosen']['ad_holonomy'],'edges',[z['edge'] for z in x['chosen']['word']])
out={'version':'v50','status':'PASS-EXACT-COMBINATORIAL-LINK-EXTRACTION','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
     'vertex_count':62,'links':links,'cycle_count_histogram':dict(Counter(x['cycle_count'] for x in links)),
     'guard':'Link words are exact combinatorial cycles of the regular primal quotient. Nontrivial coefficient holonomy, if present, must be treated as local orbifold stabilizer data; it is not silently set to identity.'}
json.dump(out,open(B/'selling_orbifold_links_v50.json','w'),indent=2)
