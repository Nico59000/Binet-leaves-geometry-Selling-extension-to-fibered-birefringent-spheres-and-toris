#!/usr/bin/env python3
from __future__ import annotations
import json, sys, hashlib
from pathlib import Path
import sympy as sp
HERE=Path('/mnt/data/v50_work'); sys.path.insert(0,str(HERE))
from selling_global_wall_linearization_v50 import GT,GU
q,r=sp.symbols('q r')

E=json.load(open(HERE/'selling_edge_orbits_v50.json'))
V=json.load(open(HERE/'selling_vertex_orbits_v50.json'))
F=json.load(open(HERE/'selling_face_cycles_v50.json')); faces={x['field']:x for x in F['faces']}
BFS=json.load(open(HERE/'exact_TU_field_bfs_v50.json')); recs={x['id']:x for x in BFS['records']}
VC=json.load(open(HERE/'selling_face_vertex_exact_certificates_v50.json'))
certs=VC['certificates']
cert_by_local={(c['field'],c['local_vertex_id']):c for c in certs}
TU=GT*GU
DECK={'I':sp.eye(3),'T':GT,'U':GU,'TU':TU}

def cycle_dart(fid,si,ci=0):
    f=faces[fid]; cy=f['cycles'][ci]; s=cy['segments'][si]; co=f['components'][s['component']]
    return {'field':fid,'cycle':ci,'si':si,'v':tuple(s['v']),'factor':sp.Poly(sp.sympify(co['factor']),q,r,domain=sp.QQ).as_expr(),
            'factor_str':co['factor'],'kind':co['kind'],'walls':tuple(co['walls'])}

def proj_formula(D):
    vv=sp.Matrix([1,-r,q]); u=D.T*vv
    den=sp.expand(u[0])
    return sp.cancel(u[2]/den), sp.cancel(-u[1]/den), den

def prim(poly):
    P=sp.Poly(sp.expand(poly),q,r,domain=sp.QQ)
    den, Pz = P.clear_denoms(convert=True)
    P=Pz
    cont,pr=sp.polys.polytools.primitive(P.as_expr(),q,r)
    P=sp.Poly(pr,q,r,domain=sp.ZZ)
    if P.LC()<0: P=-P
    return P

def curve_map_proof(src_factor,tgt_factor,D):
    qp,rp,den=proj_formula(D)
    expr=sp.together(tgt_factor.subs({q:qp,r:rp}, simultaneous=True))
    num,dd=expr.as_numer_denom()
    Ps=prim(src_factor); Pn=prim(num)
    quo,rem=sp.div(Pn,Ps,domain=sp.ZZ)
    ok=rem.is_zero and quo.total_degree()==0 and quo.as_expr()!=0
    return {'ok':bool(ok),'source_primitive':str(Ps.as_expr()),'target_pullback_primitive':str(Pn.as_expr()),
            'quotient':str(quo.as_expr()) if rem.is_zero else None,'projective_denominator':str(den),
            'q_prime':str(qp),'r_prime':str(rp)}

def source_param(c):
    t=sp.symbols('t')
    if c['representation_variable']=='rational':
        qq=sp.Rational(c['coordinates_exact'][0]); rr=sp.Rational(c['coordinates_exact'][1])
        return {'rational':True,'q':qq,'r':rr,'t':None,'P':None,'interval':None}
    rv=c['representation_variable']; sym=q if rv=='q' else r
    Pfull=sp.Poly(sp.sympify(c['root_polynomial']).subs(sym,t),t,domain=sp.QQ)
    oth=sp.sympify(c['other_coordinate_expr']).subs(sym,t)
    qq=t if rv=='q' else oth; rr=oth if rv=='q' else t
    iv=[sp.Rational(x) for x in c['root_interval']]
    # Groebner elimination polynomials can contain extraneous factors belonging to other intersections.
    # Select the unique irreducible factor carrying the certified isolated root interval.
    facs=sp.factor_list(Pfull.as_expr(),t)[1]
    carriers=[]
    for fe,mult in facs:
        Pf=sp.Poly(fe,t,domain=sp.QQ)
        if iv[0]==iv[1]:
            if Pf.eval(iv[0])==0: carriers.append(Pf)
        else:
            for riv,rm in sp.intervals(Pf,eps=sp.Rational(1,10)**30):
                if not (riv[1] < iv[0] or iv[1] < riv[0]):
                    carriers.append(Pf); break
    # de-duplicate equal factors
    uniq=[]
    for z in carriers:
        if all(z.monic()!=u.monic() for u in uniq): uniq.append(z)
    if len(uniq)!=1:
        raise RuntimeError(('root-factor-selection',c['field'],c['local_vertex_id'],str(Pfull.as_expr()),c['root_interval'],[str(z.as_expr()) for z in uniq]))
    P=uniq[0].monic()
    return {'rational':False,'q':sp.cancel(qq),'r':sp.cancel(rr),'t':t,'P':P,'Pfull':Pfull,'interval':iv,'rep':rv}

def zero_mod(expr,param):
    expr=sp.cancel(expr)
    if param['rational']:
        return bool(sp.simplify(expr)==0), {'num_rem':'0' if sp.simplify(expr)==0 else str(sp.simplify(expr)),'den_gcd_degree':0}
    num,den=sp.together(expr).as_numer_denom(); t=param['t']; P=param['P']
    Num=sp.Poly(sp.expand(num),t,domain=sp.QQ); Den=sp.Poly(sp.expand(den),t,domain=sp.QQ)
    rem=Num.rem(P)
    gd=sp.gcd(Den,P)
    return bool(rem.is_zero and gd.degree()==0), {'num_rem':str(rem.as_expr()),'den_gcd_degree':int(gd.degree())}

def point_map_proof(src_cert,tgt_cert,D):
    ps=source_param(src_cert)
    qp0,rp0,den=proj_formula(D)
    qp=sp.cancel(qp0.subs({q:ps['q'],r:ps['r']}, simultaneous=True))
    rp=sp.cancel(rp0.subs({q:ps['q'],r:ps['r']}, simultaneous=True))
    if tgt_cert['representation_variable']=='rational':
        tq=sp.Rational(tgt_cert['coordinates_exact'][0]); tr=sp.Rational(tgt_cert['coordinates_exact'][1])
        okq,pq=zero_mod(qp-tq,ps); okr,pr=zero_mod(rp-tr,ps); ok=okq and okr
        return {'ok':ok,'target_kind':'rational','q_check':pq,'r_check':pr}
    rv=tgt_cert['representation_variable']; sym=q if rv=='q' else r
    tcoord=qp if rv=='q' else rp; other=rp if rv=='q' else qp
    P2=sp.sympify(tgt_cert['root_polynomial']).subs(sym,tcoord)
    E2=sp.sympify(tgt_cert['other_coordinate_expr']).subs(sym,tcoord)
    okP,pP=zero_mod(P2,ps); okE,pE=zero_mod(other-E2,ps)
    # isolate identity: transformed numeric approx must fall in the target root interval in the representative coordinate.
    # This is only root selection after exact algebraic vanishing; target interval itself is exact-rational.
    qa,ra=src_cert['approx']; vv=sp.Matrix([1,-sp.Float(ra,60),sp.Float(qa,60)]); uu=D.T*vv
    qn=sp.N(uu[2]/uu[0],50); rn=sp.N(-uu[1]/uu[0],50); tn=qn if rv=='q' else rn
    ia,ib=[sp.Rational(x) for x in tgt_cert['root_interval']]
    inside = (sp.N(ia,50)-sp.Float('1e-40') <= tn <= sp.N(ib,50)+sp.Float('1e-40')) if ia!=ib else abs(float(tn-sp.N(ia,50)))<1e-12
    return {'ok':bool(okP and okE),'target_kind':rv,'root_equation':pP,'other_coordinate_equation':pE,
            'target_root_interval':tgt_cert['root_interval'],'interval_selection':bool(inside)}

def prove_link(src_key,tgt_key,D):
    sd=cycle_dart(*src_key); td=cycle_dart(*tgt_key)
    cp=curve_map_proof(sd['factor'],td['factor'],D)
    if not cp['ok']:
        return {'ok':False,'curve':cp}
    sverts=sd['v']; tverts=td['v']
    # test both endpoint pairings exactly
    pairings=[[(sverts[0],tverts[0]),(sverts[1],tverts[1])],[(sverts[0],tverts[1]),(sverts[1],tverts[0])]]
    results=[]
    for pairs in pairings:
        pp=[]; good=True
        for sv,tv in pairs:
            pr=point_map_proof(cert_by_local[(src_key[0],sv)], cert_by_local[(tgt_key[0],tv)], D)
            good &= pr['ok']; pp.append({'source_vertex':sv,'target_vertex':tv,'proof':pr})
        results.append({'ok':bool(good),'pairs':pp})
    oks=[i for i,x in enumerate(results) if x['ok']]
    return {'ok':len(oks)==1,'curve':cp,'endpoint_pairings':results,'orientation':'same' if oks==[0] else ('reversed' if oks==[1] else 'ambiguous/fail')}

# Build exact generator links over canonical cycle0 occurrences.
ordinary_certs=[]; direct_certs=[]; failures=[]
for x in E['ordinary_transport_links']:
    src=tuple(x['from']); tgt=tuple(x['to']); meta=x['meta']
    if meta.get('new'): D=sp.eye(3)
    else:
        a,b=meta['deck_exponents']; D=(GT**a)*(GU**b)
    # If numeric matching first landed in lobe B of field55, canonicalize B -> A by TU.
    if x.get('raw_target_cycle')==1 and tgt[0]==55: D=D*TU
    pr=prove_link(src,tgt,D)
    rec={'source':list(src),'target':list(tgt),'wall':x['wall'],'deck_matrix':[[str(z) for z in row] for row in D.tolist()],
         'raw_target_cycle':x.get('raw_target_cycle',0),'proof':pr}
    ordinary_certs.append(rec)
    if not pr['ok']: failures.append(('ordinary',src,tgt,pr))

for x in E['direct_deck_links']:
    src=tuple(x['from']); tgt=tuple(x['to']); D=DECK[x['deck']]
    if x.get('raw_cycle')==1 and tgt[0]==55: D=D*TU
    pr=prove_link(src,tgt,D)
    rec={'source':list(src),'target':list(tgt),'deck':x['deck'],'deck_matrix':[[str(z) for z in row] for row in D.tolist()],
         'raw_target_cycle':x.get('raw_cycle',0),'proof':pr}
    direct_certs.append(rec)
    if not pr['ok']: failures.append(('deck',src,tgt,pr))

# Exact lobeB -> lobeA edge mapping under TU, independently of prior numeric map.
lobe_map=[]
for si in range(len(faces[55]['cycles'][1]['segments'])):
    sd=cycle_dart(55,si,1)
    good=[]
    for tj in range(len(faces[55]['cycles'][0]['segments'])):
        td=cycle_dart(55,tj,0)
        cp=curve_map_proof(sd['factor'],td['factor'],TU)
        if cp['ok']: good.append((tj,cp))
    lobe_map.append({'source_cycle1_segment':si,'targets':[x[0] for x in good],'proofs':[x[1] for x in good]})
    if len(good)!=1: failures.append(('lobe-factor',si,good))

# Rebuild orbit DSU only from exact-certified generator links.
allkeys=[]
for fid in range(62):
    for si in range(len(faces[fid]['cycles'][0]['segments'])): allkeys.append((fid,si))
parent={k:k for k in allkeys}
def find(x):
    while parent[x]!=x:
        parent[x]=parent[parent[x]]; x=parent[x]
    return x
def union(a,b):
    a,b=find(a),find(b)
    if a!=b: parent[max(a,b)]=min(a,b)
for rec in ordinary_certs+direct_certs:
    if rec['proof']['ok']: union(tuple(rec['source']),tuple(rec['target']))
orb={}
for k in allkeys:orb.setdefault(find(k),[]).append(k)
exact_orbits=sorted([sorted(v) for v in orb.values()],key=lambda x:x[0])
expected=sorted([[tuple(y) for y in O] for O in E['edge_orbits']],key=lambda x:x[0])
orbits_match=(exact_orbits==expected)
if not orbits_match: failures.append(('orbit_partition_mismatch',len(exact_orbits),len(expected)))

# Loop orientation exact: locate direct U self-identification and require reversed endpoint pairing.
loop_certs=[]
for ee in V['edge_endpoints']:
    if not ee['loop']: continue
    ei=ee['edge']; rep=tuple(ee['representative'])
    candidates=[x for x in direct_certs if tuple(x['source'])==rep and tuple(x['target'])==rep and x.get('deck')=='U']
    if not candidates:
        # maybe another member has U self-link
        for k in map(tuple,E['edge_orbits'][ei]):
            candidates += [x for x in direct_certs if tuple(x['source'])==k and tuple(x['target'])==k and x.get('deck')=='U']
    good=[x for x in candidates if x['proof']['ok']]
    orient=None if not good else good[0]['proof']['orientation']
    status='PASS-EXACT-ORIENTATION-REVERSING' if orient=='reversed' else 'FAIL'
    loop_certs.append({'edge_orbit':ei,'representative':list(rep),'deck':'U','orientation':orient,'tangent_sign':-1 if orient=='reversed' else (1 if orient=='same' else None),'status':status})
    if status!='PASS-EXACT-ORIENTATION-REVERSING': failures.append(('loop',ei,orient))

out={
 'version':'v50','status':'PASS-EXACT' if not failures else 'FAIL','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
 'edge_occurrences':len(allkeys),'edge_orbits':len(exact_orbits),'ordinary_generator_links':len(ordinary_certs),'direct_deck_generator_links':len(direct_certs),
 'orbit_partition_matches_prior_candidate':orbits_match,'lobe55_cycle1_to_cycle0_under_TU':lobe_map,
 'loop_orientations':loop_certs,'failures':failures,
 'ordinary_transport_certificates':ordinary_certs,'direct_deck_certificates':direct_certs,
 'guard':'Exact polynomial pullback + exact algebraic endpoint transport. No H1/H2 computed here.'
}
p=HERE/'selling_edge_orbit_exact_certificates_v50.json'
p.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'status':out['status'],'edge_orbits':out['edge_orbits'],'failures':len(failures),'loop_orientations':loop_certs,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()},indent=2))
