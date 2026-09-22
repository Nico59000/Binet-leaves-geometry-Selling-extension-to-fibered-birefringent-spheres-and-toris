#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import json,csv,re
ROOT=Path(__file__).resolve().parent
arc={e['id']:e for e in json.load(open(ROOT/'fig1_v50.arc_geometry.v03.auto.json'))['edges']}
old={int(r['edge_id']):r for r in csv.DictReader(open(ROOT/'selling_global_curvilinear_123edge_audit_v52.csv',encoding='utf-8'))}
exact={int(r['edge_id']):r for r in csv.DictReader(open(ROOT/'selling_exact_triangulated_BSB_P084_123edge_audit_K12_v52.csv',encoding='utf-8'))}
hax_chains=json.load(open(ROOT/'selling_global_curvilinear_123edge_certificate_v52.json'))['hax_exact_edge_chains']
hax={i for xs in hax_chains.values() for i in xs}
inc=defaultdict(list)
for i,e in arc.items():
    for v in e['endpoints']:inc[v].append(i)
def norm(s):return re.sub(r'\s+','',str(s))
pairs=[]
for v,ids in inc.items():
    for ai,a in enumerate(ids):
        for b in ids[ai+1:]:
            if norm(arc[a]['implicit_factor'])==norm(arc[b]['implicit_factor']):pairs.append((a,b,v))
def close(R):
    R=set(R);trace=[]
    while True:
        change=False;n=set();w=[]
        for a,b,v in pairs:
            if a in R and b not in R:n.add(b);w.append({'new':b,'from':a,'vertex':v,'factor':arc[b]['implicit_factor']})
            if b in R and a not in R:n.add(a);w.append({'new':a,'from':b,'vertex':v,'factor':arc[a]['implicit_factor']})
        if n:
            R|=n;trace.append({'rule':'SAME-EXACT-IMPLICIT-FACTOR-CONTINUATION','new_edges':sorted(n),'witnesses':[x for x in w if x['new'] in n]});change=True
        n=set();w=[]
        for i,e in arc.items():
            if i in R:continue
            for v in e['endpoints']:
                if all(j in R or j==i for j in inc[v]):
                    n.add(i);w.append({'new':i,'vertex':v,'full_link':sorted(inc[v]),'already_resolved':sorted(j for j in inc[v] if j in R)});break
        if n:
            R|=n;trace.append({'rule':'UNIQUE-MISSING-EXACT-LINK-INCIDENCE','new_edges':sorted(n),'witnesses':[x for x in w if x['new'] in n]});change=True
        if not change:return R,trace
# Old checkpoint 53 -> source closure 117.
old_bsb={i for i,r in old.items() if r['final_positive_status']!='NT-GEOMETRIC-REALIZATION'}
base53=old_bsb|hax
assert len(base53)==53
R117,tr117=close(base53)
res6=sorted(set(range(123))-R117)
assert len(R117)==117 and res6==[34,53,54,78,96,97],(len(R117),res6)
# Exact P084 seed must be strong in final exact triangulation.
e96=exact[96]
assert e96['status']=='STRONG-EXACT-TRIANGULATED-STROKE-SUPPORT'
R123,trseed=close(R117|{96})
assert len(R123)==123 and not (set(range(123))-R123)
# The seed closure must contain exactly the other five residual edges.
seed_added=sorted(R123-(R117|{96}))
assert seed_added==[34,53,54,78,97],seed_added
# Compact rule trace expected from prior gate.
compact=[(x['rule'],x['new_edges']) for x in trseed]
# final provenance table: must be valid in the final exact P084 triangulated realization itself.
exact_pos={i for i,r in exact.items() if r['status']!='NT-GEOMETRIC-REALIZATION'}
final_base=exact_pos|hax
Rfinal,trfinal=close(final_base)
assert len(Rfinal)==123
final_prop=Rfinal-final_base
seed_prop=set(seed_added)
rows=[]
for i in range(123):
    pv=[]
    if i in exact_pos: pv.append(exact[i]['status'])
    if i in hax: pv.append('PROVEN-SOURCE-HAX-INCIDENCE')
    if i in final_prop: pv.append('PROVEN-SOURCE-EXACT-PROPAGATION-IN-FINAL-P084-REALIZATION')
    if i==96: pv.append('FINAL-SEED-BSB-P084-STRONG-EXACT')
    if i in seed_prop: pv.append('PROVEN-SOURCE-EXACT-PROPAGATION-FROM-e96-COUNTERFACTUAL-117-GATE')
    assert pv
    rows.append({'edge_id':i,'kind':arc[i]['kind'],'walls':'|'.join(arc[i].get('walls') or []),'endpoints':str(arc[i]['endpoints']),'final_historical_incidence':'PROVEN','provenance':' ; '.join(pv)})
with open(ROOT/'selling_fig1_final_123_historical_incidence_v52.csv','w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
cert={
 'version':'v52','phase':'final-Fig1-gate','status':'PASS-HISTORICAL-FIG1-123/123',
 'input_old_checkpoint':{'bsb_positive_plus_HAX_count':len(base53),'source_exact_closure_count':len(R117),'residual_edges_before_P084':res6},
 'source_rules':{
  'same_factor':'At a common exact source vertex, equal normalized exact implicit factor continues the same source boundary carrier.',
  'unique_link':'In the exact regular source link, once every incident edge except one has historical incidence, the remaining labelled incidence is forced.',
  'guard':'These rules transport already source-typed v50 incidences only; they do not infer a new edge from raster morphology.'},
 'preseed_closure_trace':tr117,
 'final_exact_base_count':len(final_base),
 'final_exact_closure_trace':trfinal,
 'P084_exact_seed':{'edge':96,'status':e96['status'],'geom_support_eps':float(e96['geom_support_eps']),'dark_support_eps':float(e96['dark_support_eps']),'P084_support_eps':float(e96['P084_support_eps']),'certificate':'selling_exact_triangulated_BSB_P084_homeomorphism_certificate_K12_v52.json'},
 'seed_closure_trace':trseed,
 'seed_implication':{'statement':'e96 => {e53,e34,e97,e54,e78} by exact source propagation','new_edges':seed_added,'compact_trace':compact},
 'final':{'historical_edge_orbits_proven':123,'historical_edge_orbits_total':123,'remaining_NT':0,'promotion':'PROVEN-HISTORICAL-GLOBAL-FIG1-REALIZATION'},
 'provenance_csv':'selling_fig1_final_123_historical_incidence_v52.csv',
 'guard':'BSB-P084 supplies only the e96 seed through the one global exact triangulated realization. The other residual incidences are source-exact consequences. Partial BSB support on e78/e97 is retained as corroborating provenance, not used as the logical seed.'}
(ROOT/'selling_fig1_final_123_certificate_v52.json').write_text(json.dumps(cert,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'status':cert['status'],'base53':len(base53),'preseed117':len(R117),'residual6':res6,'e96':cert['P084_exact_seed'],'seed_added':seed_added,'seed_trace':compact,'final':cert['final']},indent=2,ensure_ascii=False))
