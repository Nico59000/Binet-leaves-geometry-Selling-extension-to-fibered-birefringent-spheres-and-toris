#!/usr/bin/env python3
from pathlib import Path
import json, hashlib, sys
B=Path(__file__).resolve().parent
ROOT=B.parent / 'base' if (B.parent/'base').exists() else B
checks={}
def ck(k,v): checks[k]=bool(v); assert v,k
def J(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
f=J(B/'selling_fig2_6_common_typed_functor_v53.json')
fig2_path=(B/'selling_fig2_local_compiler_v52.json') if (B/'selling_fig2_local_compiler_v52.json').exists() else (B/'05_SELLING/v53/dependencies/selling_fig2_local_compiler_v52.json')
fig2=J(fig2_path)
fig36=J(ROOT/'selling_fig3_6_typed_fusion_v52.json')
ck('status',f['status']=='PROVEN-EXPLICIT-COPRODUCT-FUNCTOR-FIG2-6/PROVENANCE-PRESERVED')
ck('decision',f['decision']=='PASS')
ck('source_hash_fig2',sha(fig2_path)==f['source_hashes']['selling_fig2_local_compiler_v52.json'])
ck('source_hash_fig36',sha(ROOT/'selling_fig3_6_typed_fusion_v52.json')==f['source_hashes']['selling_fig3_6_typed_fusion_v52.json'])
objs=f['objects']; morph=f['morphisms']; ords=f['ordered_structures']
ck('object_count',len(objs)==f['object_count']==42)
ck('morphism_count',len(morph)==f['morphism_count']==43)
ck('ordered_count',len(ords)==f['ordered_structure_count']==3)
ids=[o['id'] for o in objs]; ck('object_ids_unique',len(ids)==len(set(ids)))
mid=[m['id'] for m in morph]; ck('morphism_ids_unique',len(mid)==len(set(mid)))
S=set(ids)
ck('all_endpoints_exist',all(m['source'] in S and m['target'] in S for m in morph))
ck('no_cross_figure_morphism',all(m['source'].split(':')[1]==m['figure']==m['target'].split(':')[1] for m in morph))
ck('all_five_figures',set(o['figure'] for o in objs)=={'Fig2','Fig3','Fig4','Fig5','Fig6'})
ck('all_five_adapters',set(f['adapters'])=={'Fig2','Fig3','Fig4','Fig5','Fig6'})
ck('no_fig1_import',all(not a['imports_Fig1'] for a in f['adapters'].values()))
ck('fig2_source_recipe',all(o.get('source_specific_recipe') is True for o in objs if o['figure']=='Fig2'))
ck('fig2_compiler_pass',fig2['status']=='PASS-TYPED-LOCAL-CARRIER')
ck('fig36_fusion_pass',fig36['status'].startswith('PASS-TYPED-DISJOINT-FUSION'))
ck('provenance_present',all(o.get('provenance') for o in objs) and all(m.get('provenance') for m in morph))
# No literal collision is used as an identification.
for c in f['coincident_literal_labels']:
    ck('collision_'+c['label'],c['status']=='SEPARATED/NO-CROSS-FIGURE-IDENTIFICATION' and len({x.split(':')[1] for x in c['objects']})>1)
ck('numeric_seed_gap_typed',all(a['numeric_seed_matrix'].startswith('NT-') for a in f['adapters'].values()))
ck('retained_R2PURE',f['retained_separations']['R2PURE13_to_SellingDeck17']=='SEPARATED/PENDING-MORPHISM')
ck('retained_Noether',f['retained_separations']['Selling_Noether']=='CONTEXT/PENDING-MORPHISM')
ck('retained_3D5D',f['retained_separations']['Selling_3D_to_5D']=='CONTEXT/PENDING-MORPHISM')
out={'version':'v53','status':'PASS','passed':sum(checks.values()),'total':len(checks),'checks':checks}
(B/'v53_SELLING_FUNCTOR_REPLAY_REPORT.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2,sort_keys=True))
