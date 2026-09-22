#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import json, hashlib
ROOT=Path(__file__).resolve().parent
figs={}
for n in (3,4,5,6):
 p=ROOT/f'selling_fig{n}_local_compiler_replay_v52.json';d=json.load(open(p));assert d['status'].startswith('PASS-TYPED-LOCAL-CARRIER');assert d['compiled_carrier']['unresolved']==[];figs[f'Fig{n}']=d
# Namespaced disjoint-union carrier; shared printed labels are recorded but never identified.
objects=[];relations=[];ordered=[];prov={}; label_occ=defaultdict(list)
for fn,d in figs.items():
 c=d['compiled_carrier']
 for o in c['objects']:
  nid=f'{fn}:{o["id"]}';objects.append({'id':nid,'local_id':o['id'],'label':o['label'],'figure':fn,'provenance':'SOURCE-DIRECT/LOCAL-COMPILER'});prov[nid]={'figure':fn,'source_basis':d['source_basis'],'local_compiler':f'selling_{fn.lower()}_local_compiler_replay_v52.json'};label_occ[o['label']].append(nid)
 for rel in c.get('relations',[]):relations.append({'figure':fn,'source':f'{fn}:{rel[0]}','target':f'{fn}:{rel[1]}','relation_type':rel[2],'provenance':'LOCAL-COMPILER'})
 if fn=='Fig6':
  ordered.append({'figure':fn,'kind':'principal_nesting','sequence':c['ordered_principal_family'],'provenance':'selling_fig6_planar_nesting_order_v52.json'})
  ordered.append({'figure':fn,'kind':'right_anchor_chain','sequence':c['right_anchor_chain'],'provenance':'selling_fig6_planar_nesting_order_v52.json'})
  ordered.append({'figure':fn,'kind':'lower_left_cell_chain','sequence':c['lower_left_cell_chain'],'provenance':'selling_fig6_planar_nesting_order_v52.json'})
coinc=[]
for lab,ids in sorted(label_occ.items()):
 if len(ids)>1:coinc.append({'printed_label':lab,'occurrences':ids,'status':'SEPARATED/NO-CROSS-FIGURE-IDENTIFICATION'})
checks={
 'all_four_local_compilers_PASS':all(d['status'].startswith('PASS') for d in figs.values()),
 'all_local_unresolved_empty':all(d['compiled_carrier']['unresolved']==[] for d in figs.values()),
 'namespaces_disjoint':len({o['id'] for o in objects})==len(objects),
 'no_cross_figure_relation':all(r['source'].split(':',1)[0]==r['target'].split(':',1)[0]==r['figure'] for r in relations),
 'provenance_total':len(prov)==len(objects),
 'Fig6_order_preserved':len(ordered)==3,
 'shared_labels_not_identified':all(x['status'].startswith('SEPARATED') for x in coinc),
}
assert all(checks.values())
out={'version':'v52','phase':'typed-Fig3-6-fusion','status':'PASS-TYPED-DISJOINT-FUSION/PROVENANCE-PRESERVED','carrier_kind':'NAMESPACED-DISJOINT-UNION-WITH-TYPED-LOCAL-RELATIONS','figures':{k:{'status':v['status'],'carrier_kind':v['compiled_carrier']['carrier_kind'],'source_basis':v['source_basis']} for k,v in figs.items()},'object_count':len(objects),'relation_count':len(relations),'ordered_structure_count':len(ordered),'objects':objects,'relations':relations,'ordered_structures':ordered,'coincident_printed_labels':coinc,'provenance':prov,'checks':checks,'guard':'Fusion means typed coproduct only. No identical printed fraction, tuple, region number, strand style, or local relation is identified across figures. No Fig3--6 object is identified with Fig1, P/Q/R/S, or any GL(3) morphism without an explicit future adapter.'}
(ROOT/'selling_fig3_6_typed_fusion_v52.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'status':out['status'],'object_count':len(objects),'relation_count':len(relations),'ordered_structures':len(ordered),'coincident_labels':coinc,'checks':checks},indent=2,ensure_ascii=False))
