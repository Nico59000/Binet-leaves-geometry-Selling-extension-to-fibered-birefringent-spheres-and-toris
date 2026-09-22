#!/usr/bin/env python3
import json,hashlib,subprocess,sys
from pathlib import Path
B=Path('/mnt/data/v50_work')
files=['selling_four_sheet_deck_cover_v50.json','selling_F29_strict_cover_functor_v50.json','selling_F29_descent_square_v50.json']
def h(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
before={f:h(B/f) for f in files}
subprocess.run([sys.executable,str(B/'build_four_sheet_deck_cover_v50.py')],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
after={f:h(B/f) for f in files}
C=json.load(open(B/files[0]));F=json.load(open(B/files[1]));D=json.load(open(B/files[2]))
checks=[]
def ck(name,cond):
 if not cond: raise AssertionError(name)
 checks.append(name)
ck('cover_status',C['status'].startswith('PROVEN-FOUR-SHEET'))
ck('cover_248_objects',C['counts']['cover_objects']==248)
ck('cover_492_arrows',C['counts']['cover_arrows']==492)
ck('twelve_base_deck_loops',C['counts']['base_nontrivial_deck_loops']==12)
ck('fortyeight_intersheet_arrows',C['counts']['lifted_intersheet_arrows']==48)
ck('cover_248_based_2relations',C['counts']['cover_based_2relations']==248)
ck('cover_242_geometric_2cells',C['counts']['cover_geometric_2cells']==242)
ck('strict_upper_status',F['status']=='PASS-STRICT-UPPER-FUNCTOR-TO-F29-CAN')
ck('strict_248_objects',F['checks']['objects_total']==248)
ck('strict_492_arrows',F['checks']['arrows_strict_endpoint_matches']==492)
ck('strict_248_2relations',F['checks']['based_2relations_strictly_closed']==248)
ck('fillable_248_2relations',F['checks']['based_2relations_fillable_by_841_fundamental_cells']==248)
ck('basis_841',F['checks']['F29_fundamental_cell_count']==841)
ck('descent_status',D['status']=='PASS-DESCENT-SQUARE-OBJECT/ARROW/2RELATION')
ck('quotient_objects_62',D['checks']['cover_object_orbits']==62)
ck('quotient_arrows_123',D['checks']['cover_arrow_orbits']==123)
ck('quotient_2relations_62',D['checks']['cover_based_2relation_orbits']==62)
ck('target_orbits_42',D['checks']['target_object_orbits']==42)
ck('square_object_arrow_2relation',D['checks']['object_square_commutes'] and D['checks']['arrow_square_commutes'] and D['checks']['upper_functor_strict_objects_arrows_2relations'])
ck('three_nontrivial_links',[(x['vertex'],x['holonomy']) for x in D['nontrivial_base_links']]==[(48,'TU'),(59,'T'),(61,'U')])
ck('deterministic_replay',before==after)
out={'version':'v50','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE','status':'PASS',
     'checks_passed':len(checks),'checks_total':len(checks),'checks':checks,'sha256':after,
     'cohomology_gate':'RECONSIDERABLE BUT NOT OPENED: first audit V4-equivariant pullback/descent naturality of the Binet and Ad local systems.'}
(B/'v50_four_sheet_descent_verification.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(f"PASS {len(checks)}/{len(checks)}")
for f,s in after.items(): print(f,s)
