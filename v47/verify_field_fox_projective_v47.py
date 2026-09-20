#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys,json,csv,hashlib,collections
import sympy as sp
B=Path(__file__).resolve().parent
checks={}
def ok(k,v):
 checks[k]=bool(v)
 if not v: print('FAIL',k,file=sys.stderr)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
# Parent v46 closed + byte-exact sealed ledger.
tok=json.loads((B/'v46_ATOMIC_COMPLETION_TOKEN.json').read_text(encoding='utf-8'))
ok('parent_v46_closed',tok.get('publication_state')=='CLOSED/APPEND-ONLY')
ok('parent_v46_atomic',tok.get('atomic_completion') is True and all(v=='PASS' for v in tok['atomic_conjunction'].values()))
ok('parent_v46_checks_8034',tok['verification']['combined']==8034 and tok['verification']['status']=='PASS_8034_OF_8034')
entries=[];mis=[]
for line in (B/'v46_SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
 if not line.strip():continue
 h,n=line.split('  ',1);entries.append(n);p=B/n
 if not p.exists() or sha(p)!=h:mis.append(n)
ok('parent_v46_ledger_276',len(entries)==276)
ok('parent_v46_ledger_byte_exact',not mis)
ok('parent_v46_pdf_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v46.pdf')==tok['pdf']['sha256'])
ok('parent_v46_tex_hash',sha(B/'formalisation_feuillets_binet_birefringence_selling_polaire_v46.tex')==tok['tex']['sha256'])
# A. Source-curve transcription / DAX graph no-go.
p=subprocess.run([sys.executable,str(B/'analyze_selling_transcription_v47.py')],cwd=B,capture_output=True,text=True)
ok('sell47_transcription_helper_exit_zero',p.returncode==0)
tr=json.loads((B/'selling_axis_curve_transcription_v47.json').read_text(encoding='utf-8'))
routes=json.loads((B/'selling_PQRS_route_semantics_v47.json').read_text(encoding='utf-8'))
c2=json.loads((B/'selling_degree1_C2_gate_v47.json').read_text(encoding='utf-8'))
ok('sell47_abstract_C6_proven',tr['historical_abstract_axis_centre_carrier']['status']=='PROVEN_SOURCE_SEMANTIC_INCIDENCE')
ok('sell47_six_HC',len(tr['historical_abstract_axis_centre_carrier']['vertices'])==6)
ok('sell47_six_HAX',len(tr['historical_abstract_axis_centre_carrier']['edges'])==6)
ok('sell47_dax_7',len(tr['diagnostic_DAX_graph']['vertices'])==7)
ok('sell47_dax_edges_7',len(tr['diagnostic_DAX_graph']['edges'])==7)
ok('sell47_dax00_isolated',tr['diagnostic_DAX_graph']['isolated_traces']==['DAX00'])
ok('sell47_no_dax_C6',tr['diagnostic_DAX_graph']['has_six_cycle'] is False)
ok('sell47_only_4_5_cycles',set(map(int,tr['diagnostic_DAX_graph']['simple_cycle_lengths'].keys()))=={4,5})
ok('sell47_raster_centres_unset',tr['raster_HC_coordinates'] is None)
ok('sell47_raster_axes_unset',tr['raster_HAX_curves'] is None)
ok('sell47_status_nt',tr['status'].startswith('NT/'))
ok('sell47_PQRS_semantics_proven',routes['status'].startswith('PROVEN_SOURCE_ROUTE_SEMANTICS'))
ok('sell47_P_raster_nt',routes['routes']['P']['raster_polyline'] is None)
ok('sell47_Q_raster_nt',routes['routes']['Q']['raster_polyline'] is None)
ok('sell47_TU_relation_retained',routes['relations']['TU_equals_UT']=='PROVEN_SOURCE_RELATION')
ok('sell47_C2_not_activated',c2['activated'] is False and c2['fixed_points'] is None)
ok('sell47_89_count',c2['degree_one_candidates']==89)
# H1/H2 remains fail-closed. Materialize explicit v47 gate.
nodes=list(csv.DictReader(open(B/'selling_fig1_nodes_v40.csv',encoding='utf-8')))
segs=list(csv.DictReader(open(B/'selling_fig1_segments_v40.csv',encoding='utf-8')))
fields=list(csv.DictReader(open(B/'selling_fig1_fields_v40.csv',encoding='utf-8')))
missing={
 'historical_node_assignments':sum(r['historical_node_type']=='UNRESOLVED' or not r['historical_node_id'] for r in nodes),
 'historical_segment_assignments':sum(r['historical_assignment_status'].startswith('NT_') for r in segs),
 'historical_field_assignments':sum(r['field_assignment_status'].startswith('NT_') for r in fields),
 'ordered_boundary_cycles':sum(not r['ordered_boundary_cycle'] or r['ordered_boundary_cycle']=='UNRESOLVED' for r in fields),
 'repetition_pairs':sum(not r['repetition_of'] or r['repetition_of']=='UNRESOLVED' for r in fields),
 'stabilizer_words':sum(not r['stabilizer_word'] or r['stabilizer_word']=='UNRESOLVED' for r in fields)}
ok('sell47_C0_136_missing',len(nodes)==136 and missing['historical_node_assignments']==136)
ok('sell47_C1_105_missing',len(segs)==105 and missing['historical_segment_assignments']==105)
ok('sell47_C2_34_missing',len(fields)==34 and missing['historical_field_assignments']==34)
ok('sell47_boundary_34_missing',missing['ordered_boundary_cycles']==34)
ok('sell47_repeat_34_missing',missing['repetition_pairs']==34)
ok('sell47_stabilizer_34_missing',missing['stabilizer_words']==34)
h12={
 'version':'v47','table_sizes':{'C0':136,'C1':105,'C2':34},'missing':missing,
 'source_semantic_axis_centre_C6':'PROVEN','source_raster_HC_HAX_embedding':'NT',
 'PQRS_route_semantics':'PROVEN_SOURCE_TEXT','PQRS_raster_arc_binding':'NT',
 'historical_reflections':0,'C2_89_census':'NOT_ACTIVATED_FAIL_CLOSED',
 'd0':'NOT_FORMED','d1':'NOT_FORMED','d1d0_test':'NOT_FORMED','face_holonomy_stabilizer_test':'NOT_FORMED',
 'F29_object_edge_face_functor':'PENDING-MORPHISM','H1_hist':'NT_NOT_RECOMPUTED','H2_hist':'NT_NOT_RECOMPUTED',
 'decision':'ABSTRACT SOURCE INCIDENCE IS KNOWN, BUT THE RASTER/CELL TRANSCRIPTION IS NOT TOTAL; COCHAIN MATRICES REMAIN FAIL-CLOSED.'}
(B/'selling_h12_materialization_gate_v47.json').write_text(json.dumps(h12,indent=2,ensure_ascii=False),encoding='utf-8')
ok('sell47_h12_blocked',h12['H1_hist'].startswith('NT_') and h12['H2_hist'].startswith('NT_'))
# B. Radius-two exact group-ring lifts.
p=subprocess.run([sys.executable,str(B/'analyze_gamma3_radius2_v47.py')],cwd=B,capture_output=True,text=True)
ok('g47_radius2_helper_exit_zero',p.returncode==0)
r2=json.loads((B/'gamma3_character_radius2_lifts_v47.json').read_text(encoding='utf-8'))
partial=json.loads((B/'gamma3_partial_resolution_v47.json').read_text(encoding='utf-8'))
ok('g47_ball_163',r2['radius2_ball']['unique_group_elements']==163)
ok('g47_depth_profile',r2['radius2_ball']['depth_profile']=={'0':1,'1':15,'2':147} or r2['radius2_ball']['depth_profile']=={0:1,1:15,2:147})
ok('g47_cols_7009',r2['exact_boundary']['columns']==7009)
ok('g47_rows_13690',r2['exact_boundary']['rows']==13690)
ok('g47_rank_5430',r2['exact_boundary']['rank_Q']==5430)
ok('g47_kernel_1579',r2['exact_boundary']['kernel_dimension']==1579)
ok('g47_residual_image_13',r2['residual_image_dimension_of_exact_radius2_kernel_mod_M43']==13)
ok('g47_relation_e14',r2['residual_annihilator_relation']['primitive_relation']==[0]*13+[1])
ok('g47_13_pure_lifts',len(r2['pure_exact_lifts'])==13)
ok('g47_targets_1_13',r2['lifted_directions']==[f'CHIRES_{i:02d}' for i in range(1,14)])
ok('g47_unlifted_14',r2['unlifted_direction_radius2']=='CHIRES_14' and r2['unlifted_direction_support']=='R3b_123')
ok('g47_all_exact_zero',all(x['exact_Fox_boundary']=='ZERO' for x in r2['pure_exact_lifts']))
ok('g47_all_radius2',all(x['max_group_word_radius']==2 for x in r2['pure_exact_lifts']))
ok('g47_all_amp2',all(x['residual_amplitude']==2 for x in r2['pure_exact_lifts']))
ok('g47_char_after_dim1',partial['character_quotient_after_adjoining_new_cycles']==1)
ok('g47_module_completeness_nt',partial['complete_resolution']=='NT')
# Independent direct check of the CSV exact Fox boundaries in faithful integer matrix group ring.
g3=json.load(open(B/'gamma3_level2_relative_3cells_v36.json',encoding='utf-8'));cells=g3['cells'];gens=g3['generator_order'];ci={c['name']:i for i,c in enumerate(cells)}
STD={}
for i in range(1,4):
 F=sp.eye(3);F[i-1,i-1]=-1;STD[f'F{i}']=F
for i in range(1,4):
 for j in range(1,4):
  if i!=j:
   E=sp.eye(3);E[i-1,j-1]=2;STD[f'E{i}{j}']=E
def mt(M):return tuple(tuple(int(M[i,j]) for j in range(3)) for i in range(3))
def mm(A,C):return tuple(tuple(sum(A[i][k]*C[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def mi(A):return mt(sp.Matrix(A).inv())
I=mt(sp.eye(3));MAT={g:mt(M) for g,M in STD.items()};MATI={g:mi(MAT[g]) for g in MAT}
def grterm(d,g,c):
 if c:
  d[g]=d.get(g,0)+c
  if d[g]==0:del d[g]
def fox(word):
 out={x:{} for x in gens};pref=I
 for t in word:
  inv=t.endswith('^-1');b=t[:-3] if inv else t
  if inv:
   q=MATI[b];grterm(out[b],mm(pref,q),-1);pref=mm(pref,q)
  else:grterm(out[b],pref,1);pref=mm(pref,MAT[b])
 assert pref==I
 return out
FOX=[fox(c['word']) for c in cells]
def parse_group_word(s):
 if s=='1':return I
 out=I
 for t in s.split('*'):
  inv=t.endswith('^-1');b=t[:-3] if inv else t
  out=mm(out,MATI[b] if inv else MAT[b])
 return out
liftrows=list(csv.DictReader(open(B/'gamma3_radius2_pure_lifts_v47.csv',encoding='utf-8')))
by=collections.defaultdict(list)
for r in liftrows:by[r['lift_id']].append(r)
ok('g47_csv_13_lifts',len(by)==13)
for lid,rr in sorted(by.items()):
 bd={}
 for r in rr:
  k=ci[r['relator']];h=parse_group_word(r['group_word']);a=int(r['coefficient'])
  for x in gens:
   for q,c in FOX[k][x].items():
    key=(x,mm(h,q));bd[key]=bd.get(key,0)+a*c
 bd={k:v for k,v in bd.items() if v}
 ok(f'g47_{lid}_fox_zero',not bd)
 ok(f'g47_{lid}_uses_radius2',max(int(r['group_radius']) for r in rr)==2)
# C. Typed guards retained.
r5=json.loads((B/'selling_gamma_R5_adapter_v43.json').read_text(encoding='utf-8'))
ok('g47_R5_strict_refuted_retained','REFUTED-TYPED' in r5['strict_incidence_preserving_16_to_16_adapter'])
sing=json.loads((B/'selling_singer_injection_gate_v45.json').read_text(encoding='utf-8'))
ok('g47_singer_pending',sing['injection_to_J_line'].startswith('NT_'))
ok('g47_lorentz_NT',sing['Lorentz_orientation']=='NT')
registry={'version':'v47','items':[
 {'id':'V47-01','object':'historical abstract HC/HAX carrier','status':'PROVEN_SOURCE_SEMANTIC_C6'},
 {'id':'V47-02','object':'DAX-only center-to-center stitching','status':'REFUTED-TYPED_FOR_CURRENT_TRACE_GRAPH_NO_C6'},
 {'id':'V47-03','object':'direct source-plate HC/HAX raster transcription','status':'NT/PENDING-SOURCE-PLATE-TRANSCRIPTION'},
 {'id':'V47-04','object':'P,Q,R,S route semantics','status':'PROVEN_SOURCE_TEXT / RASTER_BINDING_NT'},
 {'id':'V47-05','object':'89-node historical C2 action','status':'NT_NOT_ACTIVATED_FAIL_CLOSED'},
 {'id':'V47-06','object':'historical C0/C1/C2 totalization','status':'NT'},
 {'id':'V47-07','object':'F_Sell_field -> F29_can','status':'PENDING-MORPHISM'},
 {'id':'V47-08','object':'historical H1/H2','status':'NT_NOT_RECOMPUTED'},
 {'id':'V47-09','object':'13 exact radius-two residual lifts','status':'PROVEN_EXACT_FOX_CYCLES_MINIMAL_RADIUS_2'},
 {'id':'V47-10','object':'CHIRES_14 = R3b_123 radius<=2 lift','status':'REFUTED-TYPED/BOUNDED-NO-GO'},
 {'id':'V47-11','object':'full R-module structure after adding 13 cycles','status':'NT'},
 {'id':'V47-12','object':'complete Gamma3 group-ring resolution','status':'NT; ONE CHARACTER DIRECTION REMAINS AT RADIUS>=3'},
 {'id':'V47-13','object':'strict R5-Gamma incidence adapter','status':'REFUTED-TYPED_RETAINED'},
 {'id':'V47-14','object':'Singer odd injection','status':'SEPARATED/PENDING-MORPHISM'},
 {'id':'V47-15','object':'MXM^-1 covariant use','status':'OPPOSITE_VARIANCE_RETAINED'},
 {'id':'V47-16','object':'roof reversal = transpose','status':'REFUTED-TYPED_RETAINED'},
 {'id':'V47-17','object':'Lorentz orientation','status':'NT'},
 {'id':'V47-18','object':'physical phase','status':'SEPARATED/PENDING-MORPHISM'}]}
(B/'pending_morphism_registry_v47.json').write_text(json.dumps(registry,indent=2,ensure_ascii=False),encoding='utf-8')
ok('g47_registry_18',len(registry['items'])==18)
failed=[k for k,v in checks.items() if not v]
new=len(checks);report={'version':'v47','predecessor_checks':8034,'new_checks':new,'combined_checks':8034+new,'new_passed':new-len(failed),'status':'PASS' if not failed else 'FAIL','failed':failed}
(B/'v47_VERIFICATION_REPORT.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(report,indent=2,ensure_ascii=False))
sys.exit(0 if not failed else 1)
