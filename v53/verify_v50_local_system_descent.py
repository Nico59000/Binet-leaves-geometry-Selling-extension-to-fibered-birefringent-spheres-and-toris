#!/usr/bin/env python3
import json,hashlib,subprocess,sys
from pathlib import Path
B=Path('/mnt/data/v50_work')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
artifact=B/'selling_cover_local_system_descent_v50.json'
builder=B/'build_cover_local_system_descent_v50.py'
before=sha(artifact)
subprocess.run([sys.executable,str(builder)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
after=sha(artifact)
assert before==after
D=json.load(open(artifact))
checks=[]
def ck(name,cond):
 assert cond,name; checks.append(name)
ck('status',D['status']=='PASS-EXACT-V4-LOCAL-SYSTEM-PULLBACK/DESCENT')
for k,up,inv,down in [
 ('Binet',(248,418,150),(62,107,37),(62,107,37)),
 ('Ad',(744,1382,702),(186,348,174),(186,348,174))]:
 S=D['systems'][k]
 ck(k+'_dims_up',tuple(S['dimensions_upstairs'][x] for x in ['C0','C1','C2'])==up)
 ck(k+'_dims_inv',tuple(S['dimensions_V4_invariant'][x] for x in ['C0','C1','C2'])==inv)
 ck(k+'_dims_down',tuple(S['dimensions_downstairs'][x] for x in ['C0','C1','C2'])==down)
 ck(k+'_upper_d1d0',S['upper_d1d0_zero'] is True)
 ck(k+'_d0_descent',S['d0_descent_square'] is True)
 ck(k+'_d1_descent',S['d1_descent_square'] is True)
 ck(k+'_L_full_rank',tuple(S['descent_injection_ranks'][x] for x in ['C0','C1','C2'])==down)
 ck(k+'_objects_equivariance',S['equivariance_object_fiber_checks']==992)
 ck(k+'_arrows_equivariance',S['equivariance_arrow_transport_checks']==1968)
# exact special cells
Bsys=D['systems']['Binet']['C2_transfer_blocks']; Asys=D['systems']['Ad']['C2_transfer_blocks']
bn=[x for x in Bsys if x['deck_monodromy']!='I']; an=[x for x in Asys if x['deck_monodromy']!='I']
ck('Binet_nontrivial_geometric_lifts',len(bn)==6)
ck('Ad_nontrivial_geometric_lifts',len(an)==6)
ck('Binet_TU_transfer_two',all(x['matrix']==[['2']] for x in bn if x['vertex']==48))
ck('Binet_T_U_zero_C2',all(x['cols']==0 for x in bn if x['vertex'] in (59,61)))
ck('Ad_nontrivial_transfer_rank_one',all(x['rank']==1 and x['cols']==1 and x['rows']==3 for x in an))
report={
 'version':'v50','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE','status':'PASS',
 'checks_total':len(checks),'checks_passed':len(checks),'checks':checks,
 'deterministic_replay_byte_identical':True,
 'sha256':{'selling_cover_local_system_descent_v50.json':after,'build_cover_local_system_descent_v50.py':sha(builder)},
 'decision_gate':'Rational cohomology descent may now be opened; integral descent remains separately guarded.'
}
(B/'v50_local_system_descent_verification.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2))
