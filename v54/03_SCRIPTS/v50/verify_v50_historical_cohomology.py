#!/usr/bin/env python3
import json, hashlib
from pathlib import Path
B=Path('/mnt/data/v50_work')
Q=json.load(open(B/'selling_historical_cohomology_Q_v50.json'))
I=json.load(open(B/'selling_integral_orbifold_cochains_v50.json'))
Z=json.load(open(B/'selling_historical_cohomology_Z_2primary_audit_v50.json'))
checks={}
checks['Q_status']=Q['status']=='PASS-EXACT-DOUBLE-RATIONAL-HISTORICAL-COHOMOLOGY'
checks['Q_Binet_H1_23']=Q['systems']['Binet']['downstairs']['H_dimensions']['H1']==23
checks['Q_Binet_H2_15']=Q['systems']['Binet']['downstairs']['H_dimensions']['H2']==15
checks['Q_Ad_H1_3']=Q['systems']['Ad']['downstairs']['H_dimensions']['H1']==3
checks['Q_Ad_H2_15']=Q['systems']['Ad']['downstairs']['H_dimensions']['H2']==15
for s in ['Binet','Ad']:
 checks[f'{s}_Q_class_H1_match']=Q['systems'][s]['class_comparison']['H1_representatives_match_in_invariant_coordinates']
 checks[f'{s}_Q_class_H2_match']=Q['systems'][s]['class_comparison']['H2_representatives_match_in_invariant_coordinates']
 checks[f'{s}_integral_d1d0']=I['systems'][s]['d1d0_zero']
 checks[f'{s}_H1_Z2']=Z['systems'][s]['groups']['H1']['torsion_invariant_factors']==[2]
 checks[f'{s}_H2_torsionfree']=Z['systems'][s]['groups']['H2']['torsion_invariant_factors']==[]
 checks[f'{s}_mod2_jump1']=Z['systems'][s]['mod_p_audit']['2']['H1_dim']==Z['systems'][s]['groups']['H1']['free_rank']+1
 checks[f'{s}_odd3_nojump']=Z['systems'][s]['mod_p_audit']['3']['H1_dim']==Z['systems'][s]['groups']['H1']['free_rank']
checks['Binet_H1_Z_rank23']=Z['systems']['Binet']['groups']['H1']['free_rank']==23
checks['Ad_H1_Z_rank3']=Z['systems']['Ad']['groups']['H1']['free_rank']==3
checks['Binet_H2_Z_rank15']=Z['systems']['Binet']['groups']['H2']['free_rank']==15
checks['Ad_H2_Z_rank15']=Z['systems']['Ad']['groups']['H2']['free_rank']==15
checks['Binet_TU_transfer_index2']=Z['transfer_indices']['Binet']['TU']['index']==2
checks['Ad_T_transfer_index1']=Z['transfer_indices']['Ad']['T']['index']==1
checks['Ad_U_transfer_index1']=Z['transfer_indices']['Ad']['U']['index']==1
checks['Ad_TU_transfer_index1']=Z['transfer_indices']['Ad']['TU']['index']==1
checks['Binet_up_inv_H2_extra_Z2']=Z['integral_descent']['upstairs_V4_invariant_cohomology']['Binet']['H2']['torsion_invariant_factors']==[2]
checks['Ad_integral_descent_identical']=Z['integral_descent']['upstairs_V4_invariant_cohomology']['Ad']['H2']['torsion_invariant_factors']==[]
assert all(checks.values()),[k for k,v in checks.items() if not v]
files=['selling_historical_cohomology_Q_v50.json','selling_integral_orbifold_cochains_v50.json','selling_historical_cohomology_Z_2primary_audit_v50.json']
report={'version':'v50','status':'PASS','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE','checks_passed':sum(checks.values()),'checks_total':len(checks),'checks':checks,'sha256':{f:hashlib.sha256((B/f).read_bytes()).hexdigest() for f in files},'decision':'Rational historical cohomology and integral downstairs orbifold cohomology are calculated exactly. Integral V4-invariant upstairs comparison agrees for Ad and differs for Binet H2 by one Z/2 measured by the TU transfer index 2.'}
(B/'v50_historical_cohomology_verification.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(f"PASS {len(checks)}/{len(checks)}")
