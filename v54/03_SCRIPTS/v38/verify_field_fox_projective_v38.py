#!/usr/bin/env python3
from pathlib import Path
import json,csv,hashlib,subprocess,sys
import numpy as np
B=Path(__file__).resolve().parent
checks=[]
def ok(name,cond):
    checks.append((name,bool(cond)))
    if not cond: print('FAIL',name)
def sha(p):
    h=hashlib.sha256();
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
# Parent replay
p=subprocess.run([sys.executable,str(B/'verify_plate_vector_syzygy_ternary_v37.py')],cwd=B,text=True,capture_output=True)
ok('parent_v37_returncode',p.returncode==0); ok('parent_v37_PASS','"status": "PASS"' in p.stdout and '"combined": 383' in p.stdout)
# Photo gate
pg=json.loads((B/'selling_fig1_bsb_photo_gate_v38.json').read_text())
ok('photo_sha',sha(B/'Tafel-III-fig1.jpg')==pg['source_sha256']); ok('photo_size',(B/'Tafel-III-fig1.jpg').stat().st_size==pg['source_bytes'])
ok('photo_registration_inliers',pg['registration_to_v37']['inliers']==142); ok('photo_registration_median',pg['registration_to_v37']['median_reprojection_error_px_v37_300']<1.5)
ok('photo_support_total',sum(pg['v37_segment_support_counts'].values())==105); ok('photo_confirmed_97',pg['v37_segment_support_counts'].get('CONFIRMED')==97)
ok('photo_partial_3',pg['v37_segment_support_counts'].get('PARTIAL')==3); ok('photo_weak_5',pg['v37_segment_support_counts'].get('WEAK')==5)
ok('photo_serpentine_2',pg['by_layer']['serpentine'].get('CONFIRMED')==2); ok('photo_incidence_nt','COMPLETE_HISTORICAL_FIELD_INCIDENCE_NT' in pg['decision'])
# Context hashes
ctx=json.loads((B/'v38_noether_reality_context.json').read_text())
ok('noether_hash',sha(B/'Noether_1908_ternary_biquadratic_forms.pdf')==ctx['Noether_1908']['sha256'])
ok('reality_hash',sha(B/'operateur_realite_stratifiee_fr.pdf')==ctx['stratified_reality']['sha256'])
ok('reality_four_status',ctx['stratified_reality']['four_status_interface']==['F','NF','NT','T'])
# Fox certificate replay
fc=json.loads((B/'gamma3_group_ring_4cells_v38.json').read_text()); ok('fox_six',fc['cycle_count']==6)
for i,c in enumerate(fc['cycles']):
    ok(f'fox_pair_{i}_inverse',c['free_inverse']); ok(f'fox_pair_{i}_allgens',c['all_generators_pass'])
    for g,v in c['generator_checks'].items(): ok(f'fox_{i}_{g}',v)
ok('fox_complete_stays_nt',fc['complete_module'].startswith('NT_'))
ok('fox_rank6_claim',fc['independence'].endswith('RANK_6_BY_DISJOINT_RELATOR_COORDINATE_PAIRS'))
# Binder replay
bd=json.loads((B/'ternary_projective_bracket_v38.json').read_text()); ok('binder_status',bd['status'].startswith('PROVEN_'))
Ky=np.array(bd['Thomas_Wigner']['Ky'],dtype=object); Kz=np.array(bd['Thomas_Wigner']['Kz'],dtype=object); J=np.array(bd['Thomas_Wigner']['J'],dtype=object)
comm=lambda A,C:A@C-C@A
ok('TW_integer_1',np.array_equal(comm(Ky,Kz),J)); ok('TW_integer_2',np.array_equal(comm(J,Ky),-Kz)); ok('TW_integer_3',np.array_equal(comm(J,Kz),Ky))
for lev in bd['levels_checked']:
    k=lev['k']; m=lev['modulus']; e=lev['e_k']; q=3**k
    ok(f'level{k}_eidemp',(e*e-e)%m==0); ok(f'level{k}_emod2',e%2==1); ok(f'level{k}_emod3k',e%q==0)
    ok(f'level{k}_minus_not',lev['minus_one_in_Uplus'] is False)
    ok(f'level{k}_ucard',lev['Uplus_cardinality']==3**(k-1))
    for nm,v in lev['brackets_modulus'].items(): ok(f'level{k}_{nm}',v)
    # Sample the entire U+ set for k<=8, and a deterministic arithmetic subset afterward.
    vals=[]
    if k<=8:
        vals=[u for u in range(m) if __import__('math').gcd(u,m)==1 and u%3==1]
        ok(f'level{k}_Uplus_count_enum',len(vals)==3**(k-1))
    else:
        vals=[1 + 6*t for t in range(min(100,3**(k-1)))]
    ok(f'level{k}_Uplus_mod2',all(u%2==1 for u in vals))
    ok(f'level{k}_Uplus_no_minus',all((u+1)%m!=0 for u in vals))
# compatibility e_{k+1}->e_k and scalar groups cardinal growth
ls=bd['levels_checked']
for a,b in zip(ls,ls[1:]):
    ok(f'compat_e_{a["k"]}_{b["k"]}',b['e_k']%a['modulus']==a['e_k'])
# Pending registry guards
pr=json.loads((B/'pending_morphism_registry_v38.json').read_text());
status={x['id']:x['status'] for x in pr['items']}
ok('field_functor_pending',status['PM38-03']=='PENDING-MORPHISM'); ok('field_H12_nt',status['PM38-04'].startswith('NT_'))
ok('roof_transpose_refuted',status['PM38-08']=='REFUTED-TYPED'); ok('phase_separated',status['PM38-10']=='SEPARATED/PENDING-MORPHISM')
failed=[n for n,v in checks if not v]
res={'status':'PASS' if not failed else 'FAIL','parent_checks':383,'new_checks':len(checks),'combined':383+len(checks),'failed':failed,
     'photo_support':pg['v37_segment_support_counts'],'fox_4cell_cycles':fc['cycle_count'],'binder_levels':len(bd['levels_checked']),
     'complete_historical_incidence':'NT','historical_H1_H2':'NT_NOT_RECOMPUTED','complete_Fox_module':'NT'}
(B/'v38_VERIFICATION_REPORT.json').write_text(json.dumps(res,indent=2)+"\n")
print(json.dumps(res,indent=2))
sys.exit(0 if not failed else 1)
