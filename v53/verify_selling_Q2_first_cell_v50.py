#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,subprocess,sys
H=Path(__file__).resolve().parent
files=[H/'selling_Q2_first_cell_v50.json',H/'selling_Q2_first_cell_boundary_v50.csv']
subprocess.run([sys.executable,str(H/'selling_Q2_first_cell_v50.py')],check=True,capture_output=True,text=True)
h1={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
subprocess.run([sys.executable,str(H/'selling_Q2_first_cell_v50.py')],check=True,capture_output=True,text=True)
h2={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
assert h1==h2
j=json.loads((H/'selling_Q2_first_cell_v50.json').read_text())
assert j['status']=='PASS-LOCAL-CELL'
assert all(j['checks'].values())
assert j['decision']['first_source_constructed_2cell']=='PROVEN-LOCAL'
assert j['decision']['global_C0_C1_C2'].startswith('NT/')
print(json.dumps({'status':'PASS','checks':len(j['checks']),'deterministic_replay':True,'sha256':h2},indent=2))
