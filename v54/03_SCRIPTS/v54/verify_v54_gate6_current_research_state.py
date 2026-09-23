#!/usr/bin/env python3
from pathlib import Path
import json,subprocess,sys
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
C=json.loads((ROOT/"v54_GATE6_CURRENT_RESEARCH_STATE.json").read_text(encoding="utf-8"))
checks=[]
def ok(n,v):
    checks.append((n,bool(v)))
    if not v: raise AssertionError(n)

# Replay both current component verifiers.
p1=subprocess.run([sys.executable,str(BASE/"v54_gate6_consolidated"/"verify_v54_gate6_consolidated.py")],
                  cwd=BASE/"v54_gate6_consolidated",capture_output=True,text=True)
p2=subprocess.run([sys.executable,str(BASE/"v54_gate6_noether_laneA"/"verify_v54_gate6_pascal_noether_typed_bridge.py")],
                  cwd=BASE/"v54_gate6_noether_laneA",capture_output=True,text=True)
ok("consolidated_replay_exit0",p1.returncode==0)
ok("laneA_replay_exit0",p2.returncode==0)
R1=json.loads((BASE/"v54_gate6_consolidated"/"v54_GATE6_CONSOLIDATED_VERIFICATION_REPORT.json").read_text())
R2=json.loads((BASE/"v54_gate6_noether_laneA"/"v54_GATE6_PASCAL_NOETHER_TYPED_BRIDGE_VERIFICATION.json").read_text())
ok("consolidated_59",R1["status"]=="PASS" and R1["combined_checkpoint_checks"]==59)
ok("laneA_102",R2["status"]=="PASS" and R2["checks"]==102 and R2["passed"]==102)
ok("same_gate",C["policy"]=="SAME-GATE-CONTINUATION/NO-SUBGATE")
ok("v53_sealed",C["sealed_predecessor"]=="v53 CLOSED/APPEND-ONLY")
ok("v54_nonpublishable",C["publication_state"]=="RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE")
ok("fig3_2_91",C["Fig3"]["direct_raster_bound"]=="2/91")
ok("fig2_unique_binding_no_go",C["Fig2"]["point_level_binding_status"].startswith("REFUTED-TYPED"))
ok("laneA_12",C["Pascal_Noether_lane_A"]["type_preserving_dictionary_candidates"]==12)
ok("laneA_termwise_NT",C["Pascal_Noether_lane_A"]["termwise_historical_dictionary"]=="OPEN/NT")
ok("p46_coeff_NT",C["Pascal_Noether_lane_A"]["exact_wrong_p46_coefficient"]=="OPEN/NT")
ok("r2_periodic_refuted",C["R2PURE_13"]["periodic_tail_only_explanation"]=="REFUTED-TYPED")
ok("r4_closed",C["retained_closed_guards"]["r4"]=="NOT-OPENED")
ok("B6_closed",C["retained_closed_guards"]["B6"]=="NOT-OPENED")
ok("gate_open",C["gate6_status"]=="KEEP-OPEN")

rep={
 "version":"v54","gate":"Gate 6","status":"PASS",
 "component_check_events":{"consolidated":59,"pascal_noether_bridge":102,"total":161},
 "orchestrator_checks":len(checks),"orchestrator_passed":sum(v for _,v in checks),
 "failed":[n for n,v in checks if not v],
 "publication_state":C["publication_state"],"gate_closed":False
}
(ROOT/"v54_GATE6_CURRENT_RESEARCH_VERIFICATION.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n")
print(json.dumps(rep,indent=2))
