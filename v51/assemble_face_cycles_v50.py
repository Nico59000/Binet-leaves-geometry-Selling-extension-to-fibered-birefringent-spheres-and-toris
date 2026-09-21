#!/usr/bin/env python3
import json,math
from pathlib import Path
HERE=Path('/mnt/data/v50_work')
parts=[]
for name in ['selling_face_segment_candidates_0_25_v50.json','selling_face_segment_candidates_25_45_v50.json','selling_face_segment_candidates_45_50_v50.json','selling_face_segment_candidates_tail_v50.json']:
 d=json.load(open(HERE/name));parts.extend(d['faces'])
by={f['field']:f for f in parts}
assert set(by)==set(range(62))

def order_cycle(face,allowed_segments=None):
 segs=face['segments'] if allowed_segments is None else allowed_segments
 # adjacency vertices to segment ids
 vad={}
 for si,s in enumerate(segs):
  for v in s['v']:vad.setdefault(v,[]).append(si)
 usedv=sorted(vad)
 assert all(len(vad[v])==2 for v in usedv), (face['field'],{v:vad[v] for v in usedv})
 start=min(usedv,key=lambda v:(face['vertices'][v]['q'],face['vertices'][v]['r']))
 # two possible first neighbors, choose one then orient by area later
 s0=vad[start][0];curv=start;curs=s0;seqv=[start];seqs=[]
 while True:
  s=segs[curs];seqs.append(curs)
  nxt=s['v'][1] if s['v'][0]==curv else s['v'][0]
  if nxt==start:break
  seqv.append(nxt)
  opts=[x for x in vad[nxt] if x!=curs];assert len(opts)==1
  curv=nxt;curs=opts[0]
  if len(seqs)>len(segs)+2:raise RuntimeError('cycle')
 assert len(seqs)==len(segs)
 # signed polygon area using vertices; reverse if clockwise
 pts=[(face['vertices'][v]['q'],face['vertices'][v]['r']) for v in seqv]
 area=sum(pts[i][0]*pts[(i+1)%len(pts)][1]-pts[(i+1)%len(pts)][0]*pts[i][1] for i in range(len(pts)))/2
 if area<0:
  seqv=[seqv[0]]+list(reversed(seqv[1:]))
  # recompute segment sequence from consecutive vertices
  seqs=[]
  for a,b in zip(seqv,seqv[1:]+seqv[:1]):
   mm=[si for si,s in enumerate(segs) if set(s['v'])=={a,b}];assert len(mm)==1;seqs.append(mm[0])
  area=-area
 return seqv,seqs,area

out=[]
for i in range(62):
 f=by[i]
 if i!=55:
  vv,ss,area=order_cycle(f)
  cycles=[{'component_occurrence':0,'vertices':vv,'segments':[f['segments'][k] for k in ss],'area_approx':area}]
 else:
  # two exact graph cycles sharing vertex 5; choose lobe A as quotient representative.
  # identify by segment sets through vertices {0,1,4,5} and {2,3,5,6}
  lobes=[]
  for Vset in [{0,1,4,5},{2,3,5,6}]:
   seg=[s for s in f['segments'] if set(s['v']).issubset(Vset)]
   ff=dict(f);ff['segments']=seg
   vv,ss,area=order_cycle(ff)
   lobes.append({'vertices':vv,'segments':[seg[k] for k in ss],'area_approx':area})
  cycles=lobes
 out.append({'field':i,'vertices':f['vertices'],'components':f['components'],'cycles':cycles})
res={'version':'v50','status':'PASS-NUMERIC-CYCLE-STRUCTURE','form_classes':62,'prequotient_field_occurrences':63,'quotient_face_orbits_expected':62,'faces':out,
     'field55_note':'two quadrilateral components sharing one closure vertex; component orbit exchanged by certified U repetition; first lobe chosen as quotient representative pending exact transport verification'}
(HERE/'selling_face_cycles_v50.json').write_text(json.dumps(res,indent=2),encoding='utf-8')
print('faces',len(out),'pre cycles',sum(len(x['cycles']) for x in out),'quotient reps',61+1)
print('cycle length hist')
from collections import Counter
print(Counter(len(c['segments']) for x in out for c in (x['cycles'] if x['field']!=55 else [x['cycles'][0]])))
