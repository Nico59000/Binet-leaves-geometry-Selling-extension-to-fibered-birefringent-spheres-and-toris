from pathlib import Path
import json,collections,math,itertools
B=Path('/mnt/data/v45build'); g3=json.load(open(B/'gamma3_level2_relative_3cells_v36.json'));cells=g3['cells']
def invtok(t):return t[:-3] if t.endswith('^-1') else t+'^-1'
def invword(w):return tuple(invtok(t) for t in reversed(w))
def reduce_word(w):
 s=[]
 for t in w:
  if s and invtok(t)==s[-1]:s.pop()
  else:s.append(t)
 return tuple(s)
def bc(a,b):
 k=0
 for j in range(1,min(len(a),len(b))+1):
  if tuple(a[-j:])==invword(b[:j]):k=j
  else:break
 return k
inst=[]
for i,c in enumerate(cells):
 w=tuple(c['word'])
 for sign,sw in [(1,w),(-1,invword(w))]:
  for sh in range(len(sw)):inst.append({'id':len(inst),'idx':i,'word':sw[sh:]+sw[:sh]})
pairs=[]; pair_by_res=collections.defaultdict(list)
for a in inst:
 for b in inst:
  if a['idx']==b['idx']:continue
  k=bc(a['word'],b['word'])
  if k<1:continue
  r=reduce_word(a['word']+b['word']); pairs.append((a['id'],b['id'],r));pair_by_res[r].append((a['id'],b['id']))
print('pairs',len(pairs),'keys',len(pair_by_res))
triple_by_res=collections.defaultdict(list);nt=0
for ai,bi,r in pairs:
 a,b=inst[ai],inst[bi]
 for c in inst:
  if c['idx'] in (a['idx'],b['idx']):continue
  if bc(r,c['word'])<1:continue
  rr=reduce_word(r+c['word']);triple_by_res[rr].append((ai,bi,c['id']));nt+=1
print('triples',nt,'keys',len(triple_by_res))
raw=0; keys=0; combos=[]
for r,T in triple_by_res.items():
 P=pair_by_res.get(invword(r))
 if not P:continue
 keys+=1
 c=0
 for t in T:
  ti={inst[x]['idx'] for x in t}
  for p in P:
   if len(ti|{inst[p[0]]['idx'],inst[p[1]]['idx']})==5:c+=1
 raw+=c
 if c:combos.append((c,r,len(T),len(P)))
print('identity keys',keys,'raw distinctbase5',raw,'top',sorted(combos,reverse=True)[:20])
