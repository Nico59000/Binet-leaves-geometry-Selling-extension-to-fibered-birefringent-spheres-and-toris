from pathlib import Path
import json, itertools, collections
import sympy as sp
B=Path('/mnt/data/v45build')
g3=json.loads((B/'gamma3_level2_relative_3cells_v36.json').read_text());cells=g3['cells'];ci={c['name']:i for i,c in enumerate(cells)}

def invtok(t): return t[:-3] if t.endswith('^-1') else t+'^-1'
def invword(w): return tuple(invtok(t) for t in reversed(w))
def reduce_word(w):
 s=[]
 for t in w:
  if s and invtok(t)==s[-1]: s.pop()
  else:s.append(t)
 return tuple(s)
STD={}
for i in range(1,4):
 F=sp.eye(3);F[i-1,i-1]=-1;STD[f'F{i}']=F
for i in range(1,4):
 for j in range(1,4):
  if i!=j:
   E=sp.eye(3);E[i-1,j-1]=2;STD[f'E{i}{j}']=E
mt=lambda M:tuple(tuple(int(M[i,j]) for j in range(3)) for i in range(3))
def mm(A,C):return tuple(tuple(sum(A[i][k]*C[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def mi(A):return mt(sp.Matrix(A).inv())
I=mt(sp.eye(3));MAT={g:mt(M) for g,M in STD.items()};MATI={g:mi(MAT[g]) for g in MAT}
def wordmat(w):
 A=I
 for t in w:
  inv=t.endswith('^-1');b=t[:-3] if inv else t;A=mm(A,MATI[b] if inv else MAT[b])
 return A
inst=[]
for i,c in enumerate(cells):
 w=tuple(c['word'])
 for sign,sw in [(1,w),(-1,invword(w))]:
  for sh in range(len(sw)):
   q=sw[:sh]
   inst.append({'id':len(inst),'idx':i,'name':c['name'],'sign':sign,'shift':sh,'word':sw[sh:]+sw[:sh],'qinv':mi(wordmat(q))})
print('inst',len(inst))
# M43 funcs
fc=json.loads((B/'gamma3_group_ring_4cells_v38.json').read_text())
seedpairs=[]
for cyc in fc['cycles']:
 a=ci[cyc['relator_a']];b=ci[cyc['relator_b']];seedpairs.append((a,b,wordmat(tuple(cells[a]['word']))))
seedcoords={i for a,b,_ in seedpairs for i in (a,b)}
square={}
for k,c in enumerate(cells):
 w=tuple(c['word'])
 if len(w)%2==0 and w[:len(w)//2]==w[len(w)//2:]:
  half=w[:len(w)//2];g=wordmat(half)
  if mm(g,g)==I:square[k]=(half,g)
def grterm(d,g,c):
 if c:
  d[g]=d.get(g,0)+c
  if d[g]==0:del d[g]
def vterm(v,k,g,c):
 d=v.setdefault(k,{});grterm(d,g,c)
 if not d:v.pop(k,None)
def fvec(ids):
 v={}
 for ii in ids:
  x=inst[ii];vterm(v,x['idx'],x['qinv'],x['sign'])
 return v
def right_translate(d,r):
 out={}
 for h,c in d.items():grterm(out,mm(h,r),c)
 return out
def in_ideal(d,g):
 seen=set()
 for h in list(d):
  if h in seen:continue
  hg=mm(h,g);seen|={h,hg}
  if d.get(h,0)+d.get(hg,0)!=0:return False
 return True
def inM(v):
 allowed=seedcoords|set(square)
 if any(k not in allowed and d for k,d in v.items()):return False
 for a,b,r in seedpairs:
  if v.get(b,{})!=right_translate(v.get(a,{}),r):return False
 for k,(half,g) in square.items():
  if not in_ideal(v.get(k,{}),g):return False
 return True
# all ordered pair reduced products, distinct base relators, nonempty
byres=collections.defaultdict(list)
for a in inst:
 for b in inst:
  if a['idx']==b['idx']:continue
  r=reduce_word(a['word']+b['word'])
  if not r:continue
  byres[r].append((a['id'],b['id']))
print('pair residual keys',len(byres))
# match r and inverse r; canonical key min(r,inv) lexical to avoid double
seen=set();ident=[]
for r,A in byres.items():
 ri=invword(r)
 if ri not in byres:continue
 key=min(r,ri)
 if key in seen:continue
 seen.add(key)
 Bp=byres[ri]
 # enumerate but cap? collect all with four distinct relators and no pair inverse
 for a,b in A:
  ia,ib=inst[a],inst[b]
  for c,d in Bp:
   ic,id_=inst[c],inst[d]
   if len({ia['idx'],ib['idx'],ic['idx'],id_['idx']})<4:continue
   ids=(a,b,c,d)
   # product should reduce empty
   w=reduce_word(ia['word']+ib['word']+ic['word']+id_['word'])
   if w:continue
   ident.append(ids)
print('quad identities raw',len(ident))
# canonicalize by cyclic rotation and reversal of occurrence IDs (reverse order with inverse occurrence mapping needed)
# map each occurrence word inverse to exact occurrence id same base
word_to_ids=collections.defaultdict(list)
for x in inst:word_to_ids[(x['idx'],x['word'])].append(x['id'])
def inv_occ(i):
 x=inst[i]; target=invword(x['word']);
 cand=word_to_ids[(x['idx'],target)]
 return min(cand)
def canon(ids):
 seq=list(ids);vars=[]
 for k in range(4):vars.append(tuple(seq[k:]+seq[:k]))
 rev=[inv_occ(i) for i in reversed(seq)]
 for k in range(4):vars.append(tuple(rev[k:]+rev[:k]))
 return min(vars)
classes={}
for ids in ident:classes.setdefault(canon(ids),ids)
print('quad classes',len(classes))
absorbed=0;surv=[]
for can,ids in classes.items():
 v=fvec(ids)
 if inM(v):absorbed+=1
 else:surv.append((can,ids,v))
print('absorbed',absorbed,'survivors',len(surv))
# profile by relator names
print('first survivors')
for can,ids,v in surv[:20]:
 print([inst[i]['name'] for i in ids],[(inst[i]['sign'],inst[i]['shift']) for i in ids], 'coords', [cells[k]['name'] for k in v])
# write summary minimal samples
out={'version':'v45-probe','oriented_instances':len(inst),'pair_residual_keys':len(byres),'four_distinct_free_identity_sequences_raw':len(ident),'dihedral_classes':len(classes),'absorbed_by_M43':absorbed,'survivors_mod_M43':len(surv),'survivor_samples':[{'ids':list(ids),'relators':[inst[i]['name'] for i in ids],'sign_shift':[[inst[i]['sign'],inst[i]['shift']] for i in ids],'support_relators':[cells[k]['name'] for k in v]} for _,ids,v in surv[:50]]}
(B/'gamma3_quadruple_free_identity_probe_v45.json').write_text(json.dumps(out,indent=2))
