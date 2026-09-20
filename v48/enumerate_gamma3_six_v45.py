from pathlib import Path
import json,collections,itertools
import sympy as sp
B=Path('/mnt/data/v45build');cells=json.load(open(B/'gamma3_level2_relative_3cells_v36.json'))['cells'];ci={c['name']:i for i,c in enumerate(cells)}
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
STD={}
for i in range(1,4):F=sp.eye(3);F[i-1,i-1]=-1;STD[f'F{i}']=F
for i in range(1,4):
 for j in range(1,4):
  if i!=j:E=sp.eye(3);E[i-1,j-1]=2;STD[f'E{i}{j}']=E
mt=lambda M:tuple(tuple(int(M[i,j]) for j in range(3)) for i in range(3))
def mm(A,C):return tuple(tuple(sum(A[i][k]*C[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def mi(A):return mt(sp.Matrix(A).inv())
I=mt(sp.eye(3));MAT={g:mt(M) for g,M in STD.items()};MATI={g:mi(MAT[g]) for g in MAT}
def wordmat(w):
 A=I
 for t in w:
  iv=t.endswith('^-1');b=t[:-3] if iv else t;A=mm(A,MATI[b] if iv else MAT[b])
 return A
inst=[]
for i,c in enumerate(cells):
 w=tuple(c['word'])
 for sign,sw in [(1,w),(-1,invword(w))]:
  for sh in range(len(sw)):
   q=sw[:sh];inst.append({'id':len(inst),'idx':i,'name':c['name'],'sign':sign,'shift':sh,'word':sw[sh:]+sw[:sh],'qinv':mi(wordmat(q))})
# M43
fc=json.load(open(B/'gamma3_group_ring_4cells_v38.json'));seedpairs=[]
for cyc in fc['cycles']:
 a=ci[cyc['relator_a']];b=ci[cyc['relator_b']];seedpairs.append((a,b,wordmat(tuple(cells[a]['word']))))
seedcoords={i for a,b,_ in seedpairs for i in (a,b)};square={}
for k,c in enumerate(cells):
 w=tuple(c['word'])
 if len(w)%2==0 and w[:len(w)//2]==w[len(w)//2:]:
  half=w[:len(w)//2];g=wordmat(half)
  if mm(g,g)==I:square[k]=(half,g)
def grterm(d,g,c):
 if c:d[g]=d.get(g,0)+c; d.pop(g,None) if d.get(g)==0 else None
def vterm(v,k,g,c):
 d=v.setdefault(k,{});grterm(d,g,c)
 if not d:v.pop(k,None)
def fvec(ids):
 v={}
 for i in ids:
  x=inst[i];vterm(v,x['idx'],x['qinv'],x['sign'])
 return v
def right_translate(d,r):
 out={}
 for h,c in d.items():grterm(out,mm(h,r),c)
 return out
def inideal(d,g):
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
  if not inideal(v.get(k,{}),g):return False
 return True
pairs=[]
for a in inst:
 for b in inst:
  if a['idx']==b['idx']:continue
  if bc(a['word'],b['word'])<1:continue
  pairs.append((a['id'],b['id'],reduce_word(a['word']+b['word'])))
T=collections.defaultdict(list)
for ai,bi,r in pairs:
 a,b=inst[ai],inst[bi]
 for c in inst:
  if c['idx'] in (a['idx'],b['idx']):continue
  if bc(r,c['word'])<1:continue
  T[reduce_word(r+c['word'])].append((ai,bi,c['id']))
ids6=[];seenkey=set()
for r,A in T.items():
 ri=invword(r)
 if ri not in T:continue
 key=min(r,ri)
 if key in seenkey:continue
 seenkey.add(key)
 for a in A:
  ia={inst[x]['idx'] for x in a}
  for b in T[ri]:
   if len(ia|{inst[x]['idx'] for x in b})<6:continue
   ids=a+b
   if reduce_word(tuple(t for i in ids for t in inst[i]['word'])):continue
   ids6.append(ids)
print('raw',len(ids6))
# canonical under cyclic rotation and reversal occurrences
by_word=collections.defaultdict(list)
for x in inst:by_word[(x['idx'],x['word'])].append(x['id'])
def invocc(i):return min(by_word[(inst[i]['idx'],invword(inst[i]['word']))])
def canon(seq):
 seq=list(seq);vs=[];n=len(seq)
 for k in range(n):vs.append(tuple(seq[k:]+seq[:k]))
 rev=[invocc(i) for i in reversed(seq)]
 for k in range(n):vs.append(tuple(rev[k:]+rev[:k]))
 return min(vs)
cl={}
for ids in ids6:cl.setdefault(canon(ids),ids)
absorb=0;surv=[]
for ca,ids in cl.items():
 v=fvec(ids)
 if inM(v):absorb+=1
 else:surv.append((ids,v))
print('classes',len(cl),'absorbed',absorb,'surv',len(surv))
out={'version':'v45-probe','scope':'six-relator free identities assembled as inverse residual pairs of v44 sequential positive-cancellation triples, six distinct base relators','raw_sequences':len(ids6),'dihedral_classes':len(cl),'absorbed_by_M43':absorb,'survivors_mod_M43':len(surv),'decision':'NO_NEW_DIRECTION' if not surv else 'NEW_DIRECTIONS_FOUND'}
(B/'gamma3_six_relator_probe_v45.json').write_text(json.dumps(out,indent=2))
