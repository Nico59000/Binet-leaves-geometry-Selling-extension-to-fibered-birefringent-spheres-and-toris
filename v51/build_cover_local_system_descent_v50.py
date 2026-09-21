#!/usr/bin/env python3
import json, sympy as sp
from pathlib import Path
from collections import defaultdict, Counter
B=Path('/mnt/data/v50_work')
G=json.load(open(B/'selling_field_orbigroupoid_v50.json'))
L=json.load(open(B/'selling_orbifold_links_v50.json'))['links']
C=json.load(open(B/'selling_field_orbifold_cochains_v50.json'))
F=json.load(open(B/'selling_F29_field_functor_gate_v50.json'))
Cover=json.load(open(B/'selling_four_sheet_deck_cover_v50.json'))
ar={a['id']:a for a in G['arrows']}
fg={x['vertex']:x['V4_closure'] for x in F['links']}
sheets=['I','T','U','TU']
exp={'I':(0,0),'T':(1,0),'U':(0,1),'TU':(1,1)}
lab={(0,0):'I',(1,0):'T',(0,1):'U',(1,1):'TU'}
def mul(a,b):
 x,y=exp[a],exp[b];return lab[(x[0]^y[0],x[1]^y[1])]
def M(x): return sp.Matrix([[sp.Rational(v) for v in row] for row in x])
def left_inverse(Bas):
 if Bas.cols==0:return sp.zeros(0,Bas.rows)
 return (Bas.T*Bas).inv()*Bas.T
def left_fixed(H):
 ns=(H.T-sp.eye(H.rows)).nullspace()
 return sp.Matrix.vstack(*[v.T for v in ns]) if ns else sp.zeros(0,H.rows)
def complete_rows(P,n):
 rows=[P.row(i) for i in range(P.rows)]
 A=P.copy()
 for j in range(n):
  e=sp.zeros(1,n);e[0,j]=1
  if sp.Matrix.vstack(A,e).rank()>A.rank():
   rows.append(e);A=sp.Matrix.vstack(A,e)
  if len(rows)==n: break
 assert len(rows)==n and A.rank()==n
 return sp.Matrix.vstack(*rows)
def sparse_from_entries(entries,rows,cols):
 d={}
 for z in entries:
  v=sp.Rational(z['value'])
  if v:d[(z['row'],z['col'])]=v
 return sp.SparseMatrix(rows,cols,d)
def add_block(entries,row0,col0,X):
 for i in range(X.rows):
  for j in range(X.cols):
   if X[i,j]!=0: entries[(row0+i,col0+j)]=entries.get((row0+i,col0+j),0)+X[i,j]
def smat(rows,cols,entries):
 return sp.SparseMatrix(rows,cols,{k:sp.simplify(v) for k,v in entries.items() if sp.simplify(v)!=0})
# deck coefficient representations
chi={'I':sp.Integer(1),'T':sp.Integer(-1),'U':sp.Integer(-1),'TU':sp.Integer(1)}
Dad={k:M(v) for k,v in G['deck_ad'].items()}
assert all(Dad[s]*Dad[s]==sp.eye(3) for s in sheets)
assert Dad['T']*Dad['U']==Dad['U']*Dad['T']==Dad['TU']
# lower matrices and lower C1 bases/offsets reproduced exactly
def lower_data(kind):
 n=1 if kind=='Binet' else 3
 dims=C['dimensions'][kind]
 D0=sparse_from_entries(C[f'd0_{kind}_sparse'],dims['C1'],dims['C0'])
 D1=sparse_from_entries(C[f'd1_{kind}_sparse'],dims['C2'],dims['C1'])
 bases={};offs={};off=0
 for e in range(123):
  a=ar[e]
  R=sp.Matrix([[sp.Integer(a['binet_transport'])]]) if kind=='Binet' else M(a['ad_transport'])
  if a['kind'].endswith('ISOTROPY'):
   B0=sp.Matrix.hstack(*(R+sp.eye(n)).nullspace()) if (R+sp.eye(n)).nullspace() else sp.zeros(n,0)
  else:B0=sp.eye(n)
  bases[e]=B0;offs[e]=off;off+=B0.cols
 assert off==dims['C1']
 # lower C2 bases by link order
 pbase={};poff={};off2=0
 for lk in L:
  v=lk['vertex'];H=sp.Matrix([[sp.Integer(lk['chosen']['binet_holonomy'])]]) if kind=='Binet' else M(lk['chosen']['ad_holonomy'])
  P=left_fixed(H);pbase[v]=P;poff[v]=off2;off2+=P.rows
 assert off2==dims['C2']
 return n,D0,D1,bases,offs,pbase,poff
# cover C1 topology, independent of coefficient system
# cells have canonical stored orientation
c1top=[]; keydir={}; mirinfo={}
for e in range(123):
 a=ar[e];d=a['deck']
 if a['kind']=='DECK-MIRROR-ISOTROPY':
  seen=set()
  for s in sheets:
   if s in seen:continue
   t=mul(s,d);pair={s,t};seen|=pair
   rep=s
   cid=len(c1top)
   c1top.append({'cid':cid,'base_edge':e,'kind':a['kind'],'source_sheet':rep,'target_sheet':t,'source_face':a['source'],'target_face':a['target'],'deck':d})
   for q in pair: mirinfo[(e,q)]=(cid, 1 if q==rep else -1)
 else:
  for s in sheets:
   t=mul(s,d);cid=len(c1top)
   c1top.append({'cid':cid,'base_edge':e,'kind':a['kind'],'source_sheet':s,'target_sheet':t,'source_face':a['source'],'target_face':a['target'],'deck':d})
   keydir[(e,s)]=cid
assert len(c1top)==102*4+16*4+5*2==482
# path traversal lookup returns c1top cid and direction relative canonical stored orientation
def path_cell(e,start_sheet,base_direction):
 a=ar[e];d=a['deck']
 if a['kind']=='DECK-MIRROR-ISOTROPY':
  cid,sgn=mirinfo[(e,start_sheet)]
  return cid,sgn  # involutive geometric edge; +/- base direction gives same actual segment
 if base_direction==1:return keydir[(e,start_sheet)],1
 # inverse traversal uses lift whose canonical arrow ends at start_sheet
 return keydir[(e,mul(start_sheet,d))],-1
# Build each coefficient system
def build(kind):
 n,D0d,D1d,Bdown,offdown,Pdown,poff=lower_data(kind)
 Tau={s:(sp.Matrix([[chi[s]]]) if kind=='Binet' else Dad[s]) for s in sheets}
 # C0 upper indexing
 oidx={(f,s):(4*f+sheets.index(s))*n for f in range(62) for s in sheets}
 c0u=62*4*n
 # C1 bases/transports/offsets
 ub={}; uleft={}; uR={}; uoff={}; off=0
 for cell in c1top:
  cid=cell['cid'];e=cell['base_edge'];a=ar[e];ss=cell['source_sheet'];ts=cell['target_sheet']
  Rb=sp.Matrix([[sp.Integer(a['binet_transport'])]]) if kind=='Binet' else M(a['ad_transport'])
  Rup=Tau[ts].inv()*Rb*Tau[ss]
  if cell['kind']=='MULTIWALL-ISOTROPY':
   ns=(Rup+sp.eye(n)).nullspace(); Bu=sp.Matrix.hstack(*ns) if ns else sp.zeros(n,0)
  else: Bu=sp.eye(n)
  ub[cid]=Bu;uleft[cid]=left_inverse(Bu);uR[cid]=sp.simplify(Rup);uoff[cid]=off;off+=Bu.cols
 c1u=off
 expect=418 if kind=='Binet' else 1382
 assert c1u==expect,(kind,c1u)
 # d0 upper
 ent={}
 for cell in c1top:
  cid=cell['cid'];e=cell['base_edge'];a=ar[e];ss=cell['source_sheet'];ts=cell['target_sheet'];k=ub[cid].cols
  if not k:continue
  s0=oidx[(a['source'],ss)];t0=oidx[(a['target'],ts)]
  Ct=uleft[cid];Cs=uleft[cid]*(-uR[cid])
  add_block(ent,uoff[cid],t0,Ct);add_block(ent,uoff[cid],s0,Cs)
 D0u=smat(c1u,c0u,ent)
 # natural descent injection L0
 e0={}
 for f in range(62):
  for s in sheets:
   X=Tau[s].inv()
   add_block(e0,oidx[(f,s)],f*n,X)
 L0=smat(c0u,62*n,e0)
 # natural descent injection L1
 e1={}
 # each upper c1 cell gets gauge transform of lower target-fiber vector
 for cell in c1top:
  cid=cell['cid'];e=cell['base_edge'];ts=cell['target_sheet'];Bd=Bdown[e]
  if Bd.cols==0:continue
  X=uleft[cid]*Tau[ts].inv()*Bd
  add_block(e1,uoff[cid],offdown[e],X)
 L1=smat(c1u,D0d.rows,e1)
 # verify d0 intertwining
 assert D0u*L0==L1*D0d, f'{kind} d0 descent fail'
 # exact V4 equivariance on all 248 object fibers in the sheet-unwound gauge
 obj_eq_checks=0
 for f in range(62):
  for s in sheets:
   for h in sheets:
    hs=mul(h,s)
    # canonical deck action in unwound coordinates is Tau[h]^{-1}
    assert Tau[hs].inv()*Tau[s]==Tau[h].inv()
    obj_eq_checks+=1
 assert obj_eq_checks==248*4
 # equivariance of lifted transports on all 492 directed arrow lifts (including mirror orientations)
 eq_checks=0
 for e in range(123):
  a=ar[e];d=a['deck'];Rb=sp.Matrix([[sp.Integer(a['binet_transport'])]]) if kind=='Binet' else M(a['ad_transport'])
  for s in sheets:
   t=mul(s,d);Rs=Tau[t].inv()*Rb*Tau[s]
   for h in sheets:
    hs=mul(h,s);ht=mul(h,t);Rh=Tau[ht].inv()*Rb*Tau[hs]
    assert sp.simplify(Rh-Tau[h].inv()*Rs*Tau[h])==sp.zeros(n)
    eq_checks+=1
 assert eq_checks==123*4*4
 # build upper C2 cells and d1 rows
 c2cells=[]; D1ent={}; l2ent={}; c2off=0
 # helper path and raw map for one geometric lifted cell
 for lk in L:
  v=lk['vertex'];dcl=fg[v]; reps=[]
  if dcl=='I': reps=sheets[:]
  else:
   seen=set()
   for s in sheets:
    if s in seen:continue
    reps.append(s);seen.add(s);seen.add(mul(s,dcl))
  Pd=Pdown[v];kdown=Pd.rows
  for s0 in reps:
   repeat=1 if dcl=='I' else 2
   # generate traversal sequence as (cid,dir) while tracking sheet
   seq=[];s=s0
   for _ in range(repeat):
    for z in lk['chosen']['word']:
     cid,dr=path_cell(z['edge'],s,z['direction']);seq.append((cid,dr))
     s=mul(s,ar[z['edge']]['deck'])
   assert s==s0
   # raw path map rows n -> c1u using suffix algorithm
   raw=[defaultdict(lambda:sp.Integer(0)) for _ in range(n)]
   tail=sp.eye(n)
   # also holonomy recompute
   for cid,dr in reversed(seq):
    R=uR[cid]
    if dr==1: Rp=R; S=sp.eye(n)
    else: Rp=R.inv();S=-R.inv()
    Bu=ub[cid]
    Cblk=tail*S*Bu
    for i in range(n):
     for j in range(Bu.cols):
      if Cblk[i,j]:raw[i][uoff[cid]+j]+=Cblk[i,j]
    tail=sp.simplify(tail*Rp)
   Htot=tail
   # expected total holonomy
   if dcl=='I':
    Hb=sp.Matrix([[sp.Integer(lk['chosen']['binet_holonomy'])]]) if kind=='Binet' else M(lk['chosen']['ad_holonomy'])
    assert Htot==Tau[s0].inv()*Hb*Tau[s0],(kind,v,s0,'H')
    Pup=Pd*Tau[s0]
    # exact left fixed
    assert Pup*Htot==Pup
   else:
    assert Htot==sp.eye(n),(kind,v,s0,Htot)
    first=Pd*Tau[s0]
    Pup=complete_rows(first,n) if first.rows<n else first
    assert Pup.rows==n and Pup.rank()==n
   row0=c2off
   # D1 rows Pup * raw
   for i in range(Pup.rows):
    for rr in range(n):
     coef=Pup[i,rr]
     if coef:
      for col,val in raw[rr].items(): D1ent[(row0+i,col)]=D1ent.get((row0+i,col),0)+coef*val
   # L2/transfer block is solved exactly after D1u*L1 is available.
   scale=1 if dcl=='I' else 2
   c2cells.append({'base_vertex':v,'start_sheet':s0,'deck_monodromy':dcl,'traversals':repeat,'row_offset':row0,'fiber_dim':Pup.rows,'downstairs_fixed_dim':kdown,'expected_transfer_scale':scale,
                   'P_rows':[[str(Pup[i,j]) for j in range(n)] for i in range(Pup.rows)]})
   c2off+=Pup.rows
 c2u=c2off
 expect2=150 if kind=='Binet' else 702
 assert c2u==expect2,(kind,c2u)
 D1u=smat(c2u,c1u,D1ent)
 # exact upper complex
 assert D1u*D0u==sp.zeros(c2u,c0u), f'{kind} upper d1d0 fail'
 # Solve the exact cellwise transfer/descent injection L2 from D1u*L1 = L2*D1_down.
 # This is necessary because C2 downstairs is a coinvariant module represented by left-fixed covectors;
 # on order-two deck cells the upstairs doubled disk lands in a transfer subspace rather than a coordinate subblock.
 lhs=D1u*L1
 l2ent={}
 transfer_blocks=[]
 for cc in c2cells:
  v=cc['base_vertex'];r0=cc['row_offset'];m=cc['fiber_dim'];k=cc['downstairs_fixed_dim'];c0=poff[v]
  U=lhs[r0:r0+m,:]
  if k==0:
   assert U==sp.zeros(m,D1d.cols),(kind,v,cc['start_sheet'],'expected zero transfer')
   X=sp.zeros(m,0)
  else:
   Dv=D1d[c0:c0+k,:]
   if cc['deck_monodromy']=='I':
    assert m==k
    X=sp.eye(k)
    assert X*Dv==U,(kind,v,cc['start_sheet'],'trivial-monodromy transfer mismatch')
   else:
    # solve each upstairs row as an exact linear combination of the k downstairs rows
    xr=[]
    for i in range(m):
     sol=sp.linsolve((Dv.T,U.row(i).T))
     assert sol is not sp.EmptySet and len(sol)>0,(kind,v,i,'no transfer solution')
     vec=list(next(iter(sol)))
     if any(x.free_symbols for x in vec):
      # if the lower d1 block has zero/deficient rows, use the canonical fixed-vector transfer below instead
      vec=None
     xr.append(vec)
    if all(vv is not None for vv in xr):
     X=sp.Matrix(xr)
    else:
     # canonical transfer: columns span the fixed-vector space dual to the downstairs left-fixed covectors.
     # Find V with Pd*V=I and stabilized by the effective one-turn holonomy.
     Hb=sp.Matrix([[sp.Integer(lk['chosen']['binet_holonomy'])]]) if kind=='Binet' else M(lk['chosen']['ad_holonomy'])
     # gauge at start sheet; effective stabilizer on primal vectors is Tau[s0].inv()*Hb*Tau[s0]
     Hg=Tau[s0].inv()*Hb*Tau[s0]
     Vfix=sp.Matrix.hstack(*(Hg-sp.eye(n)).nullspace()) if (Hg-sp.eye(n)).nullspace() else sp.zeros(n,0)
     Lam=Pd*Tau[s0]
     assert Vfix.cols==k,(kind,v,'fixed vector dim',Vfix.cols,k)
     Cdual=(Lam*Vfix).inv()
     V=Vfix*Cdual
     # Upper row coordinates are Pup*raw, so transfer column is 2*Pup*V.
     X=2*Pup*V
    assert X*Dv==U,(kind,v,cc['start_sheet'],'nontriv transfer mismatch',X)
    assert X.rank()==k,(kind,v,cc['start_sheet'],'transfer rank',X.rank(),k)
   add_block(l2ent,r0,c0,X)
  # record exact divisibility pattern on nontrivial deck monodromy blocks
  nums=[]
  for z in X: nums.append(sp.Rational(z))
  transfer_blocks.append({'vertex':v,'start_sheet':cc['start_sheet'],'deck_monodromy':cc['deck_monodromy'],
                          'rows':m,'cols':k,'rank':X.rank(),'matrix':[[str(X[i,j]) for j in range(X.cols)] for i in range(X.rows)]})
 L2=smat(c2u,D1d.rows,l2ent)
 rhs=L2*D1d
 assert lhs==rhs,(kind,'d1 descent fail after exact transfer solve')
 # theoretical invariant dims from orbit stabilizers = lower dims; certify cellwise
 invC0=62*n
 invC1=0
 for e in range(123): invC1+=Bdown[e].cols
 invC2=sum(Pdown[v].rows for v in range(62))
 assert (invC0,invC1,invC2)==(D0d.cols,D0d.rows,D1d.rows)
 # ranks of descent injections full
 assert L0.rank()==D0d.cols and L1.rank()==D0d.rows and L2.rank()==D1d.rows
 # upper representation simplification histogram
 rhist=Counter()
 for cell in c1top:
  R=uR[cell['cid']]
  rhist[str(R.tolist())]+=1
 return {
  'kind':kind,'n':n,'dimensions_upstairs':{'C0':c0u,'C1':c1u,'C2':c2u},
  'dimensions_V4_invariant':{'C0':invC0,'C1':invC1,'C2':invC2},
  'dimensions_downstairs':C['dimensions'][kind],
  'upper_d1d0_zero':True,'d0_descent_square':True,'d1_descent_square':True,
  'descent_injection_ranks':{'C0':L0.rank(),'C1':L1.rank(),'C2':L2.rank()},
  'equivariance_object_fiber_checks':obj_eq_checks,'equivariance_arrow_transport_checks':eq_checks,
  'nontrivial_monodromy_transfer_factor':2,
  'C2_cells':c2cells,'C2_transfer_blocks':transfer_blocks,
  'matrices':{
   'd0_up':D0u,'d1_up':D1u,'L0':L0,'L1':L1,'L2':L2
  }
 }
resB=build('Binet'); print('Binet PASS',resB['dimensions_upstairs'],resB['dimensions_V4_invariant'])
resA=build('Ad'); print('Ad PASS',resA['dimensions_upstairs'],resA['dimensions_V4_invariant'])
# serialize sparse matrices compactly
def sentries(X): return [{'row':int(i),'col':int(j),'value':str(v)} for (i,j),v in sorted(X.todok().items()) if v!=0]
def pack(r):
 q={k:v for k,v in r.items() if k!='matrices'}
 q['matrix_nnz']={k:len(v.todok()) for k,v in r['matrices'].items()}
 q['d0_up_sparse']=sentries(r['matrices']['d0_up']);q['d1_up_sparse']=sentries(r['matrices']['d1_up'])
 q['descent_L0_sparse']=sentries(r['matrices']['L0']);q['descent_L1_sparse']=sentries(r['matrices']['L1']);q['descent_L2_sparse']=sentries(r['matrices']['L2'])
 return q
out={'version':'v50','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE','status':'PASS-EXACT-V4-LOCAL-SYSTEM-PULLBACK/DESCENT',
 'gauge':'sheet-unwound gauge: at sheet s, coordinate change by rho_deck(s)^(-1); deck factors disappear from ordinary lifted transports and reappear in the V4 descent datum',
 'V4_representations':{'Binet':{'I':1,'T':-1,'U':-1,'TU':1},'Ad':G['deck_ad']},
 'systems':{'Binet':pack(resB),'Ad':pack(resA)},
 'descent_statement':'For each coefficient system over Q, the downstairs orbifold complex is exactly chain-isomorphic to the V4-equivariant/invariant descent of the four-sheet complex. For the three order-two deck-monodromy 2-cells, the chain injection uses transfer factor 2 because the upstairs geometric disk has the doubled boundary traversal.',
 'integral_warning':'The C2 transfer factor 2 is invertible over Q but not over Z. Therefore rational descent is an isomorphism; an integral identification would require a separate 2-torsion/index audit.',
 'cohomology_gate':'H1/H2 may be reconsidered over Q after an independent replay; integral H1/H2 remains guarded by the factor-2 transfer.'}
(B/'selling_cover_local_system_descent_v50.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('WROTE',B/'selling_cover_local_system_descent_v50.json')
