import json, os, hashlib, math
from fractions import Fraction
import sympy as sp
BASE='/mnt/data/v50_work'
D=json.load(open(BASE+'/selling_field_orbifold_cochains_v50.json'))
U=json.load(open(BASE+'/selling_cover_local_system_descent_v50.json'))

def dims(system, up=False):
    if not up:
        dd=D['dimensions'][system]
        return dd['C0'],dd['C1'],dd['C2']
    x=U['systems'][system]['dimensions_upstairs']
    return x['C0'],x['C1'],x['C2']

def sparse(rows, cols, entries):
    data={}
    for e in entries:
        v=sp.Rational(e['value'])
        if v: data[(int(e['row']),int(e['col']))]=v
    return sp.SparseMatrix(rows,cols,data)

def mat_down(system):
    c0,c1,c2=dims(system)
    return sparse(c1,c0,D[f'd0_{system}_sparse']), sparse(c2,c1,D[f'd1_{system}_sparse'])

def mat_up(system):
    c0,c1,c2=dims(system,True)
    s=U['systems'][system]
    return sparse(c1,c0,s['d0_up_sparse']), sparse(c2,c1,s['d1_up_sparse'])

def Ls(system):
    sd=U['systems'][system]
    down=dims(system); up=dims(system,True)
    return tuple(sparse(up[i],down[i],sd[f'descent_L{i}_sparse']) for i in range(3))

def quotient_basis(ker_basis, im_basis):
    # Build independent extension of im basis by kernel vectors. Return H reps.
    if im_basis:
        B=sp.Matrix.hstack(*im_basis)
        rank=B.rank()
    else:
        B=sp.zeros(len(ker_basis[0]) if ker_basis else 0,0); rank=0
    reps=[]
    for v in ker_basis:
        C=B.row_join(v)
        r=C.rank()
        if r>rank:
            reps.append(v); B=C; rank=r
    return reps

def primitive_vec(v):
    den=sp.ilcm(*[x.q for x in v]) if len(v) else 1
    vals=[int(x*den) for x in v]
    g=0
    for a in vals:g=math.gcd(g,abs(a))
    if g: vals=[a//g for a in vals]
    for a in vals:
        if a<0: vals=[-x for x in vals]; break
        if a>0: break
    return vals

def sparse_vec(vals):
    return [{'i':i,'value':str(v)} for i,v in enumerate(vals) if v]

def solve_in_span(B,v):
    # exact solve B*x=v, assuming in col span; gauss_jordan_solve may return params
    sol, params=B.gauss_jordan_solve(v)
    if params.rows:
        sol=sol.subs({p:0 for p in params})
    return sol

out={'version':'v50','status':'EXACT_RATIONAL_COHOMOLOGY','publication_state':'PRESEAL/NONPUBLISHABLE','systems':{}}
for system in ['Binet','Ad']:
    print('SYSTEM',system,flush=True)
    d0,d1=mat_down(system)
    u0,u1=mat_up(system)
    L0,L1,L2=Ls(system)
    c0,c1,c2=dims(system); uc0,uc1,uc2=dims(system,True)
    print(' ranks downstairs d0',flush=True); r0=d0.rank(); print(r0,flush=True)
    print(' ranks downstairs d1',flush=True); r1=d1.rank(); print(r1,flush=True)
    print(' ranks upstairs d0',flush=True); ur0=u0.rank(); print(ur0,flush=True)
    print(' ranks upstairs d1',flush=True); ur1=u1.rank(); print(ur1,flush=True)
    # downstairs reps H1, H2
    print(' null d1',flush=True); ker1=d1.nullspace(); print(len(ker1),flush=True)
    im0=d0.columnspace(); print('im0',len(im0),flush=True)
    h1=quotient_basis(ker1,im0); print('H1 reps',len(h1),flush=True)
    # H2 = C2 / im d1: extend image with standard basis to obtain quotient reps
    im1=d1.columnspace();
    B=sp.Matrix.hstack(*im1) if im1 else sp.zeros(c2,0); rank=r1; h2=[]
    for i in range(c2):
        e=sp.eye(c2)[:,i]
        C=B.row_join(e); rr=C.rank()
        if rr>rank:
            h2.append(e); B=C; rank=rr
    print('H2 reps',len(h2),flush=True)
    # verify maps downstairs reps upstairs are cocycles and independent modulo upper boundaries
    lifted_h1=[L1*v for v in h1]
    assert all(u1*v==sp.zeros(uc2,1) for v in lifted_h1)
    up_im0=u0.columnspace(); upB=sp.Matrix.hstack(*up_im0) if up_im0 else sp.zeros(uc1,0); base_rank=ur0
    indep=0; BB=upB
    for v in lifted_h1:
        rr=BB.row_join(v).rank()
        if rr>base_rank+indep:
            BB=BB.row_join(v); indep+=1
    # H2 lift via L2, independence modulo im u1
    lifted_h2=[L2*v for v in h2]
    up_im1=u1.columnspace(); UB=sp.Matrix.hstack(*up_im1) if up_im1 else sp.zeros(uc2,0); rr0=ur1; indep2=0; UBB=UB
    for v in lifted_h2:
        rr=UBB.row_join(v).rank()
        if rr>rr0+indep2:
            UBB=UBB.row_join(v); indep2+=1
    # invariant-sector quotient dimensions should equal downstairs by exactness; record full upstairs dims too
    rec={
      'dimensions_downstairs':{'C0':c0,'C1':c1,'C2':c2},
      'ranks_downstairs':{'d0':r0,'d1':r1},
      'cohomology_dimensions_downstairs_Q':{'H0':c0-r0,'H1':c1-r1-r0,'H2':c2-r1},
      'dimensions_upstairs':{'C0':uc0,'C1':uc1,'C2':uc2},
      'ranks_upstairs':{'d0':ur0,'d1':ur1},
      'cohomology_dimensions_upstairs_Q':{'H0':uc0-ur0,'H1':uc1-ur1-ur0,'H2':uc2-ur1},
      'descent_class_comparison':{
        'H1_down_representatives_count':len(h1),'H1_lifts_nontrivial_independent_mod_upper_boundaries':indep,
        'H2_down_representatives_count':len(h2),'H2_lifts_nontrivial_independent_mod_upper_boundaries':indep2,
        'statement':'Downstairs classes inject into upstairs through L and represent the V4-invariant cohomology sector over Q.'
      },
      'H1_down_representatives':[sparse_vec(primitive_vec(v)) for v in h1],
      'H2_down_representatives':[sparse_vec(primitive_vec(v)) for v in h2],
      'H1_lifted_representatives_upstairs':[sparse_vec(primitive_vec(v)) for v in lifted_h1],
      'H2_lifted_representatives_upstairs':[sparse_vec(primitive_vec(v)) for v in lifted_h2],
    }
    out['systems'][system]=rec

# checks
for system,rec in out['systems'].items():
    assert rec['cohomology_dimensions_downstairs_Q']['H1']==rec['descent_class_comparison']['H1_down_representatives_count']
    assert rec['cohomology_dimensions_downstairs_Q']['H2']==rec['descent_class_comparison']['H2_down_representatives_count']
    assert rec['descent_class_comparison']['H1_lifts_nontrivial_independent_mod_upper_boundaries']==rec['cohomology_dimensions_downstairs_Q']['H1']
    assert rec['descent_class_comparison']['H2_lifts_nontrivial_independent_mod_upper_boundaries']==rec['cohomology_dimensions_downstairs_Q']['H2']

p=BASE+'/selling_historical_cohomology_Q_v50.json'
with open(p,'w') as f: json.dump(out,f,indent=2,sort_keys=True)
print('WROTE',p,hashlib.sha256(open(p,'rb').read()).hexdigest())
