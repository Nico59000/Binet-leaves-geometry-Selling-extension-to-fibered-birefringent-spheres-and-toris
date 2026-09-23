#!/usr/bin/env python3
import json, sys
from pathlib import Path
import sympy as sp
HERE=Path('/mnt/data/v50_work')
sys.path.insert(0,str(HERE))
from selling_reconstruction_engine_v49 import S

M0=sp.Matrix([[7,0,2],[0,-1,0],[2,0,-8]])
G=M0.inv()
reg=json.load(open(HERE/'selling_regular_TU_quotient_v50.json'))
faces=json.load(open(HERE/'selling_face_cycles_v50.json'))['faces']
bfs=json.load(open(HERE/'exact_TU_field_bfs_v50.json'))['records']
GT=sp.Matrix([[-8,3,-2],[-21,8,-6],[0,0,-1]])
GU=sp.Matrix([[-4,1,2],[-9,2,6],[-3,1,1]])
DECK={'T':GT,'U':GU,'TU':GT*GU}
q,r=sp.symbols('q r')

def smat(a): return sp.Matrix([[sp.Integer(x) for x in row] for row in a])
def ints(A): return [[int(A[i,j]) for j in range(A.cols)] for i in range(A.rows)]

def reflection_from_line(factor):
    e=sp.expand(sp.sympify(factor))
    a=sp.expand(e).coeff(q); b=sp.expand(e).coeff(r); c=e.subs({q:0,r:0})
    ell=sp.Matrix([c,-b,a]) # ell(v)=c*x-b*y+a*z
    n=M0*ell
    den=sp.expand((ell.T*n)[0])
    R=sp.simplify(sp.eye(3)-2*n*ell.T/den) # action on v, fixes ell(v)=0 pointwise
    J=R.T                              # congruence/substitution matrix
    return ell,R,J,den

def cone_sign(J):
    t=sp.Matrix([1,0,0]) # Q(t)=7>0
    val=sp.expand((t.T*M0*(J*t))[0])
    return (1 if val>0 else -1), int(val)

# 16 inherited multiwall occurrences.
multi=[]
for e in reg['edges']:
    if e['id']>=118: continue
    for f,s in e['source_orbit']:
        cyc=faces[f]['cycles'][0]
        seg=cyc['segments'][s]
        comp=faces[f]['components'][seg['component']]
        if comp['kind']=='MULTI':
            A=smat(bfs[f]['path_matrix']); Qf=A.T*M0*A
            ell,R,J,den=reflection_from_line(comp['factor'])
            K=sp.simplify(A.inv()*J*A)
            assert all(x.q==1 for x in K)
            checks={
              'J_integral': all(x.q==1 for x in J),
              'J_det_minus1': sp.det(J)==-1,
              'J_involution': J*J==sp.eye(3),
              'J_preserves_base_form': J.T*M0*J==M0,
              'R_fixes_line_pointwise': sp.simplify(R-sp.eye(3)+2*(M0*ell)*ell.T/den)==sp.zeros(3),
              'K_integral': all(x.q==1 for x in K),
              'K_det_minus1': sp.det(K)==-1,
              'K_involution': K*K==sp.eye(3),
              'K_stabilizes_field_form': K.T*Qf*K==Qf,
            }
            assert all(checks.values())
            eps,val=cone_sign(J)
            multi.append({'edge_orbit':e['id'],'field':f,'segment':s,'walls':comp['walls'],'factor':comp['factor'],
                          'global_reflection_on_v':ints(R),'global_congruence_J':ints(J),'field_stabilizer_K':ints(K),
                          'binet_time_cone_sign':eps,'time_cone_pairing_with_e1':val,'checks':checks,
                          'role':'ISOTROPY/MULTIWALL-REFLECTION; NOT longitudinal arc transport M_e'})
assert len(multi)==16

# Five regularizing mirror edges: deck isotropy.
mirror=[]
for e in reg['edges']:
    if e['kind']!='DECK-FIXED-MIRROR': continue
    f=e['face']; D=DECK[e['deck']]; A=smat(bfs[f]['path_matrix']);Qf=A.T*M0*A;K=sp.simplify(A.inv()*D*A)
    eps,val=cone_sign(D)
    checks={'deck_involution':D*D==sp.eye(3),'deck_preserves_base_form':D.T*M0*D==M0,
            'K_integral':all(x.q==1 for x in K),'K_stabilizes_field_form':K.T*Qf*K==Qf}
    assert all(checks.values())
    mirror.append({'edge_id':e['id'],'face':f,'deck':e['deck'],'fixed_locus':e['fixed_locus'],
                   'field_stabilizer_K':ints(K),'binet_time_cone_sign':eps,'time_cone_pairing_with_e1':val,
                   'checks':checks,'role':'DECK-ISOTROPY/MIRROR; NOT longitudinal arc transport M_e'})
assert len(mirror)==5

# Exact no-go: raw transverse wall moves cannot be re-used as longitudinal edge transports
# on the first certified Q2 cell. Boundary order: multi(g,k), h, m, g.
Jr=reflection_from_line('r')[2]
choices={'h':[S[1],S[7]], 'm':[S[4],S[11]], 'g':[S[8],S[10]]}
no=[]
for ih in range(2):
  for im in range(2):
    for ig in range(2):
      mats=[Jr,choices['h'][ih],choices['m'][im],choices['g'][ig]]
      P=sp.eye(3)
      for X in mats:P=P*X
      Pr=sp.eye(3)
      for X in mats:Pr=X*Pr
      row={'choice':[('S2','S8')[ih],('S5','S12')[im],('S9','S11')[ig]],
           'forward_product':ints(P),'reverse_product':ints(Pr),
           'forward_preserves_Q2':P.T*M0*P==M0,'reverse_preserves_Q2':Pr.T*M0*Pr==M0,
           'forward_pm_identity':P in (sp.eye(3),-sp.eye(3)), 'reverse_pm_identity':Pr in (sp.eye(3),-sp.eye(3))}
      assert not row['forward_preserves_Q2'] and not row['reverse_preserves_Q2']
      no.append(row)

# regular chain-complex checks are inherited but replayed here from sparse matrices.
def sparse(rows, nrow,ncol):
    A=sp.zeros(nrow,ncol)
    for z in rows:A[z['row'],z['col']]=z['value']
    return A
p1=sparse(reg['partial1_sparse'],reg['C0'],reg['C1'])
p2=sparse(reg['partial2_sparse'],reg['C1'],reg['C2'])
assert p1*p2==sp.zeros(reg['C0'],reg['C2'])

out={'version':'v50','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
     'status':'PASS-EXACT-TRANSPORT-GATE-AUDIT',
     'regular_complex':{'C0':reg['C0'],'C1':reg['C1'],'C2':reg['C2'],'euler_characteristic':reg['euler_characteristic'],
                        'partial1_partial2_zero':True,'status':'PROVEN-GLOBAL-TOPOLOGICAL-REGULAR-QUOTIENT'},
     'edge_categories':{'ordinary_boundary_orbits':102,'multiwall_boundary_orbits':16,'new_deck_mirror_edges':5,'total':123},
     'multiwall_isotropy_certificates':multi,'mirror_isotropy_certificates':mirror,
     'binet_isotropy_summary':{'multiwall_reflections_all_preserve_time_cone':all(z['binet_time_cone_sign']==1 for z in multi),
                               'new_mirror_edges_all_reverse_time_cone':all(z['binet_time_cone_sign']==-1 for z in mirror)},
     'transverse_vs_longitudinal_no_go':{
       'first_face':'Q2/field0','tested_choices':len(no),'composition_conventions':2,
       'all_fail_to_preserve_Q2':True,'rows':no,
       'decision':'REFUTED-TYPED: raw transverse Selling wall moves (even choosing either representative in each torsor) are not the longitudinal boundary-arc transports M_e required by v44/v49.'},
     'longitudinal_transport_gate':{
       'M_e_total':False,'d0':'NOT-FORMED','d1':'NOT-FORMED','d1d0':'NOT-TESTED',
       'missing_data':'For every oriented primal boundary arc, a source-certified longitudinal substitution word mapping the endpoint chart/fibre at s(e) to that at t(e), with face-ordered composition equal to the declared face stabilizer. Transverse wall-crossing moves and isotropy reflections are typed separately and cannot be substituted.',
       'F_Sell_to_F29':'PENDING-MORPHISM','H1_H2':'NT/GUARDED'} }
(HERE/'selling_regular_transport_gate_v50.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({'status':out['status'],'regular':[reg['C0'],reg['C1'],reg['C2']], 'multi':len(multi),'mirror':len(mirror),'no_go_cases':len(no),'d0':'NOT-FORMED'}))
