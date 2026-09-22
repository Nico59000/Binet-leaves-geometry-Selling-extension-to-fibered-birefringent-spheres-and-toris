#!/usr/bin/env python3
from __future__ import annotations
import json, csv, hashlib, sys
from pathlib import Path
import sympy as sp
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from selling_reconstruction_engine_v49 import gram_from_six, six_from_gram, S, RELABEL, numerical_form_key, homogeneous

# ------------------------------------------------------------------
# Base cell and source transitions
# ------------------------------------------------------------------
Q=gram_from_six(7,-1,-8,0,2,0)
trans={'h':(2,8),'m':(5,12),'g':(9,11)}
wall_hom={'h':'h','m':'m','g':'g'}

# target representatives: first element of each source pair
H=S[1].T*Q*S[1]
M=S[4].T*Q*S[4]
G=S[8].T*Q*S[8]
REP={'Q':Q,'H':H,'M':M,'G':G}

# exact source-algebraic constants
qA=-4+2*sp.sqrt(210)/7
rB=1-sp.sqrt(195)/15
q=sp.symbols('q')
pq=469*q**4+4564*q**3+7224*q**2+1984*q-416
qD_interval=(sp.Rational(44135388035037,324004210501480),sp.Rational(506570030731057,3718803213884307))
rD_expr=-sp.Rational(469,1620)*q**3-sp.Rational(91,36)*q**2-sp.Rational(553,270)*q+sp.Rational(158,405)

u=sp.symbols('u')
pG=567*u**4+162*u**3-114*u**2-36*u-2
uG_interval=(sp.Rational(-35865665499,481325503693),sp.Rational(-4671184207244,62688369506407))
vG_expr=(16-756*u**3-90*u**2+162*u)/5

# ------------------------------------------------------------------
# Projective marking torsor: A1 and A2 differ only by target relabel + -I.
# ------------------------------------------------------------------
def unique_relabel(N1,N2):
    sol=[B for B in RELABEL if B.T*N1*B==N2]
    assert len(sol)==1
    return sol[0]

def wall_perm(B):
    aa,bb,cc,gg,hh,kk=sp.symbols('a b c g h k')
    X=gram_from_six(aa,bb,cc,gg,hh,kk)
    old=homogeneous(X); new=homogeneous(B.T*X*B)
    out={}
    for nl,e in new.items():
        matches=[ol for ol,oe in old.items() if sp.expand(e-oe)==0]
        assert len(matches)==1
        out[nl]=matches[0]
    return out

torsor={}
for w,(i,j) in trans.items():
    A1,A2=S[i-1],S[j-1]
    N1,N2=A1.T*Q*A1,A2.T*Q*A2
    B=unique_relabel(N1,N2)
    rel=sp.simplify(A2.inv()*A1*B)
    perm=wall_perm(B)
    torsor[w]={
      'candidates':[f'S{i}',f'S{j}'],
      'candidate1_target_six':[str(x) for x in six_from_gram(N1)],
      'candidate2_target_six':[str(x) for x in six_from_gram(N2)],
      'unique_target_relabelling':B.tolist(),
      'relabel_det':int(B.det()),
      'relabel_involution':bool(B*B==sp.eye(3)),
      'wall_permutation_target2_labels_from_target1':perm,
      'crossed_wall_fixed':perm[w]==w,
      'projective_identity_matrix':rel.tolist(),
      'projective_identity_is_minus_I':bool(rel==-sp.eye(3)),
      'decision':'UNIQUE-PROJECTIVE-TRANSPORT-CLASS'
    }

# ------------------------------------------------------------------
# Cell boundaries for the three neighbours, in exact charts.
# Only minpoly+isolating interval is used for quartic vertices.
# ------------------------------------------------------------------
cell_Q={
 'id':'F_Q2','representative_six':[str(x) for x in six_from_gram(Q)],
 'vertices':['E','A','D','B'],
 'boundary':[('EA','E','A','mirror(g+k)'),('AD','A','D','h'),('DB','D','B','m'),('BE','B','E','g')]
}

# H chart: x=1, coordinates (u=y/x,v=z/x)
tH=sp.Rational(3,4)-3*sp.sqrt(5)/4
sH=-sp.Rational(7,2)+3*sp.sqrt(5)/2
cell_H={
 'id':'F_H','representative_six':[str(x) for x in six_from_gram(H)],'chart':'x=1;(u=y/x,v=z/x)',
 'vertices':{
   'A':{'uv':['-1',str(-qA)],'vanishing':['g','h','l'],'kind':'AXIS-TRIPLE'},
   'C':{'uv':[str(tH),str(sH)],'vanishing':['g','m'],'kind':'CROSSING'},
   'D':{'uv':['r_D-1','-q_D'],'q_D_minpoly':str(pq),'q_D_interval':[str(x) for x in qD_interval],'r_D_formula':str(rD_expr),'vanishing':['h','m'],'kind':'FISSURE'}},
 'boundary_order':['A','C','D','A'],
 'edges':[('AC','A','C','g'),('CD','C','D','m'),('DA','D','A','h')]
}

# M chart: z=1, coordinates (u=x/z,v=y/z)
tM=-sp.Rational(2,3)+sp.sqrt(5)/3
cell_M={
 'id':'F_M','representative_six':[str(x) for x in six_from_gram(M)],'chart':'z=1;(u=x/z,v=y/z)',
 'vertices':{
   'D':{'uv':['(r_D-q_D)/(r_D-1)','r_D/(1-r_D)'],'q_D_minpoly':str(pq),'q_D_interval':[str(x) for x in qD_interval],'r_D_formula':str(rD_expr),'vanishing':['h','m'],'kind':'FISSURE'},
   'C':{'uv':[str(tM),str(tM)],'vanishing':['g','h'],'kind':'CROSSING'},
   'B':{'uv':['-r_B/(1-r_B)','r_B/(1-r_B)'],'r_B':str(rB),'vanishing':['g','m'],'kind':'CROSSING'}},
 'boundary_order':['D','C','B','D'],
 'edges':[('DC','D','C','h'),('CB','C','B','g'),('BD','B','D','m')]
}

# G chart: x=1, coordinates (u=y/x,v=z/x)
vA_G=-sp.Rational(4,3)+2*sp.sqrt(10)/3
cell_G={
 'id':'F_G','representative_six':[str(x) for x in six_from_gram(G)],'chart':'x=1;(u=y/x,v=z/x)',
 'vertices':{
   'B':{'uv':[str(-rB),'0'],'vanishing':['g','m'],'kind':'CROSSING'},
   'E':{'uv':['0','0'],'vanishing':['g','k'],'kind':'AXIS-SEAM-JUNCTION'},
   'GA':{'uv':['0',str(vA_G)],'vanishing':['g','h','k'],'kind':'AXIS-TRIPLE'},
   'GD':{'uv':['u_G','v_G(u_G)'],'u_G_minpoly':str(pG),'u_G_interval':[str(x) for x in uG_interval],'v_G_formula':str(vG_expr),'vanishing':['h','m'],'kind':'FISSURE'}},
 'boundary_order':['B','E','GA','GD','B'],
 'edges':[('BE','B','E','g'),('EGA','E','GA','mirror(g+k)'),('GAGD','GA','GD','h'),('GDB','GD','B','m')]
}

# ------------------------------------------------------------------
# Exact H<->M gluing. First candidate plus unique relabelling to the M/H reps.
# ------------------------------------------------------------------
A_HM=S[4]  # cross m from H
N=A_HM.T*H*A_HM
B_HM=unique_relabel(N,M)
C_HM=sp.simplify(A_HM*B_HM)
A_MH=S[1]
N2=A_MH.T*M*A_MH
B_MH=unique_relabel(N2,H)
C_MH=sp.simplify(A_MH*B_MH)
assert C_HM.T*H*C_HM==M
assert C_MH.T*M*C_MH==H
assert C_MH==C_HM.inv()

# exact endpoint C maps H chart -> M chart
vHC=sp.Matrix([1,tH,sH]); yC=C_HM.T*vHC
assert sp.simplify(yC[0]/yC[2]-tM)==0 and sp.simplify(yC[1]/yC[2]-tM)==0
# D map checked modulo p(q): H direction [1,rD-1,-qD] -> M expected chart
r=rD_expr
vHD=sp.Matrix([1,r-1,-q]); yD=sp.simplify(C_HM.T*vHD)
expr1=sp.together(yD[0]/yD[2]-(r-q)/(r-1));expr2=sp.together(yD[1]/yD[2]-r/(1-r))
assert sp.rem(sp.Poly(sp.together(expr1).as_numer_denom()[0],q),sp.Poly(pq,q))==0
assert sp.rem(sp.Poly(sp.together(expr2).as_numer_denom()[0],q),sp.Poly(pq,q))==0

# ------------------------------------------------------------------
# Local 4-face patch incidence.
# Global edge orientation is chosen so all shared edges cancel face-to-face.
# ------------------------------------------------------------------
V=['E','A','D','B','C','GA','GD']
E=['EA','AD','DB','BE','AC','CD','CB','BGD','GDGA','GAE']
edge_st={
 'EA':('E','A'),'AD':('A','D'),'DB':('D','B'),'BE':('B','E'),
 'AC':('A','C'),'CD':('C','D'),'CB':('C','B'),
 'BGD':('B','GD'),'GDGA':('GD','GA'),'GAE':('GA','E')}
faces={
 'Q':{'EA':1,'AD':1,'DB':1,'BE':1},
 'H':{'AC':1,'CD':1,'AD':-1},
 'M':{'CD':-1,'CB':1,'DB':-1},
 'G':{'BE':-1,'BGD':1,'GDGA':1,'GAE':1},
}
D1=sp.zeros(len(V),len(E))
for j,e in enumerate(E):
 s,t=edge_st[e];D1[V.index(s),j]=-1;D1[V.index(t),j]=1
D2=sp.zeros(len(E),len(faces))
for j,(f,bd) in enumerate(faces.items()):
 for e,c in bd.items():D2[E.index(e),j]=c
assert D1*D2==sp.zeros(len(V),len(faces))
outer=D2*sp.ones(len(faces),1)
outer_expected={'EA':1,'AC':1,'CB':1,'BGD':1,'GDGA':1,'GAE':1}
assert all(int(outer[E.index(e),0])==outer_expected.get(e,0) for e in E)

# outgoing ordinary arcs of the patch and numerical neighbours
outgoing={}
for fid,Form,wall in [('H',H,'g'),('M',M,'g'),('G_h',G,'h'),('G_m',G,'m')]:
 con={'g':'p23','h':'p13','m':'p02'}[wall]
 # first source representative
 cand={'p23':9,'p13':2,'p02':5}[con]
 A=S[cand-1]; N=A.T*Form*A
 outgoing[fid]={'wall':wall,'move_rep':f'S{cand}','target_six':[str(x) for x in six_from_gram(N)],'numerical_key':[str(x) for x in numerical_form_key(N)]}
assert len({tuple(v['numerical_key']) for v in outgoing.values()})==4

checks={
 'three_torsors_projectively_resolved':all(x['projective_identity_is_minus_I'] and x['crossed_wall_fixed'] for x in torsor.values()),
 'H_M_transport_inverse':C_MH==C_HM.inv(),
 'H_M_form_gluing_exact':C_HM.T*H*C_HM==M,
 'H_M_C_vertex_exact':True,
 'H_M_D_vertex_exact_mod_minpoly':True,
 'local_D1D2_zero':D1*D2==sp.zeros(len(V),len(faces)),
 'local_patch_euler_char_one':len(V)-len(E)+len(faces)==1,
 'local_patch_outer_boundary_six_edges':sum(1 for e in E if outer[E.index(e),0]!=0)==6,
 'four_distinct_outgoing_numerical_neighbours':len({tuple(v['numerical_key']) for v in outgoing.values()})==4,
}
assert all(checks.values())

out={
 'version':'v50','status':'PASS-LOCAL-4-FACE-PATCH','publication_state':'RESEARCH-CHECKPOINT/PRESEAL/NONPUBLISHABLE',
 'projective_marking_torsor':torsor,
 'cells':{'Q':cell_Q,'H':cell_H,'M':cell_M,'G':cell_G},
 'H_M_exact_glue':{'C_HM':C_HM.tolist(),'det':int(C_HM.det()),'C_MH':C_MH.tolist(),'inverse_exact':True},
 'local_complex':{'C0':V,'C1':E,'C2':list(faces),'edge_source_target':edge_st,'face_boundaries':faces,
                  'D1_vertices_by_edges':D1.tolist(),'D2_edges_by_faces':D2.tolist(),
                  'euler_characteristic':len(V)-len(E)+len(faces),
                  'outer_boundary':[e for e in E if outer[E.index(e),0]!=0],
                  'outer_boundary_cycle':['E','A','C','B','GD','GA','E']},
 'outgoing_gate':outgoing,
 'global_status':{
   'local_patch':'PROVEN-LOCAL-4-FACE-SUBCOMPLEX',
   'projective_transport_classes':'PROVEN-UNIQUE',
   'full_TU_quotient_closure':'NT-NOT-YET-REACHED',
   'global_C0_C1_C2':'NT-GUARDED',
   'global_d1d2':'NOT-RUN/GUARDED','twisted_d1d0':'NOT-RUN/GUARDED','historical_H1_H2':'NT/GUARDED'},
 'checks':checks
}

p=HERE/'selling_Q2_local_patch_v50.json';p.write_text(json.dumps(out,indent=2,ensure_ascii=False,default=str),encoding='utf-8')
c=HERE/'selling_Q2_local_patch_incidence_v50.csv'
with c.open('w',newline='',encoding='utf-8') as f:
 wr=csv.writer(f);wr.writerow(['edge','source','target','Q','H','M','G'])
 for e in E:wr.writerow([e,*edge_st[e],*[faces[F].get(e,0) for F in ['Q','H','M','G']]])
print(json.dumps({'status':'PASS','checks':checks,'json_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'csv_sha256':hashlib.sha256(c.read_bytes()).hexdigest()},indent=2))
