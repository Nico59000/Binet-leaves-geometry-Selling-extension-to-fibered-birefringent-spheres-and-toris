import numpy as np, math
from scipy.optimize import differential_evolution
import sympy as sp
Q0=sp.diag(12,-1,-5)
P1=sp.Matrix([[1,0,0],[3,1,0],[0,0,1]]);P2=sp.Matrix([[1,1,0],[0,1,0],[0,0,1]])
Q1=sp.Matrix([[1,0,0],[0,1,0],[1,0,1]]);Q2=sp.Matrix([[1,0,1],[0,1,0],[0,0,1]]);Q3=Q1;Q4=sp.Matrix([[1,0,2],[0,1,0],[0,0,1]])
R1=Q2;R2=Q1;R3=sp.Matrix([[-1,0,0],[0,-1,0],[-1,0,1]]);R4=sp.Matrix([[1,0,0],[0,1,0],[-1,0,1]]);R5=sp.Matrix([[1,0,-1],[0,1,0],[0,0,1]])
S1=sp.Matrix([[1,0,0],[1,1,0],[0,0,1]]);S2=sp.Matrix([[-1,1,0],[0,1,0],[0,0,-1]]);S3=sp.Matrix([[1,0,0],[-1,1,0],[0,0,1]])

def six(M):return [float(M[0,0]),float(M[1,1]),float(M[2,2]),float(M[1,2]),float(M[0,2]),float(M[0,1])]
def margin(M,seed=1):
 Mn=np.array(M.tolist(),float); Mi=np.linalg.inv(Mn)
 vals,V=np.linalg.eigh(Mi); pos=np.where(vals>1e-10)[0]; neg=np.where(vals<-1e-10)[0]
 # equation v^T M^-1 v=1; one positive eigenvalue expected
 if len(pos)!=1:return (None,None,None)
 p=pos[0]; neg=list(neg)
 def calc(uv,sgn):
  u0,u1=uv; rhs=1-vals[neg[0]]*u0*u0-vals[neg[1]]*u1*u1
  if rhs<=0:return (1e6,None,None)
  t=sgn*math.sqrt(rhs/vals[p]); q=np.zeros(3);q[p]=t;q[neg[0]]=u0;q[neg[1]]=u1;v=V@q
  x,y,z=v;s=x+y+z;a,b,c,g,h,k=six(M)
  w=np.array([2*y*z-g,2*x*z-h,2*x*y-k,a+h+k-2*x*s,b+k+g-2*y*s,c+h+g-2*z*s])
  return w.max(),v,w
 best=(1e9,None,None);B=20
 for sg in [1,-1]:
  res=differential_evolution(lambda uv:calc(uv,sg)[0],[(-B,B),(-B,B)],seed=seed,maxiter=250,popsize=12,tol=1e-10,polish=True)
  r=calc(res.x,sg)
  if r[0]<best[0]:best=r
 return best

routes={}
for name,factors in {
 'P':[('P1',P1),('P2',P2)],
 'Q':[('Q1',Q1),('Q2',Q2),('Q3',Q3),('Q4',Q4)],
}.items():
 A=sp.eye(3)
 routes[name]=[('START',Q0)]
 for lab,F in factors:
  A=A*F;routes[name].append((lab,sp.simplify(A.T*Q0*A)))
# R starts after P, S starts after Q, because printed routes connect from targets
A=P1*P2
routes['PR']=[('Ptarget',A.T*Q0*A)]
for lab,F in [('R1',R1),('R2',R2),('R3',R3),('R4',R4),('R5',R5)]:
 A=A*F;routes['PR'].append((lab,sp.simplify(A.T*Q0*A)))
A=Q1*Q2*Q3*Q4
routes['QS']=[('Qtarget',A.T*Q0*A)]
for lab,F in [('S1',S1),('S2',S2),('S3',S3)]:
 A=A*F;routes['QS'].append((lab,sp.simplify(A.T*Q0*A)))
for rn,items in routes.items():
 print('\n',rn)
 for i,(lab,M) in enumerate(items):
  val,v,w=margin(M,seed=100+i)
  print(lab,'six',list(map(int,[M[0,0],M[1,1],M[2,2],M[1,2],M[0,2],M[0,1]])),'margin',None if val is None else round(val,10),'v',None if v is None else np.round(v,5),'walls',None if w is None else np.round(w,5))
