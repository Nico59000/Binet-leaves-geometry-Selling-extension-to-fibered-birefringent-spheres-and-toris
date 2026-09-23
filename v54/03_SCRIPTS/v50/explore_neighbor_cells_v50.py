import sympy as sp, sys, math, json
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0,'/mnt/data/v50_work')
from selling_reconstruction_engine_v49 import gram_from_six,S,six_from_gram

# wall polynomials H_w(v)=homogeneous coefficient of positive correspondent 2vvT-M
W=['g','h','k','l','m','n']
def wall_exprs(M,v):
    x,y,z=v; a,b,c=M[0,0],M[1,1],M[2,2]; k,h,g=M[0,1],M[0,2],M[1,2]; ss=x+y+z
    return {
      'g':sp.expand(2*y*z-g), 'h':sp.expand(2*x*z-h), 'k':sp.expand(2*x*y-k),
      'l':sp.expand(a+h+k-2*x*ss), 'm':sp.expand(b+k+g-2*y*ss), 'n':sp.expand(c+h+g-2*z*ss)}

def proj_polys(M,chart=0):
    u,v=sp.symbols('u v', real=True)
    d=[None,None,None]; d[chart]=sp.Integer(1); rest=[i for i in range(3) if i!=chart]; d[rest[0]]=u;d[rest[1]]=v
    d=sp.Matrix(d); Minv=M.inv(); phi=sp.factor((d.T*Minv*d)[0])
    # if phi>0, s^2=1/phi; H(s d) = quad/phi - const => numerator H*phi
    # substitute scale symbol t and derive coefficient t^2 and constant
    x,y,z=sp.symbols('x y z'); t=sp.symbols('t')
    ww=wall_exprs(M,(x,y,z)); subs={x:t*d[0],y:t*d[1],z:t*d[2]}
    out={}
    for k,e in ww.items():
      et=sp.expand(e.subs(subs)); A=sp.expand(et.coeff(t,2)); C=sp.expand(et.subs(t,0))
      # at t^2=1/phi: (A+C*phi)/phi
      out[k]=sp.factor(A+C*phi)
    # clear denom common det if needed
    den=sp.lcm([sp.denom(e) for e in [phi,*out.values()]])
    phi=sp.factor(phi*den); out={k:sp.factor(e*den) for k,e in out.items()}
    return u,v,phi,out

def eval_region(phi, walls, uv):
  u,v=uv; fphi=sp.lambdify(('u','v'),phi,'numpy'); fw={k:sp.lambdify(('u','v'),e,'numpy') for k,e in walls.items()}
  ph=float(fphi(u,v)); vals={k:float(f(u,v)) for k,f in fw.items()}
  return ph,vals, ph>0 and all(x<=1e-9 for x in vals.values())

M0=gram_from_six(7,-1,-8,0,2,0)
qA=float(-4+2*sp.sqrt(210)/7); rB=float(1-sp.sqrt(195)/15)
qD=0.136218563230169494; rD=0.0634917291636040015
# outside directions near edge midpoint in original d=(1,-r,q)
def h_point(out=False):
 q=(qD+qA)/2; r=math.sqrt(max(0,(8-56*q-7*q*q)/60));
 if out: r*=1.01 # Hh decreases -> outside
 return np.array([1,-r,q],float)
def m_point(out=False):
 q=qD/2
 # solve Mm=0 lower root
 co=[60,(-120-120*q),(8+4*q-7*q*q)]
 roots=np.roots(co); r=min(x.real for x in roots if abs(x.imag)<1e-8 and x.real>=0)
 if out: r*=1.01 # need check Mm? likely lower branch domain Mm>=0, larger r enters negative initially
 return np.array([1,-r,q],float)
def g_point(out=False):
 r=rB/2; q=-1e-3 if out else 0
 return np.array([1,-r,q],float)

for wall,idx,ptfun in [('h',2,h_point),('m',5,m_point),('g',9,g_point)]:
 A=S[idx-1]; M=A.T*M0*A
 print('\n###',wall,'S',idx,'M',six_from_gram(M))
 for label in ['boundary','outside']:
  d=sp.Matrix(ptfun(label=='outside')); dp=A.T*d
  vals=np.array(dp,dtype=float).reshape(3)
  ch=int(np.argmax(np.abs(vals)))
  norm=vals/vals[ch]
  rest=[i for i in range(3) if i!=ch]
  print(label,'dp',vals,'chart',ch,'uv',norm[rest])
 u,v,phi,ww=proj_polys(M,ch)
 print('chart',ch,'phi=',phi)
 for k,e in ww.items(): print(k,e)
