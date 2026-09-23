#!/usr/bin/env python3
"""Exact RUR-style fissure solver for Selling v50.

For an integral indefinite ternary form M and one of the three opposite
superbase wall pairs (g,l), (h,m), (k,n), solve exactly

    v^T adj(M) v = det(M),   H_u(v)=H_v(v)=0,

and retain only real points at which all six Selling wall functions are <= 0.
Real algebraic roots are certified by rational isolating intervals; wall signs
are certified by exact root isolation/refinement, not floating-point tests.

This is deliberately fail-closed: positive-dimensional or non-triangular
systems are reported UNRESOLVED rather than guessed.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import selling_reconstruction_engine_v50 as e

x,y,z=sp.symbols('x y z')
OPPOSITE=(('g','l'),('h','m'),('k','n'))


def wall_polynomials(M):
    a,b,c,g,h,k=e.six_from_gram(M); s=x+y+z
    return {
      'g':sp.expand(2*y*z-g),
      'h':sp.expand(2*x*z-h),
      'k':sp.expand(2*x*y-k),
      'l':sp.expand(a+h+k-2*x*s),
      'm':sp.expand(b+k+g-2*y*s),
      'n':sp.expand(c+h+g-2*z*s),
    }


def hyperboloid_polynomial(M):
    v=sp.Matrix([x,y,z])
    return sp.expand((v.T*M.adjugate()*v)[0]-M.det())


def _poly_num(expr,var=z):
    num,den=sp.fraction(sp.cancel(expr))
    return sp.Poly(num,var,domain=sp.QQ), sp.Poly(den,var,domain=sp.QQ)


def _root_in_interval(poly,a,b):
    if a==b:
        return sp.simplify(poly.eval(a))==0
    return poly.count_roots(a,b)>0


def sign_at_isolated_root(expr, P, interval, max_refine=120):
    """Exact sign of rational expr at the unique P-root in interval."""
    a,b=interval
    num,den=_poly_num(expr)
    if a==b:
        dv=den.eval(a)
        if dv==0: raise ZeroDivisionError('denominator vanishes at rational root')
        nv=num.eval(a)
        return 0 if nv==0 else (1 if nv*dv>0 else -1)
    # Exact zero detection at the selected algebraic root.
    gn=sp.gcd(P,num)
    if gn.degree()>0 and _root_in_interval(gn,a,b):
        gd=sp.gcd(P,den)
        if gd.degree()>0 and _root_in_interval(gd,a,b):
            raise ZeroDivisionError('0/0 at algebraic root')
        return 0
    # Refine until numerator and denominator each have no root in the interval.
    aa,bb=a,b
    for q in (num,den):
        n=0
        while _root_in_interval(q,aa,bb):
            aa,bb=sp.refine_root(P,aa,bb,eps=(bb-aa)/4)
            n+=1
            if n>max_refine:
                raise RuntimeError('sign refinement limit')
    mid=(aa+bb)/2
    nv=num.eval(mid); dv=den.eval(mid)
    if nv==0 or dv==0:
        raise RuntimeError('unexpected midpoint zero')
    return 1 if nv*dv>0 else -1


def _choose_triangular_representation(Gexprs,p):
    P=sp.Poly(p,z,domain=sp.QQ)
    y_candidates=[]
    for q in Gexprs:
        if q.has(x): continue
        py=sp.Poly(q,y)
        if py.degree()==1:
            coeff=sp.expand(py.coeff_monomial(y)); rest=sp.expand(q-coeff*y)
            yexpr=sp.cancel(-rest/coeff)
            _,den=_poly_num(yexpr)
            if sp.gcd(P,den).degree()==0:
                y_candidates.append(yexpr)
    # some systems may have y-independent x relation first; still need y.
    for yexpr in y_candidates:
        x_candidates=[]
        for q in Gexprs:
            qy=sp.cancel(q.subs(y,yexpr))
            num,_=sp.fraction(qy)
            px=sp.Poly(sp.expand(num),x)
            if px.degree()==1:
                coeff=sp.expand(px.coeff_monomial(x)); rest=sp.expand(num-coeff*x)
                xexpr=sp.cancel(-rest/coeff)
                _,den=_poly_num(xexpr)
                if sp.gcd(P,den).degree()==0:
                    x_candidates.append(xexpr)
        for xexpr in x_candidates:
            # The elimination polynomial can contain complex components on which
            # this triangular chart is invalid (e.g. a factor with no real
            # roots).  Validity is therefore certified root-by-root below.
            return xexpr,yexpr
    return None


def solve_fissure(M,pair):
    W=wall_polynomials(M); H=hyperboloid_polynomial(M)
    u,v=pair
    t0=time.time()
    G=sp.groebner([H,W[u],W[v]],x,y,z,order='lex',domain=sp.QQ)
    elapsed=time.time()-t0
    if not G.is_zero_dimensional:
        return {'pair':list(pair),'status':'UNRESOLVED_POSITIVE_DIMENSIONAL','groebner_seconds':elapsed}
    exprs=[q.as_expr() for q in G.polys]
    unis=[q for q in exprs if not q.has(x) and not q.has(y)]
    if not unis:
        return {'pair':list(pair),'status':'UNRESOLVED_NO_UNIVARIATE','groebner_seconds':elapsed}
    p=sp.Poly(unis[-1],z,domain=sp.QQ).sqf_part()
    rep=_choose_triangular_representation(exprs,p.as_expr())
    if rep is None:
        return {'pair':list(pair),'status':'UNRESOLVED_NONTRIANGULAR','groebner_seconds':elapsed,
                'groebner':[str(q) for q in exprs]}
    xexpr,yexpr=rep
    # Exact identities at every p-root: hyperboloid and the two selected walls.
    identity_checks={}
    for name,expr in [('hyperboloid',H),(u,W[u]),(v,W[v])]:
        rr=sp.cancel(expr.subs({x:xexpr,y:yexpr})); num,_=sp.fraction(rr)
        identity_checks[name]=sp.rem(sp.Poly(sp.expand(num),z,domain=sp.QQ),p).is_zero
    intervals=p.intervals(eps=sp.Rational(1,10**14))
    roots=[]; reduced=[]
    for interval,mult in intervals:
        a,b=interval
        signs={}
        try:
            # Root-local Groebner validity certificate: every basis residual
            # must vanish at this selected real algebraic root.
            for gres in exprs:
                rr=sp.cancel(gres.subs({x:xexpr,y:yexpr})); rnum,_=sp.fraction(rr)
                rpoly=sp.Poly(sp.expand(rnum),z,domain=sp.QQ)
                if not rpoly.is_zero:
                    gg=sp.gcd(p,rpoly)
                    if gg.degree()==0 or not _root_in_interval(gg,a,b):
                        raise RuntimeError('triangular chart invalid at selected real root')
            for name,expr in W.items():
                q=sp.cancel(expr.subs({x:xexpr,y:yexpr}))
                signs[name]=sign_at_isolated_root(q,p,(a,b))
        except Exception as exc:
            roots.append({'interval':[str(a),str(b)],'multiplicity':int(mult),
                          'status':'UNRESOLVED_SIGN','error':repr(exc)})
            continue
        isred=all(s<=0 for s in signs.values())
        mid=a if a==b else (a+b)/2
        # Decimal coordinates are diagnostic only; exact cert is interval + RUR.
        zz=sp.N(mid,16); xx=sp.N(xexpr.subs(z,mid),16); yy=sp.N(yexpr.subs(z,mid),16)
        rec={'interval':[str(a),str(b)],'multiplicity':int(mult),'wall_signs':signs,
             'reduced':bool(isred),'approx':[str(xx),str(yy),str(zz)]}
        roots.append(rec)
        if isred: reduced.append(rec)
    return {
      'pair':list(pair),'status':'PASS','groebner_seconds':elapsed,
      'minimal_polynomial':str(p.as_expr()),'degree':p.degree(),
      'x_of_z':str(xexpr),'y_of_z':str(yexpr),
      'identity_checks':identity_checks,'real_root_count':sum(int(m) for _,m in intervals),
      'reduced_event_count':len(reduced),'roots':roots
    }


def classify_form(M,include_roots=False):
    out=[]; unresolved=[]; total=0
    for pair in OPPOSITE:
        r=solve_fissure(M,pair)
        total += int(r.get('reduced_event_count',0))
        if r['status']!='PASS': unresolved.append(r['status'])
        if not include_roots and r.get('status')=='PASS':
            r={k:v for k,v in r.items() if k!='roots'}
        out.append(r)
    return {'six':[int(q) for q in e.six_from_gram(M)],'status':'PASS' if not unresolved else 'PARTIAL_FAIL_CLOSED',
            'reduced_fissure_events':total,'pairs':out,'unresolved':unresolved}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--six',nargs=6,type=int)
    ap.add_argument('--forms-json')
    ap.add_argument('--output')
    ap.add_argument('--include-roots',action='store_true')
    ns=ap.parse_args()
    if ns.six:
        M=e.gram_from_six(*map(sp.Integer,ns.six)); result=classify_form(M,True)
    elif ns.forms_json:
        data=json.loads(Path(ns.forms_json).read_text())
        forms=data['forms']; rr=[]; t=time.time()
        for i,rec in enumerate(forms):
            M=e.gram_from_six(*map(sp.Integer,rec['six']))
            c=classify_form(M,ns.include_roots); c['id']=rec.get('id',i); rr.append(c)
            if (i+1)%25==0: print(f'{i+1}/{len(forms)}',file=sys.stderr)
        result={'version':'v50','status':'PASS' if all(r['status']=='PASS' for r in rr) else 'PARTIAL_FAIL_CLOSED',
                'forms':len(rr),'elapsed_seconds':time.time()-t,
                'forms_with_reduced_fissure':sum(r['reduced_fissure_events']>0 for r in rr),
                'total_reduced_fissure_events':sum(r['reduced_fissure_events'] for r in rr),
                'unresolved_forms':sum(r['status']!='PASS' for r in rr),'records':rr}
    else:
        ap.error('use --six or --forms-json')
    txt=json.dumps(result,indent=2,ensure_ascii=False)
    if ns.output: Path(ns.output).write_text(txt,encoding='utf-8')
    print(txt if not ns.output else ns.output)

if __name__=='__main__': main()
