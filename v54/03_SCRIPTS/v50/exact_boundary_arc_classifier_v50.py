#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp, sys, json
from fractions import Fraction
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from selling_global_wall_linearization_v50 import R, w0, D0, labs
q,r=sp.symbols('q r', real=True)

def primitive_factor_list(poly):
    P=sp.Poly(sp.expand(poly),q,r,domain=sp.QQ)
    coeff,fac=sp.factor_list(P.as_expr())
    out=[]
    for f,m in fac:
        F=sp.Poly(f,q,r,domain=sp.QQ)
        # normalize primitive with positive leading coefficient under lex tuple
        cont,prim=sp.polys.polytools.primitive(F.as_expr(),q,r)
        F=sp.Poly(prim,q,r,domain=sp.QQ)
        if F.LC()<0:F=-F
        out.append((sp.factor(F.as_expr()),m))
    return out

def divisible(poly,f):
    P=sp.Poly(poly,q,r,domain=sp.QQ);F=sp.Poly(f,q,r,domain=sp.QQ)
    return P.rem(F).is_zero

def critical_q_poly(f,polys):
    terms=[]
    F=sp.Poly(f,q,r,domain=sp.QQ)
    # degree drop / vertical tangency
    lc=sp.Poly(F.as_poly(r).LC(),q,domain=sp.QQ)
    if lc.degree()>0:terms.append(lc.as_expr())
    dr=sp.diff(f,r)
    if dr!=0:
        res=sp.resultant(f,dr,r)
        if res!=0:terms.append(sp.Poly(res,q,domain=sp.QQ).as_expr())
    for g in polys:
        if divisible(g,f):continue
        res=sp.resultant(f,g,r)
        if res!=0:
            rr=sp.Poly(res,q,domain=sp.QQ)
            if rr.degree()>0:terms.append(rr.as_expr())
    if not terms:return sp.Poly(1,q,domain=sp.QQ)
    prod=sp.Poly(sp.prod(terms),q,domain=sp.QQ).sqf_part()
    return prod

def real_root_intervals(P,eps=sp.Rational(1,10)**10):
    if P.degree()<=0:return []
    return [iv for iv,m in sp.intervals(P,eps=eps) if m>=1]

def q_samples_from_critical(P):
    ints=real_root_intervals(P)
    if not ints:return [sp.Rational(0)]
    samples=[]
    # outer left/right plus gaps. use exact rationals strictly outside isolating intervals
    a0,b0=ints[0]; samples.append(sp.floor(a0)-1)
    for (a,b),(c,d) in zip(ints,ints[1:]):
        if b<c:samples.append((b+c)/2)
    an,bn=ints[-1]; samples.append(sp.ceiling(bn)+1)
    return samples

def root_sign(poly_r, root_poly, iv):
    """Exact sign of poly_r at a real root of root_poly isolated by rational iv."""
    G=sp.Poly(poly_r,r,domain=sp.QQ)
    F=sp.Poly(root_poly,r,domain=sp.QQ)
    if G.is_zero:return 0
    gd=sp.gcd(F,G)
    if gd.degree()>0:
        # if the isolated F-root is also a gcd root, sign zero
        for giv,m in sp.intervals(gd,eps=sp.Rational(1,10)**14):
            if not (giv[1] < iv[0] or iv[1] < giv[0]):return 0
    # isolate roots of G and refine F interval until disjoint from them
    groots=[x for x,m in sp.intervals(G,eps=sp.Rational(1,10)**14)] if G.degree()>0 else []
    a,b=iv
    for _ in range(12):
        overlap=any(not (gb<a or b<ga) for ga,gb in groots)
        if not overlap:break
        a,b=sp.refine_root(F,a,b,eps=(b-a)/10)
    mid=(a+b)/2
    val=G.eval(mid)
    if val>0:return 1
    if val<0:return -1
    # Should only happen at exact common root, caught above.
    raise RuntimeError(('undecided-sign',G,F,(a,b)))

def classify_vertical_factor(f, polys, wallset):
    F=sp.Poly(f,q,r,domain=sp.QQ)
    if F.degree(r)!=0 or F.degree(q)!=1:
        return []
    qq=sp.solve(sp.Eq(f,0),q)[0]
    crit=[]
    for p in list(polys.values())+[D0]:
        pr=sp.Poly(sp.expand(p.subs(q,qq)),r,domain=sp.QQ)
        if pr.is_zero or pr.degree()<=0: continue
        crit.append(pr.as_expr())
    if crit:
        C=sp.Poly(sp.prod(crit),r,domain=sp.QQ).sqf_part()
        ints=[iv for iv,m in sp.intervals(C,eps=sp.Rational(1,10)**12)]
    else: ints=[]
    samples=[]
    if ints:
        samples.append(sp.floor(ints[0][0])-1)
        for (a,b),(c,d) in zip(ints,ints[1:]):
            if b<c:samples.append((b+c)/2)
        samples.append(sp.ceiling(ints[-1][1])+1)
    else:samples=[sp.Rational(0)]
    valid=[]
    for rr in samples:
        signs={}
        for lab,p in polys.items():
            val=sp.simplify(p.subs({q:qq,r:rr}))
            signs[lab]=0 if val==0 else (1 if val>0 else -1)
        dv=sp.simplify(D0.subs({q:qq,r:rr}))
        if dv>0 and all(v<=0 for v in signs.values()):
            valid.append({'q':str(qq),'r':str(rr),'signs':signs})
    return valid

def classify_cell(L):
    vec=L*w0
    polys={labs[i]:sp.factor(vec[i]) for i in range(6)}
    unique={}
    for lab,p in polys.items():
        for f,m in primitive_factor_list(p):
            key=sp.srepr(f)
            unique.setdefault(key,{'factor':f,'walls':set()})['walls'].add(lab)
    arcs=[]
    all_constraints=list(polys.values())+[D0]
    for item in unique.values():
        f=item['factor']; wallset=sorted(item['walls'])
        vvalid=classify_vertical_factor(f,polys,wallset)
        if vvalid:
            arcs.append({'factor':str(f),'walls':wallset,'kind':'ORDINARY' if len(wallset)==1 else 'MULTIWALL','samples':vvalid})
            continue
        C=critical_q_poly(f,all_constraints)
        qs=q_samples_from_critical(C)
        valid=[]
        for q0 in qs:
            fq=sp.Poly(sp.expand(f.subs(q,q0)),r,domain=sp.QQ)
            if fq.is_zero or fq.degree()<=0:continue
            for iv,mult in sp.intervals(fq,eps=sp.Rational(1,10)**12):
                # exact signs at this algebraic root
                signs={}
                for lab,p in polys.items():
                    signs[lab]=root_sign(sp.expand(p.subs(q,q0)),fq,iv)
                sd=root_sign(sp.expand(D0.subs(q,q0)),fq,iv)
                if sd>0 and all(v<=0 for v in signs.values()):
                    valid.append({'q':str(q0),'r_interval':[str(iv[0]),str(iv[1])],'signs':signs})
        if valid:
            arcs.append({'factor':str(f),'walls':wallset,'kind':'ORDINARY' if len(wallset)==1 else 'MULTIWALL','samples':valid})
    ordinary=sorted({a['walls'][0] for a in arcs if a['kind']=='ORDINARY'})
    multi=[a for a in arcs if a['kind']=='MULTIWALL']
    return {'polynomials':{k:str(v) for k,v in polys.items()},'arcs':arcs,'ordinary_walls':ordinary,'multiwall_components':multi}

def main():
    # exact regression cells
    cells={'Q':sp.eye(6),'H':R['h'],'M':R['m'],'G':R['g'],
           'X':R['g']*R['h'],'Y':R['g']*R['m'],'Z':R['h']*R['g'],'W':R['m']*R['g']}
    out={}
    for name,L in cells.items():
        c=classify_cell(L);out[name]=c
        print(name,'ordinary',c['ordinary_walls'],'multi',[(x['walls'],x['factor']) for x in c['multiwall_components']])
    (HERE/'exact_boundary_arc_classifier_regression_v50.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
if __name__=='__main__':main()
