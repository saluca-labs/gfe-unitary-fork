# -*- coding: utf-8 -*-
"""Robust WF finder: nsolve with realistic tol, wide seeds, dedupe roots, track the leading
NON-VOLUME exponent vs truncation order M. True WF => nu -> ~0.6496 as M grows."""
import sympy as sp, numpy as np
d = 3
def systemM(M):
    x, kap, u0 = sp.symbols('x kappa u0', real=True)
    lam = sp.symbols('l2:%d'%(M+1), real=True)
    u = u0 + sum(lam[i]*x**(i+2)/sp.factorial(i+2) for i in range(len(lam)))
    ux, uxx = sp.diff(u,x), sp.diff(u,x,2)
    F = -d*u + (kap+x)*ux + 1/(1 + ux + 2*(kap+x)*uxx)
    P = sp.Poly(sp.expand(sp.series(F,x,0,M+1).removeO()), x)
    eqs = [P.coeff_monomial(x**n) for n in range(M+1)]
    unk = [u0,kap]+list(lam)
    J = sp.Matrix([[sp.diff(eqs[i],unk[j]) for j in range(len(unk))] for i in range(len(unk))])
    return unk, eqs, J

for M in (4,5,6,7,8):
    unk, eqs, J = systemM(M); n=len(unk); found=[]
    for ks in [0.05,0.12,0.25,0.4,0.6,0.9,1.3]:
        for ls in [0.4,0.9,1.8,3.2]:
            try:
                s = sp.nsolve(eqs, unk, [0.05,ks,ls]+[0.0]*(n-2), tol=1e-9, maxsteps=80)
                kv,l2v = float(s[1]), float(s[2])
                if kv>1e-2 and l2v>1e-2 and not any(abs(kv-f)<1e-2 for f in found):
                    found.append(kv)
                    Jn=np.array(J.subs({unk[i]:s[i] for i in range(n)}).evalf(20),dtype=complex)
                    th=sorted(-np.linalg.eigvals(Jn), key=lambda z:-z.real)
                    tl=list(th); tl.pop(min(range(len(tl)),key=lambda i:abs(tl[i]-d)))  # drop volume mode
                    lead=max(t.real for t in tl)
                    nu=1/lead if lead>0 else 0
                    mark=" <== WF" if abs(nu-0.6496)<0.06 else ""
                    print(f"M={M}: kappa*={kv:6.4f} l2*={l2v:6.4f} theta1={lead:6.3f} nu={nu:6.4f}{mark}")
            except Exception:
                pass
