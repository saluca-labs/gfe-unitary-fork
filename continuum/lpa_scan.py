# -*- coding: utf-8 -*-
"""Find the TRUE Wilson-Fisher FP: wide seed scan, collect distinct roots, for each report
the leading NON-VOLUME exponent (the volume/vacuum mode is always +d=3). WF => nu~0.65."""
import sympy as sp, numpy as np
d = 3
def build(M):
    x, kap, u0 = sp.symbols('x kappa u0', real=True)
    lam = sp.symbols('l2:%d'%(M+1), real=True)
    u = u0 + sum(lam[i]*x**(i+2)/sp.factorial(i+2) for i in range(len(lam)))
    ux, uxx = sp.diff(u,x), sp.diff(u,x,2)
    F = -d*u + (kap+x)*ux + 1/(1 + ux + 2*(kap+x)*uxx)
    P = sp.Poly(sp.expand(sp.series(F,x,0,M+1).removeO()), x)
    eqs = [P.coeff_monomial(x**n) for n in range(M+1)]
    unk = [u0,kap]+list(lam)
    J = sp.Matrix([[sp.diff(eqs[i],unk[j]) for j in range(len(unk))] for i in range(len(unk))])
    fe = [sp.lambdify(unk, e, 'mpmath') for e in eqs]
    return unk, eqs, J

for M in (5,6,7):
    unk, eqs, J = build(M)
    n = len(unk)
    roots = []
    import itertools
    for ks in [0.05,0.1,0.2,0.35,0.5,0.7,1.0,1.4]:
        for ls in [0.3,0.6,1.0,1.6,2.5,4.0]:
            try:
                s = sp.nsolve(eqs, unk, [0.05,ks,ls]+[0.0]*(n-2), prec=20, tol=1e-12)
                kv, l2v = float(s[1]), float(s[2])
                if kv>1e-3 and abs(l2v)>1e-3:
                    key = (round(kv,3), round(l2v,3))
                    if key not in [r[0] for r in roots]:
                        Jn = np.array(J.subs({unk[i]:s[i] for i in range(n)}).evalf(20),dtype=complex)
                        th = sorted(-np.linalg.eigvals(Jn), key=lambda z:-z.real)
                        # non-volume: drop one eigenvalue closest to +d=3
                        thl = list(th); 
                        vidx = min(range(len(thl)), key=lambda i: abs(thl[i]-d))
                        thl.pop(vidx)
                        lead = max(t.real for t in thl)
                        roots.append((key, kv, l2v, 1/lead if lead>0 else None, lead))
            except Exception:
                pass
    print(f"--- M={M}: {len(roots)} distinct nontrivial roots ---")
    for key,kv,l2v,nu,lead in sorted(roots, key=lambda r: (abs((r[3] or 9)-0.6496))):
        tag = "  <== WF (nu~0.65)" if nu and abs(nu-0.6496)<0.05 else ""
        print(f"   kappa*={kv:7.4f}  l2*={l2v:7.4f}  theta_lead={lead:6.3f}  nu={nu if nu else 0:6.4f}{tag}")
