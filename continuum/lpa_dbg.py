# -*- coding: utf-8 -*-
import sympy as sp, numpy as np
d, M = 3, 5
x, kap, u0 = sp.symbols('x kappa u0', real=True)
lam = sp.symbols('l2:%d'%(M+1), real=True)
# u as function of x=rho-kappa ; rho=kappa+x
u = u0 + sum(lam[i]*x**(i+2)/sp.factorial(i+2) for i in range(len(lam)))
ux = sp.diff(u, x); uxx = sp.diff(u, x, 2)   # d/drho = d/dx
F = -d*u + (kap + x)*ux + 1/(1 + ux + 2*(kap + x)*uxx)
ser = sp.series(F, x, 0, M+1).removeO()
P = sp.Poly(sp.expand(ser), x)
eqs = [P.coeff_monomial(x**n) for n in range(M+1)]
unknowns = [u0, kap] + list(lam)
print("n_eqs =", len(eqs), " n_unknowns =", len(unknowns))
print("eq0 (F at min) =", sp.simplify(eqs[0]))
print("eq1 =", sp.simplify(eqs[1]))
# hunt WF
fp=None
for ks in [0.2,0.4,0.6,0.9,1.2,1.8]:
  for ls in [0.3,0.8,1.5,3.0]:
    try:
      s=sp.nsolve(eqs, unknowns, [0.1,ks,ls]+[0.0]*(len(lam)-1), prec=18, tol=1e-11)
      if float(s[1])>1e-3 and abs(float(s[2]))>1e-3:
        fp={unknowns[i]:s[i] for i in range(len(unknowns))}
        print(f"WF found seed(kap={ks},l2={ls}): kappa*={float(s[1]):.4f} l2*={float(s[2]):.4f}")
        break
    except Exception as e:
      pass
  if fp: break
if not fp: print("no WF"); raise SystemExit
J=sp.Matrix([[sp.diff(eqs[i],unknowns[j]) for j in range(len(unknowns))] for i in range(len(unknowns))])
th=sorted(-np.linalg.eigvals(np.array(J.subs(fp).evalf(20),dtype=complex)),key=lambda z:-z.real)
print("thetas:", [f"{t.real:.3f}{'+' if t.imag>=0 else '-'}{abs(t.imag):.3f}i" for t in th])
nu=1/max(t.real for t in th)
print(f"nu = {nu:.4f}   (target Litim-LPA 0.6496)  =>", "PASS" if abs(nu-0.6496)<0.03 else "CHECK")
