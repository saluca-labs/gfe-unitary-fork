# -*- coding: utf-8 -*-
"""lpa_flow.py — INTEGRATE a functional RG flow derived from scratch: scalar LPA, Litim
cutoff, d=3 Ising. Fixed point via expansion around the potential MINIMUM kappa (robust).
Validate nu against the known Litim-LPA value ~0.6496. Instrument calibration for gravity."""
import sympy as sp, numpy as np
sp.init_printing()

# Derivation: Wetterich LPA + Litim => d_t U = c_d k^{d+2}/(k^2+U''); rescale u=c_d v,
# phi~=c_d^{1/2}phi^ absorbs c_d (nu invariant). Dimensionless FP eqn (d=3), c_d->1:
#     0 = F(rho) := -3 u + rho u' + 1/(1 + u' + 2 rho u'')      ('=d/drho)
d, M = 3, 6
rho = sp.symbols('rho', positive=True)
kap = sp.symbols('kappa', positive=True)
u0  = sp.symbols('u0', real=True)
lam = sp.symbols('l2:%d'%(M+1), real=True)     # l2..lM = derivatives at the minimum
# potential expanded about the minimum kappa (no linear term => kappa is the minimum)
u = u0 + sum(lam[i]*(rho-kap)**(i+2)/sp.factorial(i+2) for i in range(len(lam)))
up  = sp.diff(u, rho); upp = sp.diff(u, rho, 2)
F = -d*u + rho*up + 1/(1 + up + 2*rho*upp)
# fixed point: F(rho)=0 for all rho => Taylor coeffs about rho=kappa, orders 0..M = 0
ser = sp.series(F, rho, kap, M+1).removeO()
eqs = [ser.coeff((rho-kap), n) for n in range(M+1)]
unknowns = [u0, kap] + list(lam)

def solve_from(seed):
    try:
        s = sp.nsolve(eqs, unknowns, seed, prec=18, tol=1e-12)
        return s
    except Exception:
        return None

fp = None
for kseed in [0.3,0.5,0.8,1.0,1.5]:
    for l2seed in [0.5,1.0,2.0]:
        seed = [0.1, kseed, l2seed] + [0.0]*(len(lam)-1)
        s = solve_from(seed)
        if s is not None and float(s[1])>1e-3 and abs(float(s[2]))>1e-3:
            fp = {unknowns[i]: s[i] for i in range(len(unknowns))}; 
            print(f"WF found: seed kappa={kseed}, l2={l2seed}")
            break
    if fp: break
if fp is None:
    print("!! WF not found"); raise SystemExit

print("Wilson-Fisher fixed point (LPA, d=3, Litim), minimum expansion:")
print(f"   kappa*  = {float(fp[kap]):.5f}   (dimensionless potential minimum)")
print(f"   u(kappa)= {float(fp[u0]):.5f}")
for i in range(len(lam)):
    print(f"   lambda_{i+2} = {float(fp[lam[i]]):+.5f}")

# stability matrix in the (u0,kappa,lam...) basis: theta = -eig(d eqs/d unknown)
J = sp.Matrix([[sp.diff(eqs[i], unknowns[j]) for j in range(len(unknowns))]
               for i in range(len(unknowns))])
Jn = np.array(J.subs(fp).evalf(20), dtype=complex)
th = sorted(-np.linalg.eigvals(Jn), key=lambda z:-z.real)
print("\nCritical exponents theta (Re>0 = relevant):")
for t in th: print(f"   {t.real:+.4f} {'+' if t.imag>=0 else '-'} {abs(t.imag):.4f} i")
theta1 = max(t.real for t in th)
nu = 1/theta1
print(f"\n   leading theta_1 = {theta1:.4f}  =>  nu = 1/theta_1 = {nu:.4f}")
print(f"   LITERATURE (Litim-LPA d=3 Ising): nu = 0.6496 ; exact ~ 0.6301")
print(f"   VALIDATION |nu-0.6496| = {abs(nu-0.6496):.4f}  =>",
      "PASS" if abs(nu-0.6496)<0.03 else "CHECK")
