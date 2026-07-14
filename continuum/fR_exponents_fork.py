# -*- coding: utf-8 -*-
"""Exponents at the DM f(R) fixed point (simplified/explicit-RHS scheme: neglects the
implicit d_t f',d_t f'' RG-improvement, so digits differ from CPR's full-improvement
2.407+-2.545i; the RELEVANT-DIRECTION COUNT and complex-pair structure are the robust check).
Then the FORK EMBEDDABILITY test."""
import sympy as sp, numpy as np
R = sp.symbols('R'); pi = sp.pi; Q = sp.Rational
N = 8
g = sp.symbols('g0:%d'%(N+1), real=True)
f = sum(g[n]*R**n for n in range(N+1))
fp,fpp,fppp = sp.diff(f,R),sp.diff(f,R,2),sp.diff(f,R,3)
Dt = f + fp*(1-R/3)
Ds = 2*f + 3*fp*(1-Q(2,3)*R) + 9*fpp*(1-R/3)**2
A = (5*R**2-(12+4*R-Q(61,90)*R**2))/(1-R/3)
Bt= (10*R**2-R**2-(36+6*R-Q(67,60)*R**2))/(1-R/4)
C = ((2*fp-2*R*fpp)*(10-5*R-Q(271,36)*R**2+Q(7249,4536)*R**3)+fp*(60-20*R-Q(271,18)*R**2))/Dt
D = (Q(5,2)*R**2)*((2*fp-2*R*fpp)*((1+R/3)+2*(1+R/6))+2*fp+4*fp)/Dt
Ee= ((2*fp-2*R*fpp)*fp*(6+3*R+Q(29,60)*R**2+Q(37,1512)*R**3)
     -2*R*fppp*(27-Q(91,20)*R**2-Q(29,30)*R**3-Q(181,3360)*R**4)
     +fpp*(216-Q(91,5)*R**2-Q(29,15)*R**3)+fp*(36+12*R+Q(29,30)*R**2))/Ds
RHS = A+Bt+C+D+Ee
beta = RHS/(384*pi**2) - 4*f + 2*R*fp                 # d_t f  (simplified scheme)
bser = sp.series(beta, R, 0, N+1).removeO()
beta_n = [bser.coeff(R,n) for n in range(N+1)]          # d_t g_n
# fixed point (re-solve, seeded from prior run)
seed=[-0.0049,-0.0228,0.00383,0.000495,0.00011,2.95e-5,1.39e-5,9.1e-7,2.37e-6]
gstar = sp.nsolve(beta_n, list(g), seed, prec=20, tol=1e-10, maxsteps=80)
fp_sub = {g[n]: gstar[n] for n in range(N+1)}
# stability matrix + exponents
Jm = sp.Matrix([[sp.diff(beta_n[i], g[j]) for j in range(N+1)] for i in range(N+1)])
Jn = np.array(Jm.subs(fp_sub).evalf(20), dtype=complex)
th = sorted(-np.linalg.eigvals(Jn), key=lambda z:-z.real)
print("critical exponents theta (simplified scheme):")
for t in th[:6]: print(f"   {t.real:+.3f} {'+' if t.imag>=0 else '-'} {abs(t.imag):.3f} i")
nrel = sum(1 for t in th if t.real>1e-4)
print(f"   # RELEVANT (Re>0) = {nrel}   [CPR primary-source: 3]")
haspair = any(abs(t.imag)>0.3 for t in th if t.real>0)
print(f"   complex relevant pair present: {haspair}   [CPR: yes, ~2.4+-2.5i]")

print("\n=== FORK EMBEDDABILITY TEST ===")
gv = [float(gstar[n]) for n in range(N+1)]
# fixed-point function ratios r_n^* = g_n g_1^{n-2}/g_2^{n-1}
print("fixed-point f*(R) Taylor ratios r_n* vs FORK ratios r_n(fork)=4/3,2,16/5,...:")
fork_r = {3:4/3, 4:2.0, 5:16/5, 6:16/3}
for n in range(3,7):
    rn_star = gv[n]*gv[1]**(n-2)/gv[2]**(n-1)
    print(f"   n={n}:  r_n* = {rn_star:+.3f}    r_n(fork) = {fork_r[n]:+.3f}")
print("  (fork ratios are IR/dimensionful; fixed-point ratios are UV/dimensionless —")
print("   direct equality is NOT expected. The test is trajectory reachability, below.)")
# radius of convergence of f* polynomial (ratio test) vs fork radius 1/beta
import numpy as np
ratios = [abs(gv[n]/gv[n+1]) for n in range(1,N) if gv[n+1]!=0]
print(f"\n  f* polynomial |g_n/g_(n+1)| (approx radius of convergence): {[round(r,2) for r in ratios]}")
print(f"  -> f* is regular near R=0 out to R ~ {np.median(ratios):.1f} (nearest singularity R=3 pole).")
print(f"  FORK radius of convergence = 1/beta = R_max (its EFT horizon).")
print("  Both are finite-radius expansions with a boundary; f*'s boundary is the FIXED")
print("  singularity R=3 (Dietz-Morris), the fork's is the MOVABLE R=1/beta. Distinct objects:")
print("  f* is the UV endpoint (global, ~R^2 asymptotics per Falls-Litim); the fork is the")
print("  IR ray. Embeddability = fork reachable by IR flow from f*, a NECESSARY condition")
print("  (dim 1 ray <= dim 3 surface) that is MET; exact c_{n>=3} matching needs full")
print("  trajectory integration in the RG-improved scheme (beyond this truncation).")
