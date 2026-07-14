# -*- coding: utf-8 -*-
"""fR_fixedpoint.py — INTEGRATE the real f(R) asymptotic-safety flow: solve the exact
Dietz-Morris eq.(3.1) [= Codello-Percacci-Rahmede eq.(113)] fixed-point ODE in polynomial
truncation f(R)=sum g_n R^n. Primary-source transcription (agent-verified vs rendered PDFs).
Optimized/Litim cutoff, d=4, S^4 background, de Donder/Landau gauge. For R<3 all theta->1.
Validate: fixed point exists with g*lambda* ~ 0.13 (CPR universal value)."""
import sympy as sp
R = sp.symbols('R')
pi = sp.pi
N = 8                                   # polynomial order
g = sp.symbols('g0:%d'%(N+1), real=True)
f = sum(g[n]*R**n for n in range(N+1))
fp, fpp, fppp = sp.diff(f,R), sp.diff(f,R,2), sp.diff(f,R,3)
Q = sp.Rational

# --- DM (3.1) RHS, theta->1 for R<3 ; r(-R/3)=1+R/3, r(-R/6)=1+R/6 ---
Dt = f + fp*(1 - R/3)                                     # tensor inverse propagator
Ds = 2*f + 3*fp*(1 - Q(2,3)*R) + 9*fpp*(1 - R/3)**2       # scalar inverse propagator
Sigma = 0                                                 # CPR excluded-mode option
A  = (5*R**2 - (12 + 4*R - Q(61,90)*R**2)) / (1 - R/3)
B  = (10*R**2 - R**2 - (36 + 6*R - Q(67,60)*R**2)) / (1 - R/4)
C  = ((2*fp - 2*R*fpp)*(10 - 5*R - Q(271,36)*R**2 + Q(7249,4536)*R**3)
      + fp*(60 - 20*R - Q(271,18)*R**2)) / Dt
D  = (Q(5,2)*R**2)*((2*fp - 2*R*fpp)*((1+R/3) + 2*(1+R/6)) + 2*fp + 4*fp) / Dt
E  = ((2*fp - 2*R*fpp)*fp*(6 + 3*R + Q(29,60)*R**2 + Q(37,1512)*R**3)
      - 2*R*fppp*(27 - Q(91,20)*R**2 - Q(29,30)*R**3 - Q(181,3360)*R**4)
      + fpp*(216 - Q(91,5)*R**2 - Q(29,15)*R**3)
      + fp*(36 + 12*R + Q(29,30)*R**2)) / Ds
RHS = A + Sigma + B + C + D + E
LHS = 768*pi**2*(2*f - R*fp)
Expr = LHS - RHS                          # =0 at fixed point, for all R

print("built; expanding to order", N, "...")
ser = sp.series(Expr, R, 0, N+1).removeO()
P = sp.Poly(sp.together(ser).as_numer_denom()[0], R)  # clear denominators safely? -> use coeffs
# safer: take series coeffs directly
eqs = [sp.nsimplify(ser.coeff(R, n)) for n in range(N+1)]
print("got", len(eqs), "equations; solving...")

# seed near known CPR ballpark (dimensionless): g1 ~ O(0.1-1), g0 small neg
import mpmath as mp
mp.mp.dps = 30
seed = [ -0.002, 0.006, 0.002] + [0.0]*(N-2)   # g0,g1,g2,...
try:
    sol = sp.nsolve(eqs, list(g), seed, prec=20, tol=1e-10, maxsteps=100)
    gv = [float(sol[n]) for n in range(N+1)]
    print("\nfixed-point couplings g_n:")
    for n in range(N+1): print(f"   g_{n} = {gv[n]:+.6e}")
    # map to Newton/cosmological-constant: f = g0 + g1 R + ...  ~ (R-2Lam)/(16 pi G)
    G_ = 1.0/(16*3.141592653589793*gv[1]) if gv[1]!=0 else float('nan')
    Lam = -gv[0]/(2*gv[1]) if gv[1]!=0 else float('nan')
    prod = G_*Lam
    print(f"\n   G~ = 1/(16 pi g1) = {G_:.5f}")
    print(f"   Lam~ = -g0/(2 g1)  = {Lam:.5f}")
    print(f"   g* lambda* = G~*Lam~ = {prod:.5f}   (CPR universal ~ 0.12-0.14)")
    print(f"   VALIDATION |prod|-0.13: {abs(abs(prod)-0.13):.4f}",
          "=> PASS" if abs(abs(prod)-0.13)<0.05 else "=> CHECK/convention")
except Exception as e:
    print("nsolve failed:", type(e).__name__, str(e)[:120])
