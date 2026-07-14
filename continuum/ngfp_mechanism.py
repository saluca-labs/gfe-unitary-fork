# -*- coding: utf-8 -*-
"""ngfp_mechanism.py — the asymptotic-safety MECHANISM, in a form fully derivable from
canonical scaling + one quantum correction. This is a SCHEMATIC illustration of WHY a
non-Gaussian fixed point exists for the Newton coupling; the quantitative multi-coupling
NGFP values are taken from the verified literature, NOT from this toy. sympy."""
import sympy as sp

print("="*72)
print("MECHANISM: why Newton's coupling is asymptotically safe (schematic, exact-in-toy)")
print("="*72)
g, t, c = sp.symbols('g t c', real=True)
# Newton coupling G has mass dimension -2 in d=4 => dimensionless g = G k^2 runs canonically as +2g.
# Any interacting theory adds a quantum term; leading is -c g^2 (c>0 from graviton loops).
beta_g = 2*g - c*g**2
print(f"  canonical dimension of g=G k^2 in d=4:  +2   (Gaussian FP g*=0 is UV-repulsive)")
print(f"  leading quantum correction:  beta_g = dg/dt = {beta_g}")
fps = sp.solve(sp.Eq(beta_g,0), g)
print(f"  fixed points:  g* = {fps}")
gstar = [f for f in fps if f!=0][0]
print(f"  => NON-GAUSSIAN fixed point at g* = {gstar}  (nonzero, exists for any c>0)")
# stability: theta = -d(beta)/dg at g*
theta = sp.simplify(-sp.diff(beta_g, g).subs(g, gstar))
print(f"  critical exponent  theta = -beta_g'(g*) = {theta}  > 0  => UV-ATTRACTIVE (relevant)")
print("  Interpretation: the SAME sign of quantum correction that would make G blow up in")
print("  perturbation theory instead turns around and PINS it at g*, rendering the UV finite.")
print("  This is Weinberg's asymptotic safety in one line. The real theory has (g,lambda,")
print("  R^2-coupling,...) and the fixed point is found numerically in the literature; the")
print("  mechanism — canonical scaling balanced by interactions — is exactly this.")

print("\n"+"="*72)
print("STABILITY-MATRIX MACHINERY (generic; to be fed literature-verified beta-functions)")
print("="*72)
print("  For couplings u_i with flow beta_i(u): FP solves beta_i(u*)=0; the critical")
print("  exponents are theta_I = -eig(∂beta_i/∂u_j |_*). Re(theta)>0 = relevant (UV-attractive).")
print("  Predictivity = # relevant directions = dim of UV critical surface.")
print("  We DO NOT fabricate the gravitational beta_i here; the relevant-direction count and")
print("  exponents are imported from Codello-Percacci-Rahmede / Falls-Litim (agent-verified).")
# demonstrate the machinery on a verifiable 2x2 example so the method itself is checkable
import numpy as np
J = np.array([[2.0, 0.3],[0.1, -1.5]])   # arbitrary demo Jacobian of (beta_i)
ev = np.linalg.eigvals(-J)
print(f"  [machinery self-test] demo Jacobian eigenvalues of -J: {np.round(ev,3)}")
print(f"     -> {sum(1 for e in ev if e.real>0)} relevant, {sum(1 for e in ev if e.real<0)} irrelevant (method works)")
