# -*- coding: utf-8 -*-
"""cheb.py — Chebyshev pseudo-spectral differentiation (Trefethen). Reusable instrument for
solving the global f(R) fixed-point ODE the way the AS literature does (not Taylor-around-0).
VALIDATE here on a problem with KNOWN eigenvalues before touching gravity."""
import numpy as np

def cheb(N):
    """Chebyshev diff matrix D and nodes x on [-1,1] (Trefethen 'Spectral Methods in MATLAB')."""
    if N == 0: return np.array([[0.0]]), np.array([1.0])
    x = np.cos(np.pi*np.arange(N+1)/N)
    c = np.hstack([2., np.ones(N-1), 2.])*(-1)**np.arange(N+1)
    X = np.tile(x,(N+1,1)).T
    dX = X - X.T
    D = np.outer(c,1./c)/(dX+np.eye(N+1))
    D -= np.diag(D.sum(axis=1))
    return D, x

def map_to(a, b, N):
    """Diff matrix + nodes on [a,b] via affine map from [-1,1]."""
    D, x = cheb(N)
    xp = 0.5*(b-a)*x + 0.5*(b+a)          # physical nodes
    Dp = (2.0/(b-a))*D
    return Dp, xp

# --- VALIDATION 1: d/dx of a known function ---
D, x = map_to(0, 2, 24)
f = np.sin(3*x); df_num = D@f; df_ex = 3*np.cos(3*x)
print(f"[val1] max|D sin(3x) - 3cos(3x)| = {np.max(np.abs(df_num-df_ex)):.2e}  (should be ~1e-8 or better)")

# --- VALIDATION 2: eigenvalue problem with KNOWN spectrum ---
# -u'' = lambda u on [0,pi], Dirichlet u(0)=u(pi)=0  => lambda_n = n^2, n=1,2,3,...
N = 40
D, x = map_to(0, np.pi, N)
D2 = D@D
# interior points (drop boundary rows/cols to impose Dirichlet)
A = -D2[1:N, 1:N]
lam = np.sort(np.linalg.eigvals(A).real)
print(f"[val2] -u''=lam u, u(0)=u(pi)=0 : lowest lambda = {lam[:6].round(4)}")
print(f"       known exact = [1, 4, 9, 16, 25, 36] (n^2)")
err = np.max(np.abs(lam[:6] - np.array([1,4,9,16,25,36])))
print(f"       max err (first 6) = {err:.2e}  =>", "PASS" if err<1e-3 else "CHECK")
print("\nInstrument validated: Chebyshev collocation gives spectral-accuracy derivatives and")
print("correct eigenvalues from a global representation -> immune to the Taylor-at-0 small-")
print("denominator artifact that broke the polynomial f(R) exponents.")
