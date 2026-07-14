# -*- coding: utf-8 -*-
"""fR_spectral_fp.py — solve the GLOBAL f(R) fixed-point function by Chebyshev collocation
(the AS-literature method), on the pole-free Benedetti-Caravelli type-II equation (DM eqs
3.2-3.4, which DM endorse for admitting genuine global solutions). Domain [0,Rmax<singularity].
Validate: solution exists, ~R^2 at large R, g*lambda* in ~0.05-0.15."""
import numpy as np
from scipy.optimize import root
import sys; sys.path.insert(0,'continuum'); from cheb import map_to
pi = np.pi

Rmax = 5.0                     # below fixed singularity (R(R^4-54R^2-54)=0 at R~7.41) & pole-free
Nn = 40
D, R = map_to(0.0, Rmax, Nn)   # nodes R (descending from Rmax to 0), diff matrix
D2 = D@D; D3 = D@D@D

def resid(f):
    fp = D@f; fpp = D2@f; fppp = D3@f
    T2  = 40*(R*fpp - 4*fp)/((R-2)*fp - 2*f)
    T1  = -36.0
    T0n = -12 - 5*R**2
    T0h = ( R*(R**4 - 54*R**2 - 54)*fppp
            - (R**3 + 18*R**2 + 12)*(R*fpp - fp)
            + 18*(R**2 + 2)*(fp + 6*fpp) ) / (9*fpp - (R-3)*fp + 2*f)
    LHS = 768*pi**2*(2*f - R*fp)
    return LHS - (T2 + T1 + T0n + T0h)

# seed: small quadratic (Einstein-Hilbert + R^2 shape)
f0 = -0.005 - 0.023*R + 0.004*R**2
sol = root(resid, f0, method='lm', tol=1e-11)
f = sol.x
print(f"solver: success={sol.success}, |residual|_max = {np.max(np.abs(resid(f))):.2e}, nfev={sol.nfev}")
if not sol.success and np.max(np.abs(resid(f)))>1e-6:
    print("  (retry hybr)"); sol=root(resid,f0,method='hybr',tol=1e-11); f=sol.x
    print(f"  retry: |res|_max={np.max(np.abs(resid(f))):.2e}")

# report f* at R=0 and its low-order Taylor (fit) -> g0,g1,g2
# fit a polynomial to f near R=0 for the couplings
order=6
# interpolate f on nodes -> Chebyshev; evaluate f and derivatives at R=0 via D matrices
i0 = np.argmin(np.abs(R-0.0))            # node nearest 0 (should be exactly Rmax*... endpoint)
fp = D@f; fpp=D2@f
g0 = f[np.argmin(np.abs(R))]             # f(0)
# value/deriv at R=0: use the node at R=0 (endpoint of [0,Rmax] -> last node)
idx0 = np.argmin(np.abs(R))
g0 = f[idx0]; g1 = (D@f)[idx0]; g2 = 0.5*(D2@f)[idx0]
print(f"\nf*(0)=g0={g0:+.5f}   f*'(0)=g1={g1:+.5f}   f*''(0)/2=g2={g2:+.5f}")
if g1!=0:
    G_=1/(16*pi*g1); Lam=-g0/(2*g1); print(f"g*lambda* = {G_*Lam:+.4f}  (band 0.05-0.15)")
# large-R behaviour: f*(Rmax)/Rmax^2 should approach ~const if ~R^2
print(f"\nlarge-R check: f*(R)/R^2 at R={R[0]:.2f} : {f[0]/R[0]**2:+.4f}")
print(f"               f*(R)/R^2 at R={R[3]:.2f} : {f[3]/R[3]**2:+.4f}   (flat => R^2 asymptotics)")
np.save('continuum/fstar_spectral.npy', np.vstack([R,f]))
print("\nsaved f*(R) on grid -> continuum/fstar_spectral.npy")
