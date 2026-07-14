# -*- coding: utf-8 -*-
"""fR_spectral_v2.py — global f(R) FP by collocation on [0, R_c], R_c=sqrt(27+3 sqrt87)=7.415,
the exact zero of the f''' coefficient R(R^4-54R^2-54). At BOTH endpoints the equation
degenerates (f''' drops) -> natural regularity BCs -> well-posed. Newton from an R^2 seed."""
import numpy as np
from scipy.optimize import root
import sys; sys.path.insert(0,'continuum'); from cheb import map_to
pi = np.pi
Rc = np.sqrt(27 + 3*np.sqrt(87))          # = 7.4149...  exact root of R^4-54R^2-54
print(f"R_c = sqrt(27+3sqrt87) = {Rc:.5f}  (verify R^4-54R^2-54 = {Rc**4-54*Rc**2-54:.2e})")
for Nn in (30, 44, 60):
    D, R = map_to(0.0, Rc, Nn); D2=D@D; D3=D@D@D
    def resid(f):
        fp=D@f; fpp=D2@f; fppp=D3@f
        with np.errstate(all='ignore'):
            T2=40*(R*fpp-4*fp)/((R-2)*fp-2*f)
            T0h=(R*(R**4-54*R**2-54)*fppp-(R**3+18*R**2+12)*(R*fpp-fp)
                 +18*(R**2+2)*(fp+6*fpp))/(9*fpp-(R-3)*fp+2*f)
        return 768*pi**2*(2*f-R*fp) - (T2 - 36 - 12 - 5*R**2 + T0h)
    f0 = 0.005*(R**2) - 0.02*R - 0.005
    sol = root(resid, f0, method='lm', tol=1e-12, options={'maxiter':4000})
    f=sol.x; res=np.max(np.abs(resid(f)))
    idx0=np.argmin(np.abs(R)); g0=f[idx0]; g1=(D@f)[idx0]; g2=0.5*(D2@f)[idx0]
    glam = (1/(16*pi*g1))*(-g0/(2*g1)) if g1 else float('nan')
    conv = "CONVERGED" if res<1e-5 else "no"
    print(f"N={Nn:3d}: |res|={res:.2e} {conv} | g0={g0:+.4f} g1={g1:+.4f} g2={g2:+.4f} | g*lam*={glam:+.4f} | f/R^2@Rc={f[0]/R[0]**2:+.4f}")
    if res<1e-5: np.save('continuum/fstar_v2.npy', np.vstack([R,f]))
