# -*- coding: utf-8 -*-
"""fR_spectral_seeded.py — genuine attempt: solve global f*(R) by Chebyshev collocation of the
DM(3.1) type-Ib eq on [0,2.5] (BELOW the R=3 pole, so clean), SEEDED by the polynomial fixed
point (g0..g8, known good), and solved with a robust trust-region least-squares (not naive
Newton). If it converges to small residual with g*lambda* in band -> we have the global f*."""
import numpy as np
from scipy.optimize import least_squares
import sys; sys.path.insert(0,'continuum'); from cheb import map_to
pi=np.pi; Q=lambda a,b: a/b

# polynomial fixed point (from fR_fixedpoint.py, N=8) -> good seed
gpoly = np.array([-4.900838e-03,-2.277792e-02,3.833061e-03,4.947071e-04,
                  1.097310e-04,2.953745e-05,1.392327e-05,9.116159e-07,2.374868e-06])

Rmax=2.5
for Nn in (24,32,40):
    D,R = map_to(0.0,Rmax,Nn); D2=D@D; D3=D@D@D
    def rhs_resid(f):
        fp=D@f; fpp=D2@f; fppp=D3@f
        Dt=f+fp*(1-R/3); Ds=2*f+3*fp*(1-Q(2,3)*R)+9*fpp*(1-R/3)**2
        A=(5*R**2-(12+4*R-Q(61,90)*R**2))/(1-R/3)
        B=(10*R**2-R**2-(36+6*R-Q(67,60)*R**2))/(1-R/4)
        C=((2*fp-2*R*fpp)*(10-5*R-Q(271,36)*R**2+Q(7249,4536)*R**3)+fp*(60-20*R-Q(271,18)*R**2))/Dt
        Dd=(Q(5,2)*R**2)*((2*fp-2*R*fpp)*((1+R/3)+2*(1+R/6))+2*fp+4*fp)/Dt
        E=((2*fp-2*R*fpp)*fp*(6+3*R+Q(29,60)*R**2+Q(37,1512)*R**3)
           -2*R*fppp*(27-Q(91,20)*R**2-Q(29,30)*R**3-Q(181,3360)*R**4)
           +fpp*(216-Q(91,5)*R**2-Q(29,15)*R**3)+fp*(36+12*R+Q(29,30)*R**2))/Ds
        return 768*pi**2*(2*f-R*fp) - (A+B+C+Dd+E)
    f0 = sum(gpoly[k]*R**k for k in range(len(gpoly)))     # polynomial-FP seed on grid
    sol = least_squares(rhs_resid, f0, method='lm', xtol=1e-14, ftol=1e-14, max_nfev=20000)
    f=sol.x; res=np.max(np.abs(rhs_resid(f)))
    i0=np.argmin(np.abs(R)); g0=f[i0]; g1=(D@f)[i0]; g2=0.5*(D2@f)[i0]
    glam=(1/(16*pi*g1))*(-g0/(2*g1)) if g1 else float('nan')
    tag="CONVERGED" if res<1e-4 else "no"
    print(f"N={Nn:3d} Rmax={Rmax}: |res|={res:.2e} {tag} | g0={g0:+.5f} g1={g1:+.5f} g2={g2:+.5f} | g*lam*={glam:+.4f}")
    if res<1e-4:
        np.save('continuum/fstar_seeded.npy', np.vstack([R,f]))
        print(f"   -> saved global f*(R) on [0,{Rmax}] grid (N={Nn})")
