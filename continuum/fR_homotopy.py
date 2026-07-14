# -*- coding: utf-8 -*-
"""fR_homotopy.py — homotopy-continuation solver for the global f(R) fixed point.
E_tau = 768pi^2(2f-Rf') - (A+B) - tau*(C+D+E).  tau=0: linear (A,B are R-only) => solvable.
Continue tau 0->1 (each solution seeds the next) to the full DM(3.1) type-Ib eq. Domain
[0,2.8] (below the R=3 pole). If tau=1 residual -> 0, we have the global f* the direct solves
could not reach. Then (part 2) the critical spectrum on the grid."""
import numpy as np
from scipy.optimize import least_squares
import sys; sys.path.insert(0,'continuum'); from cheb import map_to
pi=np.pi; Q=lambda a,b:a/b
Rmax=2.8; Nn=28
D,R=map_to(0.0,Rmax,Nn); D2=D@D; D3=D@D@D

def parts(f):
    fp=D@f; fpp=D2@f; fppp=D3@f
    Dt=f+fp*(1-R/3); Ds=2*f+3*fp*(1-Q(2,3)*R)+9*fpp*(1-R/3)**2
    A=(5*R**2-(12+4*R-Q(61,90)*R**2))/(1-R/3)
    B=(10*R**2-R**2-(36+6*R-Q(67,60)*R**2))/(1-R/4)
    C=((2*fp-2*R*fpp)*(10-5*R-Q(271,36)*R**2+Q(7249,4536)*R**3)+fp*(60-20*R-Q(271,18)*R**2))/Dt
    Dd=(Q(5,2)*R**2)*((2*fp-2*R*fpp)*((1+R/3)+2*(1+R/6))+2*fp+4*fp)/Dt
    E=((2*fp-2*R*fpp)*fp*(6+3*R+Q(29,60)*R**2+Q(37,1512)*R**3)
       -2*R*fppp*(27-Q(91,20)*R**2-Q(29,30)*R**3-Q(181,3360)*R**4)
       +fpp*(216-Q(91,5)*R**2-Q(29,15)*R**3)+fp*(36+12*R+Q(29,30)*R**2))/Ds
    return 768*pi**2*(2*f-R*fp), A, B, C, Dd, E

def resid_tau(f, tau):
    LHS,A,B,C,Dd,E = parts(f)
    return LHS - (A+B) - tau*(C+Dd+E)

# seed: polynomial FP evaluated on grid
gpoly=np.array([-4.900838e-03,-2.277792e-02,3.833061e-03,4.947071e-04,1.097310e-04,
                2.953745e-05,1.392327e-05,9.116159e-07,2.374868e-06])
f=sum(gpoly[k]*R**k for k in range(len(gpoly)))
taus=np.concatenate([np.linspace(0,0.85,12),np.linspace(0.9,1.0,6)])
print("homotopy continuation tau: 0 -> 1")
for tau in taus:
    sol=least_squares(resid_tau,f,args=(tau,),method='trf',xtol=1e-13,ftol=1e-13,max_nfev=3000)
    f=sol.x; res=np.max(np.abs(resid_tau(f,tau)))
    if abs(tau-round(tau,2))<1e-9 and (abs(tau*20-round(tau*20))<1e-6):
        pass
    if True:
        i0=np.argmin(np.abs(R)); g0=f[i0]; g1=(D@f)[i0]
        glam=(1/(16*pi*g1))*(-g0/(2*g1)) if g1 else float('nan')
        print(f"  tau={tau:5.3f}: |res|={res:.2e}  g0={g0:+.5f} g1={g1:+.5f} g*lam*={glam:+.4f}")
res1=np.max(np.abs(resid_tau(f,1.0)))
print(f"\nFINAL tau=1: |res|={res1:.2e}  =>", "CONVERGED (global f* obtained)" if res1<1e-4 else "NOT converged")
if res1<1e-4: np.save('continuum/fstar_homotopy.npy', np.vstack([R,f]))
