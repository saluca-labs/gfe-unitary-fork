# -*- coding: utf-8 -*-
"""Physical-branch selection: among shooting solutions that pass R_c, the true fixed point is
the one with correct far-field f ~ A R^2 asymptotics. Integrate to R=25, report f/R^2 stability
and g*lambda* for survivors. A converged f/R^2 + physical-band g*lambda* = the NGFP."""
import numpy as np
from scipy.integrate import solve_ivp
pi=np.pi; C=768*pi**2
def rhs_known(R,f,fp,fpp):
    return C*(2*f-R*fp) - 40*(R*fpp-4*fp)/((R-2)*fp-2*f) + 36 + 12 + 5*R**2
def fppp(R,f,fp,fpp):
    Ds=9*fpp-(R-3)*fp+2*f
    return (rhs_known(R,f,fp,fpp)*Ds+(R**3+18*R**2+12)*(R*fpp-fp)-18*(R**2+2)*(fp+6*fpp))/(R*(R**4-54*R**2-54))
def fpp0_reg(f0,fp0):
    rk0=C*2*f0 - 40*(-4*fp0)/(-2*fp0-2*f0)+48
    A=9*rk0-216; B=rk0*(3*fp0+2*f0)-48*fp0
    return -B/A if abs(A)>1e-30 else np.nan
def integrate(f0,fp0,Rend=25.0):
    fpp0=fpp0_reg(f0,fp0)
    if not np.isfinite(fpp0): return None
    def ode(R,y): return [y[1],y[2],fppp(R,y[0],y[1],y[2])]
    def blow(R,y): return 1e4-abs(y[0])-abs(y[2])
    blow.terminal=True; blow.direction=-1
    try:
        s=solve_ivp(ode,[0.05,Rend],[f0,fp0,fpp0],method='RK45',events=blow,
                    rtol=1e-7,atol=1e-10,dense_output=True)
        return s
    except Exception: return None

Rc=np.sqrt(27+3*np.sqrt(87))
f0s=np.linspace(-0.03,0.02,22); f1s=np.linspace(-0.035,0.015,22)
cands=[]
for f0 in f0s:
    for f1 in f1s:
        s=integrate(f0,f1)
        if s is None or s.t[-1]<Rc+2: continue          # must pass R_c and reach far
        R2=s.t[-1]
        if R2<20: continue
        a1=s.sol(18)[0]/18**2; a2=s.sol(22)[0]/22**2      # f/R^2 at two large R
        gl=(1/(16*pi*f1))*(-f0/(2*f1)) if f1 else np.nan
        conv=abs(a1-a2)/(abs(a2)+1e-12)
        cands.append((conv,f0,f1,a1,a2,gl,R2))
cands.sort()
print(f"pass-through solutions reaching R>=20: {len(cands)}")
print("best by f/R^2 convergence (conv=|a(18)-a(22)|/|a(22)|):")
for conv,f0,f1,a1,a2,gl,R2 in cands[:12]:
    flag=" <-- physical-band g*lam*" if 0.03<abs(gl)<0.2 else ""
    print(f"  f0={f0:+.4f} f1={f1:+.4f} | f/R^2: {a1:+.4f}->{a2:+.4f} conv={conv:.3f} | g*lam*={gl:+.3f} | reach={R2:.1f}{flag}")
