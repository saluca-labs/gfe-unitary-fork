# -*- coding: utf-8 -*-
"""fR_shooting.py — spike-plot SHOOTING for the global f(R) fixed point (Dietz-Morris/Falls
method), Benedetti-Caravelli type-II eq. Solve f'''=F(R,f,f',f''); regularity at R=0 fixes
f''(0) from (f(0),f'(0)); integrate outward as an IVP; the global solution appears as a SPIKE
where the integration extends furthest before hitting a movable singularity. Scan (f0,f1)."""
import numpy as np
from scipy.integrate import solve_ivp
pi=np.pi; C=768*pi**2

def rhs_known(R,f,fp,fpp):
    T2=40*(R*fpp-4*fp)/((R-2)*fp-2*f)
    return C*(2*f-R*fp) - T2 - (-36) - (-12-5*R**2)

def fppp(R,f,fp,fpp):
    Ds=9*fpp-(R-3)*fp+2*f
    num=rhs_known(R,f,fp,fpp)*Ds + (R**3+18*R**2+12)*(R*fpp-fp) - 18*(R**2+2)*(fp+6*fpp)
    return num/(R*(R**4-54*R**2-54))

def fpp0_reg(f0,fp0):
    # regularity at R=0: numerator(0)=0, linear in fpp0
    rk0=C*(2*f0) - 40*(-4*fp0)/((-2)*fp0-2*f0) + 36 + 12   # rhs_known(0)
    # num(0)= rk0*(9 fpp0+3 fp0+2 f0) + 12*(-fp0) - 36*(fp0+6 fpp0) = 0
    A=9*rk0-216; B=rk0*(3*fp0+2*f0) -12*fp0 -36*fp0
    return -B/A if abs(A)>1e-30 else np.nan

def reach(f0,fp0,Rend=10.0):
    fpp0=fpp0_reg(f0,fp0)
    if not np.isfinite(fpp0): return 0.0
    eps=0.05
    def ode(R,y): return [y[1],y[2],fppp(R,y[0],y[1],y[2])]
    def blow(R,y): return 1e3-abs(y[0])-abs(y[2])
    blow.terminal=True; blow.direction=-1
    try:
        s=solve_ivp(ode,[eps,Rend],[f0,fp0,fpp0],method='RK45',
                    events=blow,rtol=1e-6,atol=1e-9)
        return s.t[-1]
    except Exception:
        return eps

# coarse scan near the polynomial-FP region (g0~-0.005, g1~-0.023), broadened
f0s=np.linspace(-0.03,0.02,20)
f1s=np.linspace(-0.06,0.02,20)
best=(0,None,None); grid=np.zeros((len(f0s),len(f1s)))
for i,f0 in enumerate(f0s):
    for j,f1 in enumerate(f1s):
        r=reach(f0,f1); grid[i,j]=r
        if r>best[0]: best=(r,f0,f1)
print(f"max reach R={best[0]:.3f} at f(0)={best[1]:+.4f}, f'(0)={best[2]:+.4f}  (R_c=7.415; Rend=7)")
# how many params reach far (spike sharpness)
print(f"params reaching R>5: {int((grid>5).sum())} / {grid.size};  R>6.5: {int((grid>6.5).sum())}")
# show the far-reaching band
ii,jj=np.where(grid>max(5.0,0.8*best[0]))
if len(ii):
    print("far-reaching (f0,f1,reach) samples:")
    for k in range(0,len(ii),max(1,len(ii)//8)):
        print(f"   f0={f0s[ii[k]]:+.4f} f1={f1s[jj[k]]:+.4f} -> R={grid[ii[k],jj[k]]:.3f}")
Rc=np.sqrt(27+3*np.sqrt(87))
print(f"R_c (fixed singularity) = {Rc:.4f}")
print(f"params blowing up BEFORE R_c ({Rc:.2f}): {int((grid<Rc-0.05).sum())}/{grid.size}")
print(f"params PASSING R_c (reach > {Rc+0.3:.2f}): {int((grid>Rc+0.3).sum())}/{grid.size}")
hist,edges=np.histogram(grid.ravel(),bins=[0,2,4,6,7,7.3,7.42,8,9,10.01])
print("R_blow histogram:", dict(zip([f'{edges[k]:.1f}-{edges[k+1]:.1f}' for k in range(len(hist))],hist.tolist())))
ii,jj=np.where(grid>Rc+0.3)
if len(ii):
    print("PASS-THROUGH params (candidate fixed points):")
    for k in range(len(ii)):
        print(f"   f0={f0s[ii[k]]:+.4f} f1={f1s[jj[k]]:+.4f} -> R_blow={grid[ii[k],jj[k]]:.3f}")
np.save('continuum/shoot_grid.npy', grid)
