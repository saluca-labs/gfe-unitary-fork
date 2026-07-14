# -*- coding: utf-8 -*-
"""opv_trajectory.py — integrate the FULL OPV flow from the exact NGFP (Solution 1) toward the
IR along its relevant directions, and test whether the fork's coupling ratios r_n = g_n g1^(n-2)
/g2^(n-1) = {4/3, 2, 16/5, 16/3} emerge. r_n are dimensionless & k-independent along the flow.
Flow ODEs: A(g) gdot = c(g)  (implicit from the phi_dot',phi_dot'' improvement terms)."""
import sympy as sp, numpy as np, mpmath as mp
from scipy.integrate import solve_ivp
mp.mp.dps=30
r=sp.symbols('r'); P=sp.Float(str(mp.pi**2),30)
N=4
g=sp.symbols('g0:%d'%(N+1)); gd=sp.symbols('gd0:%d'%(N+1))
F=lambda x: sp.Float(str(x),28); s=mp.sqrt(265); pi2=mp.pi**2
al,be,ga=F((11-s)/54),F((5*s-73)/216),F((67-2*s)/108)
gstar=[F((49+s)/(1536*pi2)),F(-(4141+121*s)/(82944*pi2)),F((67795+3583*s)/(4478976*pi2))]+[F(0)]*(N-2)
phi=sum(g[n]*r**n for n in range(N+1)); phid=sum(gd[n]*r**n for n in range(N+1))
pp=sp.diff(phi,r); ppp=sp.diff(phi,r,2); pppp=sp.diff(phi,r,3)
pdp=sp.diff(phid,r); pdpp=sp.diff(phid,r,2)
d1=5*(6+(6*al-1)*r)*(12+(12*al-1)*r)/(384*P); d2=5*(6+(6*al-1)*r)*(3+(3*al-2)*r)/(3456*P)
d3=(2+(2*be+3)*r)*(3+(3*be-1)*r)*(6+(6*be-5)*r)/(2304*P); d4=((2*be-1)*r+2)*((12*be+11)*r+12)/(256*P)
d5=(-72-18*r*(1+8*ga)+r**2*(19-18*ga-72*ga**2))/(192*P)
RHS=(d1/(6+(6*al+1)*r)-d2*(2*r*ppp-2*pp-pdp)/pp+(d3*(pdpp-2*r*pppp)+d4*ppp)/((3+(3*be-1)*r)*ppp+pp)+d5/(4+(4*ga-1)*r))
Fexpr=phid-2*r*pp+4*phi-RHS
print("expanding flow residual to order",N,"...")
ser=sp.series(sp.together(Fexpr),r,0,N+1).removeO(); Pn=sp.Poly(sp.numer(sp.together(ser)),r)
# Actually take series coeffs directly (denominator nonzero near r=0 at FP):
Fn=[ser.coeff(r,n) for n in range(N+1)]; print("Fn built; lambdifying...")
A=sp.Matrix([[sp.diff(Fn[i],gd[j]) for j in range(N+1)] for i in range(N+1)])
c=sp.Matrix([-Fn[i].subs({gd[j]:0 for j in range(N+1)}) for i in range(N+1)])
Af=sp.lambdify([g],A,'numpy'); cf=sp.lambdify([g],c,'numpy')
def beta(gv):
    return np.linalg.solve(np.array(Af(gv),float).reshape(N+1,N+1), np.array(cf(gv),float).ravel())
g0=np.array([float(x) for x in gstar])
# Jacobian of beta at FP -> relevant eigenvectors (mu=-theta<0 grow toward IR)
J=np.zeros((N+1,N+1)); h=1e-6
b0=beta(g0)
for j in range(N+1):
    gp=g0.copy(); gp[j]+=h; J[:,j]=(beta(gp)-b0)/h
mu,V=np.linalg.eig(J)
print("beta-Jacobian eigenvalues mu (theta=-mu):",np.round(sorted(mu.real),3))
# pick the relevant non-constant direction: mu approx -2.02 (theta=+2.02)
idx=np.argmin(np.abs(mu-(-2.02)))
v=np.real(V[:,idx]); v=v/np.linalg.norm(v)
print(f"perturbing along theta={-mu[idx].real:.2f} relevant direction")
def rk4(y,dt):
    k1=beta(y); k2=beta(y+0.5*dt*k1); k3=beta(y+0.5*dt*k2); k4=beta(y+dt*k3)
    return y+(dt/6)*(k1+2*k2+2*k3+k4)
def rn(y,n): return y[n]*y[1]**(n-2)/y[2]**(n-1)
print("fork targets: r3=1.333  r4=2.000")
for eps in (2e-3,-2e-3):
    print(f"
--- trajectory from NGFP along theta=2.08 relevant dir, eps={eps:+.0e}, toward IR ---")
    y=g0+eps*v; t=0.0; dt=-0.01
    for step in range(600):
        yn=rk4(y,dt)
        if not np.all(np.isfinite(yn)) or np.max(np.abs(yn))>1e3:
            print(f"   t={t:.2f}: flow blows up (Landau pole) -> stop"); break
        y=yn; t+=dt
        if step%50==49:
            print(f"   t={t:+.2f} (k=e^t): r3={rn(y,3):+.4f}  r4={rn(y,4):+.4f}  |g3={y[3]:+.2e} g4={y[4]:+.2e}")
