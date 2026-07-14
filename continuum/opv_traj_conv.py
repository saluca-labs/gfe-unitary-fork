# -*- coding: utf-8 -*-
"""opv_traj_conv.py — convergence study: fork trajectory test at N=6,8,10,12. For each N,
scan the 2-dim relevant surface, report closest-approach r3..r6 to the fork and the Landau-pole
location. Converge toward fork -> artifact; converge away -> physical mismatch; no converge ->
unresolved. Numeric beta(g) via Taylor division; NGFP = OPV Solution 1."""
import numpy as np, mpmath as mp
mp.mp.dps=30; pi2=float(mp.pi**2)
s=float(mp.sqrt(265)); al=(11-s)/54; be=(5*s-73)/216; ga=(67-2*s)/108
g012=[(49+s)/(1536*pi2), -(4141+121*s)/(82944*pi2), (67795+3583*s)/(4478976*pi2)]
P=lambda *c: np.array(c,float)
d1=np.convolve(P(6,6*al-1),P(12,12*al-1))*5/(384*pi2); d2=np.convolve(P(6,6*al-1),P(3,3*al-2))*5/(3456*pi2)
d3=np.convolve(np.convolve(P(2,2*be+3),P(3,3*be-1)),P(6,6*be-5))/(2304*pi2); d4=np.convolve(P(2,2*be-1),P(12,12*be+11))/(256*pi2)
d5=P(-72,-18*(1+8*ga),19-18*ga-72*ga**2)/(192*pi2); Da=P(6,6*al+1); Dg=P(4,4*ga-1)
def padd(*a):
    L=max(len(x) for x in a); o=np.zeros(L)
    for x in a: o[:len(x)]+=x
    return o
def der(c): return np.array([(k+1)*c[k+1] for k in range(len(c)-1)]) if len(c)>1 else np.array([0.0])
def shift(c,m): return np.concatenate([np.zeros(m),c]) if m>0 else c
def tdiv(num,den,N):
    nu=np.zeros(N+1); de=np.zeros(N+1); nu[:min(len(num),N+1)]=num[:N+1]; de[:min(len(den),N+1)]=den[:N+1]
    q=np.zeros(N+1); q[0]=nu[0]/de[0]
    for k in range(1,N+1): q[k]=(nu[k]-de[1:k+1]@q[k-1::-1][:k])/de[0]
    return q
def make_beta(N):
    def beta(g):
        phi=np.array(g,float); p1=der(phi); p2=der(p1); p3=der(p2)
        Ds=padd(np.convolve(P(3,3*be-1),p2),p1); a=padd(shift(2*p2,1),-2*p1)
        b=padd(np.convolve(d3,shift(-2*p3,1)),np.convolve(d4,p2))
        nonpd=padd(tdiv(d1,Da,N),-tdiv(np.convolve(d2,a),p1,N),tdiv(b,Ds,N),tdiv(d5,Dg,N))
        cpoly=padd(shift(2*p1,1),-4*phi,nonpd); c=np.array([cpoly[k] if k<len(cpoly) else 0.0 for k in range(N+1)])
        A=np.zeros((N+1,N+1))
        for m in range(N+1):
            col=np.zeros(N+1); col[m]=1.0
            if m>=1: col=padd(col,-tdiv(np.convolve(d2,shift(P(float(m)),m-1)),p1,N))[:N+1]
            if m>=2: col=padd(col,-tdiv(np.convolve(d3,shift(P(float(m*(m-1))),m-2)),Ds,N))[:N+1]
            A[:,m]=col[:N+1]
        return np.linalg.solve(A,c)
    return beta
def rn(y,n): return y[n]*y[1]**(n-2)/y[2]**(n-1)
print("N | theta_2 | min-dist | r3 r4 r5 r6 (fork:1.333 2.000 3.200 5.333) | t_pole")
for N in (6,8,10,12):
    gstar=np.array(g012+[0.0]*(N-2)); beta=make_beta(N)
    J=np.zeros((N+1,N+1)); b0=beta(gstar); h=1e-7
    for j in range(N+1):
        gp=gstar.copy(); gp[j]+=h; J[:,j]=(beta(gp)-b0)/h
    w,V=np.linalg.eig(J); th2=sorted(-w.real,reverse=True)[1]
    i4=int(np.argmin(np.abs(-w-4))); i2=int(np.argmin(np.abs(-w-2.02)))
    v4=np.real(V[:,i4])/np.linalg.norm(np.real(V[:,i4])); v2=np.real(V[:,i2])/np.linalg.norm(np.real(V[:,i2]))
    def rk4(y,dt): 
        k1=beta(y);k2=beta(y+.5*dt*k1);k3=beta(y+.5*dt*k2);k4=beta(y+dt*k3); return y+dt/6*(k1+2*k2+2*k3+k4)
    best=(1e9,None); 
    for e4 in np.linspace(-7e-3,7e-3,9):
        for e2 in np.linspace(-7e-3,7e-3,9):
            if abs(e4)+abs(e2)<1e-9: continue
            y=gstar+e4*v4+e2*v2; t=0; tp=0
            for st in range(1500):
                yn=rk4(y,-0.01)
                if not np.all(np.isfinite(yn)) or np.max(np.abs(yn))>1e4: tp=t; break
                y=yn; t-=0.01; tp=t
                dist=abs(rn(y,3)-4/3)+abs(rn(y,4)-2.0)
                if dist<best[0]: best=(dist,(rn(y,3),rn(y,4),rn(y,5),rn(y,6),t))
    r3,r4,r5,r6,tt=best[1]
    print(f"{N:2d}| {th2:+.3f} | {best[0]:.3f} | {r3:+.2f} {r4:+.2f} {r5:+.2f} {r6:+.2f} | t={tt:+.2f}")
