# -*- coding: utf-8 -*-
"""opv_traj_fast.py — FULL OPV flow trajectory, purely numeric (no sympy). beta(g) via
polynomial Taylor division. Flow from exact NGFP (Solution 1) toward IR along relevant dirs;
track fork ratios r_n = g_n g1^(n-2)/g2^(n-1) vs fork {r3=4/3, r4=2, r5=16/5, r6=16/3}."""
import numpy as np, mpmath as mp
mp.mp.dps=30; pi2=float(mp.pi**2)
N=6
s=float(mp.sqrt(265))
al=(11-s)/54; be=(5*s-73)/216; ga=(67-2*s)/108
gstar=np.array([(49+s)/(1536*pi2), -(4141+121*s)/(82944*pi2), (67795+3583*s)/(4478976*pi2)]+[0.0]*(N-2))

def padd(*a):
    L=max(len(x) for x in a); o=np.zeros(L)
    for x in a: o[:len(x)]+=x
    return o
def der(c): return np.array([(k+1)*c[k+1] for k in range(len(c)-1)]) if len(c)>1 else np.array([0.0])
def shift(c,m): return np.concatenate([np.zeros(m),c]) if m>0 else c
def tdiv(num,den,N):
    nu=np.zeros(N+1); de=np.zeros(N+1)
    nu[:min(len(num),N+1)]=num[:N+1]; de[:min(len(den),N+1)]=den[:N+1]
    q=np.zeros(N+1); q[0]=nu[0]/de[0]
    for k in range(1,N+1): q[k]=(nu[k]-sum(de[j]*q[k-j] for j in range(1,k+1)))/de[0]
    return q
P=lambda *c: np.array(c,float)
d1=np.convolve(P(6,6*al-1),P(12,12*al-1))*5/(384*pi2)
d2=np.convolve(P(6,6*al-1),P(3,3*al-2))*5/(3456*pi2)
d3=np.convolve(np.convolve(P(2,2*be+3),P(3,3*be-1)),P(6,6*be-5))/(2304*pi2)
d4=np.convolve(P(2,2*be-1),P(12,12*be+11))/(256*pi2)
d5=P(-72,-18*(1+8*ga),19-18*ga-72*ga**2)/(192*pi2)
Da=P(6,6*al+1); Dg=P(4,4*ga-1)

def beta(g):
    phi=np.array(g,float); p1=der(phi); p2=der(p1); p3=der(p2)
    Ds=padd(np.convolve(P(3,3*be-1),p2), p1)          # (3+(3b-1)r)phi'' + phi'
    a=padd(shift(2*p2,1), -2*p1)                        # 2r phi'' - 2 phi'
    b=padd(np.convolve(d3,shift(-2*p3,1)), np.convolve(d4,p2))
    nonpd=padd(tdiv(d1,Da,N), -tdiv(np.convolve(d2,a),p1,N), tdiv(b,Ds,N), tdiv(d5,Dg,N))
    cpoly=padd(shift(2*p1,1), -4*phi, nonpd)
    c=np.array([cpoly[k] if k<len(cpoly) else 0.0 for k in range(N+1)])
    A=np.zeros((N+1,N+1))
    for m in range(N+1):
        col=np.zeros(N+1); col[m]=1.0
        if m>=1: col=padd(col, -tdiv(np.convolve(d2, shift(P(float(m)),m-1)), p1, N))[:N+1]
        if m>=2: col=padd(col, -tdiv(np.convolve(d3, shift(P(float(m*(m-1))),m-2)), Ds, N))[:N+1]
        A[:,m]=col[:N+1]
    return np.linalg.solve(A,c)

print("beta(g*) (should ~0):", np.round(beta(gstar),4))
J=np.zeros((N+1,N+1)); b0=beta(gstar); h=1e-7
for j in range(N+1):
    gp=gstar.copy(); gp[j]+=h; J[:,j]=(beta(gp)-b0)/h
mu=np.linalg.eigvals(J); th=sorted(-mu.real,reverse=True)
print("critical exponents theta=-eig(J):", [f"{t:+.3f}" for t in th], " [OPV: 4, 2.02 relevant]")
w,V=np.linalg.eig(J); idx=int(np.argmin(np.abs(-w-2.02))); v=np.real(V[:,idx]); v/=np.linalg.norm(v)
def rn(y,n): return y[n]*y[1]**(n-2)/y[2]**(n-1)
def rk4(y,dt):
    k1=beta(y);k2=beta(y+.5*dt*k1);k3=beta(y+.5*dt*k2);k4=beta(y+dt*k3); return y+dt/6*(k1+2*k2+2*k3+k4)
print("\nfork targets: r3=1.333 r4=2.000 r5=3.200 r6=5.333")
for eps in (3e-3,-3e-3):
    print(f"--- eps={eps:+.0e} along theta=2.02 relevant dir, toward IR ---")
    y=gstar+eps*v; t=0.0
    for step in range(1200):
        yn=rk4(y,-0.01)
        if not np.all(np.isfinite(yn)) or np.max(np.abs(yn))>1e4: print(f"   t={t:.2f}: blows up (Landau pole)"); break
        y=yn; t-=0.01
        if step%150==149: print(f"   t={t:+.2f}: r3={rn(y,3):+.4f} r4={rn(y,4):+.4f} r5={rn(y,5):+.4f} r6={rn(y,6):+.4f}")

# --- 2-parameter scan over the full relevant surface (theta=4 and theta=2.02 dirs) ---
print("\n=== 2-parameter relevant-surface scan: closest approach to fork {1.333,2.0,3.2,5.3} ===")
order=np.argsort(-(-w.real))   # by theta descending
i4=int(np.argmin(np.abs(-w-4.0))); i2=int(np.argmin(np.abs(-w-2.02)))
v4=np.real(V[:,i4]); v4/=np.linalg.norm(v4); v2=np.real(V[:,i2]); v2/=np.linalg.norm(v2)
best=(1e9,None)
for e4 in np.linspace(-6e-3,6e-3,9):
    for e2 in np.linspace(-6e-3,6e-3,9):
        if abs(e4)<1e-9 and abs(e2)<1e-9: continue
        y=gstar+e4*v4+e2*v2; t=0
        maxr3=0
        for step in range(1500):
            yn=rk4(y,-0.01)
            if not np.all(np.isfinite(yn)) or np.max(np.abs(yn))>1e4: break
            y=yn; t-=0.01
            r3=rn(y,3); 
            if abs(r3)>abs(maxr3): maxr3=r3
            dist=abs(rn(y,3)-4/3)+abs(rn(y,4)-2.0)
            if dist<best[0]: best=(dist,(e4,e2,t,rn(y,3),rn(y,4),rn(y,5),rn(y,6)))
print(f"closest approach to fork (min |r3-4/3|+|r4-2|): dist={best[0]:.3f}")
if best[1]:
    e4,e2,t,a3,a4,a5,a6=best[1]
    print(f"  at e4={e4:+.1e} e2={e2:+.1e} t={t:+.2f}:  r3={a3:+.3f} r4={a4:+.3f} r5={a5:+.3f} r6={a6:+.3f}")
    print(f"  fork:                                  r3=+1.333 r4=+2.000 r5=+3.200 r6=+5.333")
