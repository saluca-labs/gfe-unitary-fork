# -*- coding: utf-8 -*-
"""opv_spectrum.py — critical spectrum around the OPV exact quadratic fixed point. Linearize the
spectral-sum flow (3.18); eigenperturbation phi=phi*+eps e^{-theta t} v(r) gives a generalized
eigenvalue problem B[v]=theta C[v]. Polynomial basis v=sum a_n r^n (OPV use 9th order). VALIDATE
against Table 1: Solution 1 -> (4, 2.02) relevant; sanity: v=const must give theta=4."""
import sympy as sp, numpy as np, mpmath as mp
mp.mp.dps=40
r=sp.symbols('r'); P=sp.Float(str(mp.pi**2),38)
u0,u1,u2,u3,w1,w2=sp.symbols('u0 u1 u2 u3 w1 w2')  # phi,phi',phi'',phi''',phidot',phidot''

def spectrum(al,be,ga,g0,g1,g2,N,label):
    d1=5*(6+(6*al-1)*r)*(12+(12*al-1)*r)/(384*P)
    d2=5*(6+(6*al-1)*r)*(3+(3*al-2)*r)/(3456*P)
    d3=(2+(2*be+3)*r)*(3+(3*be-1)*r)*(6+(6*be-5)*r)/(2304*P)
    d4=((2*be-1)*r+2)*((12*be+11)*r+12)/(256*P)
    d5=(-72-18*r*(1+8*ga)+r**2*(19-18*ga-72*ga**2))/(192*P)
    RHS=( d1/(6+(6*al+1)*r) - d2*(2*r*u2-2*u1-w1)/u1
          + (d3*(w2-2*r*u3)+d4*u2)/((3+(3*be-1)*r)*u2+u1) + d5/(4+(4*ga-1)*r) )
    phi=g0+g1*r+g2*r**2
    sub={u0:phi,u1:sp.diff(phi,r),u2:sp.diff(phi,r,2),u3:0,w1:0,w2:0}
    A={k:sp.diff(RHS,k).subs(sub) for k in (u0,u1,u2,u3,w1,w2)}
    # operators B[v]=4v-2r v'-(Au0 v+Au1 v'+Au2 v''+Au3 v'''),  C[v]=v-(Aw1 v'+Aw2 v'')
    def Bop(v): 
        return 4*v-2*r*sp.diff(v,r)-(A[u0]*v+A[u1]*sp.diff(v,r)+A[u2]*sp.diff(v,r,2)+A[u3]*sp.diff(v,r,3))
    def Cop(v):
        return v-(A[w1]*sp.diff(v,r)+A[w2]*sp.diff(v,r,2))
    Bm=np.zeros((N+1,N+1)); Cm=np.zeros((N+1,N+1))
    for m in range(N+1):
        v=r**m
        bs=sp.series(sp.together(Bop(v)),r,0,N+1).removeO()
        cs=sp.series(sp.together(Cop(v)),r,0,N+1).removeO()
        for n in range(N+1):
            Bm[n,m]=float(bs.coeff(r,n)); Cm[n,m]=float(cs.coeff(r,n))
    ev=np.linalg.eigvals(np.linalg.solve(Cm,Bm))
    th=sorted(ev,key=lambda z:-z.real)
    print(f"{label} (N={N}): critical exponents theta (Re>0=relevant):")
    for t in th:
        print(f"    {t.real:+.3f}{'+' if t.imag>=0 else '-'}{abs(t.imag):.3f}i" + ("  <-- relevant" if t.real>1e-3 else ""))
    rel=[t for t in th if t.real>1e-3]
    print(f"    => {len(rel)} relevant; leading = {th[0].real:+.3f}")

def F(x): return sp.Float(str(x),36)
s=mp.sqrt(265); pi2=mp.pi**2
spectrum(F((11-s)/54),F((5*s-73)/216),F((67-2*s)/108),
         F((49+s)/(1536*pi2)),F(-(4141+121*s)/(82944*pi2)),F((67795+3583*s)/(4478976*pi2)),
         9,"Solution 1 [OPV Table1: 4, 2.02]")
