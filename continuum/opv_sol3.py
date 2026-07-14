# -*- coding: utf-8 -*-
"""Second validation: OPV Solution 3 (has a COMPLEX pair). Verify residual=0 then spectrum.
OPV Table 1: critical exponents 4, 2.4 +- 0.8 i, 0.13."""
import sympy as sp, numpy as np, mpmath as mp
mp.mp.dps=40
r=sp.symbols('r'); P=sp.Float(str(mp.pi**2),38)
u0,u1,u2,u3,w1,w2=sp.symbols('u0 u1 u2 u3 w1 w2')
def build(al,be,ga,g0,g1,g2):
    d1=5*(6+(6*al-1)*r)*(12+(12*al-1)*r)/(384*P); d2=5*(6+(6*al-1)*r)*(3+(3*al-2)*r)/(3456*P)
    d3=(2+(2*be+3)*r)*(3+(3*be-1)*r)*(6+(6*be-5)*r)/(2304*P); d4=((2*be-1)*r+2)*((12*be+11)*r+12)/(256*P)
    d5=(-72-18*r*(1+8*ga)+r**2*(19-18*ga-72*ga**2))/(192*P)
    RHS=(d1/(6+(6*al+1)*r)-d2*(2*r*u2-2*u1-w1)/u1+(d3*(w2-2*r*u3)+d4*u2)/((3+(3*be-1)*r)*u2+u1)+d5/(4+(4*ga-1)*r))
    phi=g0+g1*r+g2*r**2
    return RHS, {u0:phi,u1:sp.diff(phi,r),u2:sp.diff(phi,r,2),u3:sp.Integer(0),w1:sp.Integer(0),w2:sp.Integer(0)}
F=lambda x: sp.Float(str(x),36); s=mp.sqrt(1489); pi2=mp.pi**2
al,be,ga=F((2*s-41)/270),F(-(37*s+1559)/1080),F((143+4*s)/540)
g0,g1,g2=F((156-s)/(4224*pi2)),F((391*s-101773)/(1140480*pi2)),F((2479219-59293*s)/(153964800*pi2))
RHS,sub=build(al,be,ga,g0,g1,g2)
# residual check (FP form)
phi=g0+g1*r+g2*r**2
res=sp.simplify(sp.together((4*phi-2*r*sp.diff(phi,r)) - RHS.subs(sub)))
print("Solution 3 residual identically zero:", res==0, " (numeric @r=1,3:", 
      [complex(sp.lambdify(r,res,'mpmath')(v)) for v in (1,3)],")")
# spectrum
A={k:sp.diff(RHS,k).subs(sub) for k in (u0,u1,u2,u3,w1,w2)}
N=9; Bm=np.zeros((N+1,N+1)); Cm=np.zeros((N+1,N+1))
for m in range(N+1):
    v=r**m
    Bv=4*v-2*r*sp.diff(v,r)-(A[u0]*v+A[u1]*sp.diff(v,r)+A[u2]*sp.diff(v,r,2)+A[u3]*sp.diff(v,r,3))
    Cv=v-(A[w1]*sp.diff(v,r)+A[w2]*sp.diff(v,r,2))
    bs=sp.series(sp.together(Bv),r,0,N+1).removeO(); cs=sp.series(sp.together(Cv),r,0,N+1).removeO()
    for n in range(N+1): Bm[n,m]=float(bs.coeff(r,n)); Cm[n,m]=float(cs.coeff(r,n))
th=sorted(np.linalg.eigvals(np.linalg.solve(Cm,Bm)),key=lambda z:-z.real)
print("Solution 3 spectrum [OPV Table1: 4, 2.4+-0.8i, 0.13]:")
for t in th[:5]: print(f"    {t.real:+.3f}{'+' if t.imag>=0 else '-'}{abs(t.imag):.3f}i"+("  <-- relevant" if t.real>1e-3 else ""))
print(f"    relevant count = {sum(1 for t in th if t.real>1e-3)}")
