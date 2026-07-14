# -*- coding: utf-8 -*-
"""fR_improved.py — RG-IMPROVED f(R) flow: DM(2.9) with the implicit d_t f', d_t f'' terms.
Replacements vs fixed-point eq (3.1):  (2f'-2Rf'') -> (d_t f' + 2f' - 2Rf'')  in C,D,E ;
  -2R f''' -> (d_t f'' - 2R f''')  in E.  Flow is LINEAR in gd_m=d_t g_m at each order:
  sum_m A_nm(g) gd_m = c_n(g)  =>  beta = A^{-1} c.  FP: c=0 (=eq 3.1). Exponents:
  S = A^{-1} dc/dg ; theta = -eig(S)... (sign fixed by matching Gaussian/CPR). Validate vs
  CPR: 3 relevant directions, leading complex pair Re~2.4."""
import sympy as sp, numpy as np
R = sp.symbols('R'); pi = sp.pi; Q = sp.Rational
N = 6                                        # order (keep modest: implicit system is heavy)
g  = sp.symbols('g0:%d'%(N+1), real=True)
gd = sp.symbols('gd0:%d'%(N+1), real=True)   # d_t g_n
f   = sum(g[n]*R**n for n in range(N+1))
fp,fpp,fppp = sp.diff(f,R),sp.diff(f,R,2),sp.diff(f,R,3)
dtf = sum(gd[n]*R**n for n in range(N+1))
dtfp, dtfpp = sp.diff(dtf,R), sp.diff(dtf,R,2)

IMP1 = dtfp + 2*fp - 2*R*fpp                  # was (2f'-2Rf'')
IMP2 = dtfpp - 2*R*fppp                        # was -2Rf'''  (as the f''' additive piece)
Dt = f + fp*(1-R/3)
Ds = 2*f + 3*fp*(1-Q(2,3)*R) + 9*fpp*(1-R/3)**2
A_ = (5*R**2-(12+4*R-Q(61,90)*R**2))/(1-R/3)
B_ = (10*R**2-R**2-(36+6*R-Q(67,60)*R**2))/(1-R/4)
C_ = (IMP1*(10-5*R-Q(271,36)*R**2+Q(7249,4536)*R**3)+fp*(60-20*R-Q(271,18)*R**2))/Dt
D_ = (Q(5,2)*R**2)*(IMP1*((1+R/3)+2*(1+R/6))+2*fp+4*fp)/Dt
E_ = (IMP1*fp*(6+3*R+Q(29,60)*R**2+Q(37,1512)*R**3)
      + IMP2*(27-Q(91,20)*R**2-Q(29,30)*R**3-Q(181,3360)*R**4)
      + fpp*(216-Q(91,5)*R**2-Q(29,15)*R**3)+fp*(36+12*R+Q(29,30)*R**2))/Ds
RHS = A_+B_+C_+D_+E_
Flow = 384*pi**2*(dtf + 4*f - 2*R*fp) - RHS       # =0
print("expanding flow to order", N, "...")
fser = sp.series(Flow, R, 0, N+1).removeO()
Fn = [fser.coeff(R,n) for n in range(N+1)]        # N+1 eqns, linear in gd
# split A_nm (coeff of gd_m) and c_n (gd->0 part)
print("extracting A(g), c(g) ...")
Amat = sp.Matrix([[sp.diff(Fn[i], gd[j]) for j in range(N+1)] for i in range(N+1)])
cvec = sp.Matrix([Fn[i].subs({gd[j]:0 for j in range(N+1)}) for i in range(N+1)])
# fixed point: c(g)=0  (identical to eq 3.1)
cf = [sp.lambdify(g, cvec[i], 'mpmath') for i in range(N+1)]
seed=[-0.0049,-0.0228,0.00383,0.000495,0.00011,2.95e-5,1.39e-5][:N+1]
gstar = sp.nsolve(list(cvec), list(g), seed, prec=20, tol=1e-10, maxsteps=100)
sub = {g[n]:gstar[n] for n in range(N+1)}
print("fixed point g*:", [f"{float(gstar[n]):+.4e}" for n in range(N+1)])
G_=1/(16*np.pi*float(gstar[1])); Lam=-float(gstar[0])/(2*float(gstar[1]))
print(f"g*lambda* = {G_*Lam:+.4f}   (CPR band 0.05-0.14)")
# stability: beta = A^{-1} c... FlowEq: A gd + c = 0 => gd = -A^{-1} c = beta ; S=dbeta/dg=-A^{-1} dc/dg
An = np.array(Amat.subs(sub).evalf(20), dtype=complex)
dcdg = np.array(sp.Matrix([[sp.diff(cvec[i],g[j]) for j in range(N+1)] for i in range(N+1)]).subs(sub).evalf(20),dtype=complex)
S = -np.linalg.solve(An, dcdg)                    # dbeta/dg
theta = sorted(-np.linalg.eigvals(S), key=lambda z:-z.real)
print("\nRG-IMPROVED critical exponents theta:")
for t in theta: print(f"   {t.real:+.3f} {'+' if t.imag>=0 else '-'} {abs(t.imag):.3f} i")
nrel = sum(1 for t in theta if t.real>1e-4)
print(f"   # relevant = {nrel}   [CPR: 3]")
lead = [t for t in theta if t.real>0.1]
pair = [t for t in lead if abs(t.imag)>0.3]
print(f"   leading complex pair Re~{pair[0].real:.2f} Im~{abs(pair[0].imag):.2f} [CPR ~2.4+-2.5i]" if pair else "   no complex pair")
np.save('continuum/gstar.npy', np.array([float(gstar[n]) for n in range(N+1)]))
