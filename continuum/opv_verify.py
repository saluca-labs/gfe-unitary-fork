# -*- coding: utf-8 -*-
"""opv_verify.py — VERIFY the OPV (arXiv:1511.09393) exact quadratic fixed-point solutions
against the spectral-sum flow equation (3.18)+(3.19). If the residual vanishes identically in r
for the exact closed-form (alpha,beta,gamma,g0,g1,g2), both the transcription and the exact
solution are confirmed. Exponential parametrization, d=4, Litim cutoff."""
import sympy as sp
r = sp.symbols('r'); P = sp.pi**2

def resid(al, be, ga, g0, g1, g2):
    phi  = g0 + g1*r + g2*r**2
    phip = sp.diff(phi, r); phipp = sp.diff(phi, r, 2)   # phi''' = 0 for quadratic
    d1 = 5*(6+(6*al-1)*r)*(12+(12*al-1)*r)/(384*P)
    d2 = 5*(6+(6*al-1)*r)*(3+(3*al-2)*r)/(3456*P)
    d4 = ((2*be-1)*r+2)*((12*be+11)*r+12)/(256*P)
    d5 = (-72-18*r*(1+8*ga)+r**2*(19-18*ga-72*ga**2))/(192*P)
    # FP form of (3.18): dots=0, phi'''=0
    LHS = 4*phi - 2*r*phip
    RHS = ( d1/(6+(6*al+1)*r)
            - d2*(2*r*phipp - 2*phip)/phip
            + d4*phipp/((3+(3*be-1)*r)*phipp + phip)
            + d5/(4+(4*ga-1)*r) )
    return sp.simplify(sp.together(LHS - RHS))

s265 = sp.sqrt(265)
sol1 = dict(al=(11-s265)/54, be=(5*s265-73)/216, ga=(67-2*s265)/108,
            g0=(49+s265)/(1536*P), g1=-(4141+121*s265)/(82944*P), g2=(67795+3583*s265)/(4478976*P))
R = resid(**sol1)
num = sp.numer(R)
print("Solution 1 residual (LHS-RHS), simplified:")
print("   =", R)
print("   numerator simplified =", sp.simplify(sp.expand(num)))
print("   IS IDENTICALLY ZERO:", sp.simplify(R) == 0)
# also numeric sanity at a few r
f = sp.lambdify(r, R, 'mpmath')
import mpmath as mp; mp.mp.dps=30
print("   numeric residual at r=0.5,1,2,5:", [sp.nsimplify(0)+complex(f(v)) for v in (0.5,1,2,5)])
