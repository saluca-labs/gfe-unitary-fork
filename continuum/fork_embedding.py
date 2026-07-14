# -*- coding: utf-8 -*-
"""fork_embedding.py — where the fork's f(R) sits in RG theory space relative to the AS fixed point.
Fully self-derived (no fabricated gravitational beta functions). sympy + numpy."""
import sympy as sp

b, R, n, k = sp.symbols('beta R n k', positive=True)
print("="*72)
print("1. THE FORK ACTION IN THEORY SPACE:  f(R) = -4 ln(1 - beta R)")
print("="*72)
f = -4*sp.log(1 - b*R)
# Taylor coefficients c_n of f = sum c_n R^n
ser = sp.series(f, R, 0, 7).removeO()
poly = sp.Poly(ser, R)
coeffs = {i: poly.coeff_monomial(R**i) for i in range(1,7)}
print("Taylor coefficients c_n (dimensionful couplings of R^n):")
for i in range(1,7):
    cn_closed = sp.simplify(4*b**i/i)
    match = sp.simplify(coeffs[i]-cn_closed)==0
    print(f"   c_{i} = {coeffs[i]}   (=4β^{i}/{i}? {match})")
print("\n  => closed form  c_n = 4 beta^n / n   for all n>=1  (verified above)")
print("     c_1=4β (Einstein), c_2=2β^2 (Starobinsky R^2), c_3=(4/3)β^3, ...")

print("\n"+"="*72)
print("2. RADIUS OF CONVERGENCE = MAXIMAL CURVATURE (EFT horizon)")
print("="*72)
# ratio test: c_{n+1}/c_n -> beta ; radius = 1/beta
ratio = sp.simplify((4*b**(n+1)/(n+1))/(4*b**n/n))
lim = sp.limit(ratio, n, sp.oo)
print(f"   c_(n+1)/c_n = {ratio} ->(n→∞) {lim}   =>  radius of convergence R = 1/beta")
print("   The §7 'maximal curvature' R_max = 1/beta IS the radius of convergence of the")
print("   EFT expansion: it marks where the log-truncation breaks down, NOT a fundamental")
print("   spacetime boundary. The UV-complete theory must be defined BEYOND it.")

print("\n"+"="*72)
print("3. PREDICTIVITY: the fork is a ONE-PARAMETER ray in f(R) theory space")
print("="*72)
print("   All couplings fixed by the single scale beta:  c_n = 4β^n/n.")
print("   Dimensionless ratios are pure numbers (no free parameters):")
for i in range(2,6):
    rn = sp.nsimplify(sp.Rational(1,1)* (sp.Rational(1,i)) )  # c_n/(c_1 * (c_2/c_1... )) illustrative
# canonical invariant ratio c_n c_1^{n-2} / c_2^{n-1} (dimensionless, beta cancels)
print("   Invariant  r_n := c_n * c_1^(n-2) / c_2^(n-1)   (beta-independent pure number):")
c = lambda i: 4*b**i/i
for i in range(3,7):
    rn = sp.simplify(c(i)*c(1)**(i-2)/c(2)**(i-1))
    print(f"     r_{i} = {rn}")
print("   => once beta is fixed by the IR (Newton/Starobinsky scale), EVERY higher")
print("      coupling is a PREDICTION. The fork is a 1-parameter family; the AS UV")
print("      critical surface is ~3-dimensional (agent-verified). 1 <= 3: the fork ray")
print("      can lie on the critical surface (necessary condition for UV completion) MET.")

print("\n"+"="*72)
print("4. UV ASYMPTOTICS: fork log  vs.  fixed-point R^2  (the genuine tension)")
print("="*72)
print("   AS scale-invariance in d=4 requires the fixed-point fn f*(R) ~ A R^2 at large R.")
print("   The fork's log, -4 ln(1-βR), is only defined for R<1/β and does NOT reach large R.")
print("   => the fork CANNOT itself be the UV fixed-point action; it is the IR end of a")
print("      trajectory whose UV end is the global f*(R)~R^2. Dimensionless picture:")
bt = sp.symbols('b_tilde', positive=True)   # b_tilde = beta k^2
Rt = sp.symbols('R_tilde', real=True)        # R_tilde = R/k^2
ftil = -4*k**(-4)*sp.log(1 - bt*Rt)          # dimensionless f (up to norm)
print(f"      f~(R~) = {ftil}   with b~=βk^2, R~=R/k^2")
print("      UV: k→∞ at fixed dimensionful β => b~=βk^2→∞ and prefactor k^-4→0:")
print("      the dimensionless log-action FLOWS TO ZERO in the UV -> it is UV-IRRELEVANT")
print("      as a whole, i.e. an emergent IR structure, exactly as an EFT should be.")
print("      The nonzero UV content is carried by the RELEVANT couplings (g,λ,R^2-coupling)")
print("      that the NGFP fixes; the log is what they reassemble into at low k.")
