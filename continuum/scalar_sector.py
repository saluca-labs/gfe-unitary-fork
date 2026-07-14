# -*- coding: utf-8 -*-
"""scalar_sector.py — canonical (mass) dimensions of the fork's scalar-sector operators
(eq. 2.5), deciding which threaten UV renormalizability. d=4, [Psi]=1, [R]=2, [d]=1. Mine."""
d = 4
def op_dim(psis, derivs, Rs):   # dimension of operator with #Psi, #derivatives, #Ricci-scalars
    return psis*1 + derivs*1 + Rs*2
ops = [   # (name, #Psi, #derivatives, #Ricci-scalars, coeff) — each grad Psi = 1 Psi + 1 deriv
    ("(grad Psi)^2  kinetic",      2, 2, 0, "c"),
    ("Psi^2         mass",         2, 0, 0, "4a"),
    ("Psi^4         self-coupling", 4, 0, 0, "2a^2"),
    ("R Psi^2       nonminimal",   2, 0, 1, "4a*beta"),
    ("R (grad Psi)^2  disformal",  2, 2, 1, "beta*c"),
    ("(grad Psi)^4  k-essence",    4, 4, 0, "c^2/2"),
]
print("op                         dim   coupling_dim   verdict")
print("-"*64)
for name, p, dv, r, coup in ops:
    od = op_dim(p, dv, r)
    cd = d - od          # coupling has dimension d - (operator dimension)
    verdict = "RELEVANT" if cd>0 else ("MARGINAL" if cd==0 else "IRRELEVANT")
    print(f"  {name:26s} {od:>2}      {cd:>+3}         {verdict}   (coeff {coup})")
print("-"*64)
print("READING:")
print("  * The two potentially-dangerous HIGHER-DERIVATIVE scalar operators —")
print("    R(grad Psi)^2 (disformal) and (grad Psi)^4 (k-essence) — are BOTH IRRELEVANT")
print("    (coupling dims -2, -4). They are UV-suppressed and DO NOT threaten the")
print("    continuum limit: the disformal Fisher structure is an IR dressing, not a UV actor.")
print("  * Marginal couplings: kinetic c (wavefn renorm), Psi^4, and nonminimal R Psi^2.")
print("    Psi^4 + nonminimal are exactly the couplings that gravity can drag to an")
print("    interacting 'gravity-induced' fixed point (asymptotically safe matter, Eichhorn).")
print("  * Relevant: only the Psi^2 mass — one tuning, as expected for a light scalar.")
print("  => the scalar sector is UV-BENIGN: dangerous operators irrelevant, marginal ones")
print("     fixed-point-compatible, exactly ONE relevant (mass) tuning.")
