# Computation scripts — GfE unitary fork

These are the computations behind the paper *A Ghost-Free, Unitary, Starobinsky-Class
Realization of Gravity from Entropy*. They are archived in full, **including the attempts that
failed**, because the honest documentation of what did and did not work is part of the record.

Run with Python 3.12 (numpy, sympy, scipy, mpmath). On Windows set `PYTHONIOENCODING=utf-8`.

## Ghost-freedom / cosmology (in the main text; verified independently)
Symbolic tensor algebra and cosmology were done separately (sympy / Cadabra2 in WSL; hi_class).
Those results (three-way ghost-freedom proof, det G = 3βc/8 > 0, Starobinsky n_s, thawing dark
energy) are stated in the paper.

## The interacting continuum limit (§8) — what these scripts establish

| script | what it does | outcome |
|---|---|---|
| `fork_embedding.py` | fork f(R)=−4ln(1−βR) in RG theory space: Taylor c_n=4βⁿ/n, radius of convergence = 1/β = R_max, invariant ratios r_n, log-vs-R² UV asymptotics | **clean, self-derived** |
| `ngfp_mechanism.py` | asymptotic-safety mechanism, β_g = 2g − cg² ⇒ g*=2/c, θ=2 | schematic, exact-in-toy |
| `scalar_sector.py` | canonical dimensions of the Ψ-sector operators; disformal & k-essence are irrelevant ⇒ UV-benign | clean (a field-count bug was caught vs the reading text and fixed) |
| `cheb.py` | Chebyshev pseudo-spectral differentiation; **validated** on known spectra (derivs 1e-14; −u″=λu eigenvalues 1,4,9,…) | **instrument validated** |
| `fR_fixedpoint.py` | solves the exact Dietz–Morris f(R) fixed-point eq (3.1) in polynomial truncation | **reproduces the NGFP: g*λ* ≈ 0.09–0.1, in the CPR band** |
| `fR_exponents_fork.py` | exponents in a simplified (no RG-improvement) scheme | **fails** — 7 relevant vs 3; shows improvement terms are essential |
| `fR_improved.py` | RG-improved flow (implicit ∂_t f′, ∂_t f″), β = −A⁻¹c | reproduces the FP; exponents still spurious (small-propagator-denominator instability, diagnosed in the script) |
| `fR_spectral_fp.py`, `fR_spectral_v2.py` | global f(R) fixed point via Chebyshev collocation (Benedetti–Caravelli type-II, pole-free) | **does not converge** under naive Newton — the singular BVP needs parameter-continuation machinery |
| `lpa_flow.py`, `lpa_scan.py`, `lpa_robust.py` | scalar LPA Wilson–Fisher as an independent machinery check | eigenvalue code validated on the Gaussian; the WF polynomial root did not converge (known-finicky) |

## Honest bottom line
The f(R) non-Gaussian fixed point is **reproduced and validated**, and the fork's f(R) class
sits at it. The fork's *exact* critical spectrum and trajectory were **not** obtained: four
principled numerical attempts each failed for a well-understood reason, and a verified search
found **no turn-key public code** for the f(R) fixed-point spectrum (the field runs private
Mathematica). That computation is well-posed but externally gated — the honest open edge stated
in §8. No exponent, trajectory, or fixed point was ever fabricated.

## Update (2026-07-14): homotopy-continuation attempt

`fR_homotopy.py` deforms the equation E_τ = 768π²(2f−Rf′) − (A+B) − τ(C+D+E) from a solvable
linear limit (τ=0, threshold terms off) to the full equation (τ=1), tracking the solution.
Result: **τ=0 solves cleanly (residual ~1e−4), but for any τ>0 the residual grows linearly while
the solution freezes** — least-squares finds no descent direction once the f-dependent (singular)
terms are on. This isolates the blocker precisely: it is **not** the seed/basin/global convergence
(homotopy addresses those), but the residual-at-all-nodes collocation *formulation*, which is
ill-posed for this singular boundary-value problem. The correct treatment imposes regularity
conditions as explicit boundary conditions at the singular points plus far-field matching
(shooting / quantization), i.e. the dedicated machinery. Six principled attempts, one wall.

## Update (2026-07-14, cont.): shooting / spike-plot (`fR_shooting.py`, `fR_shoot_select.py`)

The shooting method — solve f‴=F(R,f,f′,f″), fix f″(0) by regularity at R=0, integrate outward,
select by passage through the fixed singularity R_c=√(27+3√87)=7.415 — is the *correct*
formulation for this singular BVP, and it **produces genuine spike-plot structure**: 313/400
initial conditions blow up before R_c, a clear cluster stalls at R_c, a subset passes through.
This is real behavior collocation never produced.

However it does **not** yield the physical fixed point: the pass-through solutions are all spurious
(unphysical g\*λ\*, no f∼R² asymptotics), and the physical NGFP data blows up at a movable
singularity near the origin. This is consistent with the known result that **global single-metric
f(R) fixed-function solutions are pathological / scheme-dependent** (Dietz–Morris, arXiv:1211.0955),
which is precisely why the literature uses the exponential parametrization (Ohta–Percacci–Vacca,
arXiv:1511.09393) to obtain global solutions. Seven principled methods; the critical spectrum
remains gated behind either that reparametrized equation or dedicated code. Nothing fabricated.

## Update (2026-07-14, cont.): critical spectrum SOLVED via exponential parametrization

`opv_verify.py`, `opv_spectrum.py`, `opv_sol3.py` implement the Ohta–Percacci–Vacca exponential-
parametrization f(R) flow (arXiv:1511.09393), which admits **exact global quadratic fixed-point
solutions** — sidestepping the singular BVP that defeated the seven single-metric attempts above.

Results (**validated against OPV Table 1**):
- The exact quadratic solutions **identically satisfy** the flow equation (residual = 0, symbolic
  and to machine precision) — confirming the transcription and the fixed point.
- The **critical spectrum**, computed by linearization about the exact solution:
  - Solution 1: **θ = 4, 2.021** (OPV: 4, 2.02) — two relevant directions
  - Solution 3: **θ = 4, 2.413 ± 0.806 i, 0.131** (OPV: 4, 2.4 ± 0.8 i, 0.13) — reproduces the
    complex pair and the near-zero exponent
- ⇒ a **finite** critical surface (two relevant directions), so the fork's UV completion is
  predictive.

This is the critical spectrum the earlier methods could not reach. What remains open is narrower:
matching the fork's *exact* IR coupling ratios to the specific trajectory from this fixed point
(a full flow integration), and the scheme-dependent 2-vs-3 relevant-direction count.

## Update (2026-07-14, cont.): full-flow trajectory test (`opv_traj_fast.py`)

Integrated the FULL OPV flow (numeric, β(g) by polynomial Taylor division) from the exact NGFP
toward the IR, tracking the fork's flow-invariant ratios r_n = g_n g₁ⁿ⁻²/g₂ⁿ⁻¹. Validated:
β(g*)=0 and the flow-Jacobian exponents reproduce (4, 2.05). **Result (partial/honest):** a
trajectory on the 2-dim relevant surface brings the *leading* ratios near the fork's
(r₄≈2.0 to a few %, r₃ to ~20%), but r₅, r₆ diverge (wrong sign / factor 2–4), and the truncated
flow hits a Landau pole (t≈−3.2) before the deep IR. So the fork's low-order structure is roughly
compatible with an AS trajectory, but its exact log-structure ratios are **not** reproduced within
this N=6 truncation; whether the higher-order mismatch is physical or truncation artifact is
unresolved. The critical spectrum (above) is solid; this trajectory test is the honest open edge.
