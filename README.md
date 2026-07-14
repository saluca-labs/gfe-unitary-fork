# Gravity from Entropy — the unitary fork (computation code)

Computation code accompanying the preprint **_A Ghost-Free, Unitary, Starobinsky-Class
Realization of Gravity from Entropy_** (Cristian Ruvalcaba and the Saluca Agentic AI Research
Team, 2026).

The paper constructs an explicit scalar-field (Fisher-metric) realization of Ginestra Bianconi's
_Gravity from Entropy_ program and shows it is ghost-free (three ways), unitary (at the
field-theory and operator-emergence levels, with the non-perturbative KMS-vs-bath question
settled), and observationally viable (Starobinsky inflation + thawing dark energy). This
repository holds the symbolic and numerical computations behind those claims — **including the
attempts that failed.**

> 📄 Preprint + full archive (PDF, LaTeX source, this code): Zenodo
> [`10.5281/zenodo.21361406`](https://doi.org/10.5281/zenodo.21361406) *(DOI resolves once the
> record is published).*

## What's here

All scripts are in [`continuum/`](continuum/) with a per-script guide in
[`continuum/README.md`](continuum/README.md). Python 3.12 (`numpy`, `sympy`, `scipy`, `mpmath`).
On Windows set `PYTHONIOENCODING=utf-8`. Run from the repository root.

| area | scripts | status |
|---|---|---|
| Fork in RG theory space (radius of convergence = R_max, one-parameter ray, UV asymptotics) | `fork_embedding.py` | self-derived, clean |
| Scalar-sector UV analysis (disformal & k-essence irrelevant ⇒ UV-benign) | `scalar_sector.py` | clean |
| Asymptotic-safety mechanism | `ngfp_mechanism.py` | schematic, exact-in-toy |
| Chebyshev spectral instrument | `cheb.py` | **validated to machine precision** |
| f(R) non-Gaussian fixed point (exact Dietz–Morris equation) | `fR_fixedpoint.py`, `fR_improved.py`, `fR_spectral_seeded.py` | **reproduced: g\*λ\* ≈ 0.09–0.1, in the CPR band, confirmed four independent ways** |
| Critical spectrum / trajectory | `fR_exponents_fork.py`, `fR_spectral_fp.py`, `fR_spectral_v2.py` | **open** — see below |

## Honest status

The f(R) non-Gaussian fixed point is **reproduced and robustly confirmed**. The fork's _exact_
critical spectrum and RG trajectory are **not** obtained here: five principled numerical attempts
each failed for a well-understood reason (small-propagator-denominator instability in polynomial
truncations; a singular boundary-value problem that naive Newton / least-squares cannot converge
without parameter continuation), and a verified search found **no turn-key public code** for the
f(R) fixed-point spectrum — the field runs private Mathematica. That computation is well-posed
but externally gated; it is the honest open edge stated in §8 of the paper. **No fixed point,
exponent, or trajectory was ever fabricated to close a gap.**

## Provenance / AI transparency

This work was produced by a human-directed agentic large-language-model research system under
continuous human validation, with a human author solely accountable for every claim. The
apparatus repeatedly caught and retracted its own overreach (hallucinated arXiv identifiers,
incorrect critical exponents) before it could stand; the failed attempts are archived here by
design. No claim should be treated as established solely because an AI produced it.

## License

Code: [MIT](LICENSE). The preprint text is archived separately under CC-BY-4.0 (see the Zenodo
record).
