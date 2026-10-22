# AGENTS.md

Instructions for AI coding agents working in this repository.

## What this is

The archived computation scripts behind the preprint "A Ghost-Free, Unitary, Starobinsky-Class
Realization of Gravity from Entropy" (concept DOI 10.5281/zenodo.21361405). A research record,
not a library: 23 Python scripts in `continuum/`, failed attempts kept on purpose.

## Read in this order

1. `README.md`: the dated correction note under Status, then "How it works" and "Reproduce a
   result".
2. `continuum/README.md`: per-script outcomes, with the dated 2026-07-14 updates at the bottom
   (the later updates supersede the "Honest bottom line" section above them).
3. `continuum/opv_spectrum.py`: the script that produced the critical spectrum.
4. `continuum/opv_traj_conv.py`: the open question (trajectory ratios undecided).

## What you can and cannot see

- You can see every script and its full input: there are no data files. The only run-time link
  between scripts is `from cheb import map_to` (in `fR_spectral_fp.py`, `fR_spectral_v2.py`,
  `fR_spectral_seeded.py`, `fR_homotopy.py`). Seeds and fork ratios were copied into sources by
  hand; see the second diagram in README.md.
- Scripts that save `continuum/*.npy` write files that no script reads and that are not
  committed.
- You cannot see the Cadabra2 ghost-freedom algebra, the hi_class cosmology runs, or the preprint
  itself. Those are not in this repository; the preprint and its LaTeX are on Zenodo.

## The only numbers to trust

Each is printed by the script named; quote it from a run, not from memory.

| number | source |
|---|---|
| fork ratios r_3..r_6 = 4/3, 2, 16/5, 16/3 | `continuum/fork_embedding.py` |
| Solution 1 residual identically zero | `continuum/opv_verify.py` |
| Solution 1 exponents 4.000, 2.021 (2 relevant) | `continuum/opv_spectrum.py` |
| Solution 3 exponents 4.000, 2.413 +/- 0.806i, 0.131 | `continuum/opv_sol3.py` |
| theta_2 by truncation N = 6, 8, 10, 12 | `continuum/opv_traj_conv.py` |

Anything else (g\*λ\* bands, "agent-verified" literature values printed by `ngfp_mechanism.py`
or `fork_embedding.py`) is either from a slow script not rerun on 2026-10-22 or quoted from the
literature, and must not be restated as a result of this repository.

## Setup

```bash
python -m pip install numpy scipy sympy mpmath
```

Python 3.12 (tested 3.12.10). On Windows set `PYTHONIOENCODING=utf-8` before running any script.
Run every script from the repository root (the `cheb` import and the `.npy` paths assume it).

## Test

There is no test suite. The fast checks, run from the repository root, are:

```bash
python continuum/cheb.py          # must print a line containing "=> PASS"
python continuum/opv_verify.py    # must print "IS IDENTICALLY ZERO: True"
python continuum/opv_spectrum.py  # must print "2 relevant; leading = +4.000"
```

## Conventions

- numpy, scipy, sympy and mpmath only; do not add dependencies.
- Each script is standalone and prints its own result; keep it that way.
- No em-dashes in docs.

## Do not

- Delete or "fix" the failed attempts (`fR_*`, `lpa_*`). They are archived on purpose as part of
  the published record.
- Change a reported number in a README without rerunning the script that produces it, and record
  a changed published claim as a dated correction note, not a silent edit.
- Commit `.npy` outputs, secrets, hostnames, internal paths or personal emails.

## Done means

The three checks above print the expected lines, and if README.md or continuum/README.md changed,
`llms-full.txt` is regenerated from them in the same change.
