# AGENTS.md

Instructions for AI coding agents working in this repository.

## Setup

```bash
python -m pip install numpy scipy sympy mpmath
```

Python 3.12. On Windows set `PYTHONIOENCODING=utf-8` before running any script.

## Test

There is no test suite. The fast checks, run from the repository root, are:

```bash
python continuum/cheb.py          # must end with "=> PASS"
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
- Commit secrets, hostnames, internal paths or personal emails.

## Done means

The three checks above print the expected lines, and if README.md or continuum/README.md changed,
`llms-full.txt` is regenerated from them in the same change.
