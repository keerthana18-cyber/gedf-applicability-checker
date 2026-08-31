# GEDF Applicability Checker

This repository accompanies the manuscript currently under review:
"An Algebraic and Geometric Applicability Theory for the Generalized
Exponential Differential Function Method" (N. Keerthana, V. Umapathi,
N. Annapoorani).

A symbolic pre-check for the applicability of the Generalized Exponential
Differential Function (GEDF) method to the cubic master ODE

    u'' = a3*u^3 + a2*u^2 + a1*u + a0,   a3 != 0,

implementing conditions (C1)-(C4) of Theorem 4.1, the factorization
criterion of Theorem 5.1, and the integration-constant recovery of
Theorem 7.1 from the accompanying paper.

## Requirements

- Python 3.8+
- SymPy (`pip install sympy`)

## Usage

```python
from gedf_applicability_checker import check_applicability, report

# WBK-BK example (Section 8.1), R0 = 1
result = check_applicability(a0=0.5, a1=0.5, a2=-1.5, a3=0.5)
report(result)
```

Coefficients may be Python floats, ints, `Fraction`/`Rational`, or fully
symbolic SymPy expressions (see Example 4 in `gedf_applicability_checker.py` for a
symbolic check with `R0` left free).

Running the script directly reproduces all four examples discussed in
the paper:

```bash
python3 gedf_applicability_checker.py
```

## What it checks

Given `(a0, a1, a2, a3)`:

1. Computes `R0 = -a2/(3*a3)` (from condition C2).
2. Checks `f(R0) = 0` (condition C1).
3. Matches `f'(R0)` against `delta_phi in {-1, +1}` to identify the
   trigonometric or hyperbolic branch (condition C4).
4. Checks the admissibility sign `eps_phi * a3 > 0` and computes `H1`
   (condition C3).
5. If (C1) fails but (C2)-(C4) hold, reports the unique integration
   constant `a0 = C*` (Theorem 7.1) that restores applicability.
6. If (C4) fails for every branch, reports the value of `a1` that
   would need to hold instead.
7. Cross-checks the result against the factorization
   `f(u) = a3(u-R0)[(u-R0)^2 - mu^2]` (Theorem 5.1) when applicable.

## Additional examples

The `examples/` directory contains standalone, runnable checks against
several other well-known nonlinear evolution equations and nonlinear
oscillators from the traveling-wave literature, beyond the two detailed
applications (WBK-BK, KD) presented in the paper:

- `examples/mkdv.py` — modified KdV equation (hyperbolic branch, bounded sech soliton)
- `examples/gardner.py` — Gardner equation (hyperbolic branch, bounded sech soliton)
- `examples/klein_gordon.py` — nonlinear Klein-Gordon / phi^4-type equation (hyperbolic branch, bounded sech soliton)
- `examples/duffing.py` — undamped Duffing oscillator, both hardening and softening cases
- `examples/kdv_nonapplicable.py` — plain KdV equation, demonstrating correct **rejection**: its traveling-wave reduction gives `a3 = 0`, placing it outside the scope of the n=1 GEDF ansatz entirely (rather than a failure of conditions (C1)-(C4))

Each example independently verifies its candidate solution by direct
symbolic substitution back into the governing ODE. Run any of them with,
e.g.,

```bash
python3 examples/mkdv.py
```

## Citation

If you use this software, please cite the Zenodo record:

    <Zenodo DOI — to be added once archived>

Citation information for the associated paper will be added once it
is published.

## License

MIT License — see `LICENSE`.
