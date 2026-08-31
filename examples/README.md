# Examples

Standalone, runnable checks of the GEDF applicability test against
several nonlinear evolution equations and nonlinear oscillators from
the traveling-wave literature, beyond the two detailed applications
(WBK-BK, KD) presented in the accompanying paper.

Each script is self-contained: it performs (or states) the
traveling-wave reduction to the master ODE, calls
`check_applicability` from the parent `gedf_applicability_checker.py`,
and independently verifies the resulting candidate solution by direct
symbolic substitution back into the governing ODE.

| File | Equation | Outcome |
|---|---|---|
| `mkdv.py` | Modified KdV | Applicable — hyperbolic branch, bounded sech soliton |
| `gardner.py` | Gardner equation | Applicable — hyperbolic branch, bounded sech soliton |
| `klein_gordon.py` | Nonlinear Klein-Gordon (phi^4-type) | Applicable — hyperbolic branch, bounded sech soliton |
| `duffing.py` | Undamped Duffing oscillator | Applicable in both the hardening and softening cases |
| `kdv_nonapplicable.py` | Plain KdV | **Not applicable** — traveling-wave reduction gives `a3 = 0`, outside the scope of the n=1 GEDF ansatz by construction |

## Running an example

```bash
python3 mkdv.py
```

(run from inside this directory, or adjust `sys.path` as shown at the
top of each file if running from elsewhere.)

## Why `kdv_nonapplicable.py` is included

It demonstrates that the implementation correctly identifies and
rejects an equation lying outside the admissible class (a quadratic,
not cubic, master ODE), rather than only ever reporting success.
