# Changelog

## Version 1.0

- Corrected branch-selection logic for hyperbolic auxiliary functions:
  previously, when both `sinh` and `cosh` matched condition (C4) via
  the same `delta_phi = +1`, the checker always selected `sinh`
  regardless of whether its admissibility sign (`eps_phi * a3 > 0`)
  actually held, incorrectly reporting non-applicability in several
  genuinely applicable cases. The checker now selects the delta-matched
  branch whose admissibility condition also holds.
- Added standalone worked examples (`examples/`): modified KdV, Gardner,
  nonlinear Klein-Gordon, Duffing oscillator, and a non-applicable case
  (plain KdV).
- Added symbolic verification of every example's candidate solution by
  direct substitution back into the governing ODE.
