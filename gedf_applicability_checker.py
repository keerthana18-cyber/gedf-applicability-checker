"""
gedf_applicability_checker.py
==============================

A symbolic pre-check for the applicability of the Generalized Exponential
Differential Function (GEDF) method to the cubic master ODE

    u'' = f(u) = a3*u^3 + a2*u^2 + a1*u + a0,   a3 != 0.

Implements conditions (C1)-(C4) of Theorem 4.1, the factorization
criterion of Theorem 5.1, and the integration-constant recovery of
Theorem 7.1, for a cubic master ODE obtained from a traveling-wave
reduction of a nonlinear evolution equation.

No numerical algebraic system needs to be solved by the user: given the
four coefficients (a0, a1, a2, a3), the script reports whether the GEDF
representation

    u(Lambda) = R0 + H1 / phi(Lambda)

is an exact solution, for phi in {cos, -sin, sinh, cosh}, and if not,
diagnoses which condition fails and, where possible, the minimal
coefficient adjustment (following Theorem 7.1) that would restore
applicability.

Usage
-----
    from gedf_check import check_applicability, report

    result = check_applicability(a0=0.5, a1=0.5, a2=-1.5, a3=0.5)
    report(result)

Coefficients may be Python floats/ints or SymPy expressions (for a
fully symbolic check, e.g. with a3, a2, a1, a0 left as free symbols
subject to R0-based relations).

Author's note: this implementation accompanies the paper "A
Mathematical Characterization of the Applicability of the Generalized
Exponential Differential Function Method" and is provided so that
Theorem 4.1's applicability test can be applied to a candidate cubic
master ODE without hand computation.
"""

from dataclasses import dataclass, field
from typing import Optional, List
import sympy as sp

u = sp.symbols('u')

# Auxiliary-function branch data, from Table 2 of the paper.
# eps_phi, delta_phi, admissibility sign on a3, and a human label.
BRANCHES = [
    {"name": "cos(Lambda)",   "eps": 1,  "delta": -1, "type": "trigonometric"},
    {"name": "-sin(Lambda)",  "eps": 1,  "delta": -1, "type": "trigonometric"},
    {"name": "sinh(Lambda)",  "eps": 1,  "delta": 1,  "type": "hyperbolic"},
    {"name": "cosh(Lambda)",  "eps": -1, "delta": 1,  "type": "hyperbolic"},
]


@dataclass
class GEDFResult:
    a0: object
    a1: object
    a2: object
    a3: object
    R0: object
    fR0: object                      # f(R0), should be 0 for (C1)
    fprime_R0: object                # f'(R0)
    C1: bool
    C2_R0: object                    # R0 as implied by (C2); should equal R0 above
    C2: bool
    matched_branch: Optional[dict] = None
    H1: Optional[object] = None
    C3: bool = False
    C4: bool = False
    applicable: bool = False
    factorization_mu2: Optional[object] = None
    diagnosis: List[str] = field(default_factory=list)
    suggested_a0: Optional[object] = None   # Theorem 7.1 recovery
    suggested_a1: Optional[object] = None   # coefficient needed for (C4)


def _f(a0, a1, a2, a3, x):
    return a3 * x**3 + a2 * x**2 + a1 * x + a0


def check_applicability(a0, a1, a2, a3, numeric_tol=1e-9) -> GEDFResult:
    """
    Check conditions (C1)-(C4) of Theorem 4.1 for the cubic master ODE
    u'' = a3 u^3 + a2 u^2 + a1 u + a0.

    Returns a GEDFResult with a full diagnostic trail. Works with both
    exact SymPy values and floats; floats are compared to zero with
    `numeric_tol`.
    """
    a0, a1, a2, a3 = sp.sympify(a0), sp.sympify(a1), sp.sympify(a2), sp.sympify(a3)
    if a3 == 0:
        raise ValueError("a3 = 0: the master equation is not cubic; GEDF (n=1) does not apply.")

    def is_zero(expr):
        expr = sp.nsimplify(sp.simplify(expr))
        if expr.free_symbols:
            return sp.simplify(expr) == 0
        return abs(float(expr)) < numeric_tol

    fexpr = _f(a0, a1, a2, a3, u)
    fprime = sp.diff(fexpr, u)

    # (C2) determines R0 = -a2/(3 a3); this is the only candidate R0
    # consistent with the reduced ansatz (Theorem 3.1 / eq. 4).
    R0 = sp.simplify(-a2 / (3 * a3))

    fR0 = sp.simplify(fexpr.subs(u, R0))
    fprime_R0 = sp.simplify(fprime.subs(u, R0))

    C2_ok = True  # R0 was *defined* via (C2), so (C2) holds by construction
    C1_ok = is_zero(fR0)

    diagnosis = []
    if not C1_ok:
        diagnosis.append(
            f"(C1) fails: f(R0) = {fR0} != 0. The integration constant does not "
            f"currently place the equation on the applicability manifold."
        )

    # Try to match a branch via (C4): f'(R0) == delta_phi. Since sinh and
    # cosh share delta_phi = +1 (and cos/-sin share delta_phi = -1), more
    # than one branch can satisfy (C4); among those, prefer one whose
    # admissibility sign eps_phi * a3 > 0 also holds (C3), since that is
    # the branch that actually produces a real H1. Fall back to the first
    # delta-match (for diagnostic purposes) only if none are admissible.
    delta_matches = [br for br in BRANCHES if is_zero(fprime_R0 - br["delta"])]
    matched = None
    for br in delta_matches:
        try:
            admissible = float(br["eps"] * a3) > 0 if not a3.free_symbols else None
        except TypeError:
            admissible = None
        if admissible is not False:
            matched = br
            break
    if matched is None and delta_matches:
        matched = delta_matches[0]

    result = GEDFResult(
        a0=a0, a1=a1, a2=a2, a3=a3, R0=R0, fR0=fR0, fprime_R0=fprime_R0,
        C1=C1_ok, C2_R0=R0, C2=C2_ok, diagnosis=diagnosis,
    )

    if matched is None:
        result.diagnosis.append(
            f"(C4) fails for every auxiliary function: f'(R0) = {fprime_R0}, "
            f"which matches none of delta_phi in {{-1 (trig), +1 (hyperbolic)}}."
        )
        # Suggest the coefficient a1 needed to satisfy (C4) for the
        # trigonometric branch (the more common case, a3 > 0), per
        # Theorem 7.1's use of a1^(eff).
        a1_eff_trig = sp.simplify(3 * a3 * R0**2 - 1)
        a1_eff_hyp = sp.simplify(3 * a3 * R0**2 + 1)
        result.suggested_a1 = {"trigonometric": a1_eff_trig, "hyperbolic": a1_eff_hyp}
        result.diagnosis.append(
            f"Suggested a1 to satisfy (C4): {a1_eff_trig} (trigonometric branch, "
            f"needs a3 > 0) or {a1_eff_hyp} (hyperbolic branch, needs a3 < 0)."
        )
        result.applicable = False
        return result

    result.matched_branch = matched
    result.C4 = True

    # (C3): H1^2 = eps_phi * 2 / a3, subject to eps_phi * a3 > 0
    eps = matched["eps"]
    admissible_sign = None
    try:
        admissible_sign = float(eps * a3) > 0 if not a3.free_symbols else None
    except TypeError:
        admissible_sign = None

    H1_squared = sp.simplify(eps * 2 / a3)
    if admissible_sign is False:
        result.diagnosis.append(
            f"(C3) admissibility fails for {matched['name']}: requires "
            f"eps_phi * a3 > 0, i.e. a3 {'> 0' if eps == 1 else '< 0'}."
        )
        result.C3 = False
        result.applicable = False
        return result

    result.C3 = True
    result.H1 = sp.sqrt(H1_squared)

    result.applicable = bool(C1_ok)  # C2, C3, C4 already confirmed by this point

    if not C1_ok:
        # Theorem 7.1, Part (i): with (C2) and (C4) already satisfied,
        # a unique a0 = C* restores (C1).
        a1_eff = 3 * a3 * R0**2 - 1 if matched["type"] == "trigonometric" else 3 * a3 * R0**2 + 1
        if is_zero(a1 - a1_eff):
            C_star = sp.simplify(-(a3 * R0**3 + a2 * R0**2 + a1 * R0))
            result.suggested_a0 = C_star
            result.diagnosis.append(
                f"(C2) and (C4) already hold for the {matched['type']} branch; "
                f"Theorem 7.1 gives the required integration constant "
                f"a0 = C* = {C_star} to satisfy (C1)."
            )
        else:
            result.diagnosis.append(
                f"(C4) uses a1 = {a1}, but the {matched['type']} branch actually "
                f"requires a1 = {a1_eff} for f'(R0) to equal {matched['delta']}. "
                f"The match above was against a different branch's delta; "
                f"re-check which auxiliary function is intended before applying "
                f"Theorem 7.1's a0-only recovery."
            )
    else:
        # Full factorization check (Theorem 5.1), as a redundant cross-check.
        mu2 = sp.simplify(1 / a3) if matched["type"] == "trigonometric" else sp.simplify(-1 / a3)
        result.factorization_mu2 = mu2
        refactored = sp.expand(a3 * (u - R0) * ((u - R0)**2 - mu2))
        if sp.simplify(refactored - fexpr) == 0:
            result.diagnosis.append(
                f"Factorization confirmed: f(u) = a3 (u-R0)[(u-R0)^2 - mu^2] "
                f"with mu^2 = {mu2}."
            )

    return result


def report(result: GEDFResult):
    print("=" * 70)
    print(f"Coefficients: a3={result.a3}, a2={result.a2}, a1={result.a1}, a0={result.a0}")
    print(f"R0 (from C2)      : {result.R0}")
    print(f"f(R0)             : {result.fR0}   -> (C1) {'holds' if result.C1 else 'FAILS'}")
    print(f"f'(R0)            : {result.fprime_R0}")
    if result.matched_branch:
        print(f"Matched branch    : {result.matched_branch['name']} "
              f"({result.matched_branch['type']}) -> (C4) holds")
        print(f"(C3) admissibility: {'holds' if result.C3 else 'FAILS'}")
        if result.H1 is not None:
            print(f"H1                : +/- {result.H1}")
    else:
        print("(C4) FAILS for every auxiliary function.")
    print("-" * 70)
    print(f"APPLICABLE: {result.applicable}")
    if result.diagnosis:
        print("Diagnosis:")
        for line in result.diagnosis:
            print(f"  - {line}")
    print("=" * 70)


if __name__ == "__main__":
    R0sym = sp.symbols('R0', positive=True)

    print("\n### Example 1: WBK-BK system (Section 8.1), R0 = 1")
    R0 = 1
    a3 = sp.Rational(1, 2) / R0**2
    a2 = sp.Rational(-3, 2) / R0
    a1 = sp.Rational(1, 2)
    a0 = sp.Rational(1, 2) * R0
    res = check_applicability(a0, a1, a2, a3)
    report(res)

    print("\n### Example 2: KD system, C = 0 (Section 8.2, 'Failure' subsection)")
    a, m, b = 1, 0, sp.Rational(2, 5)
    a3 = sp.Rational(a**2, 2)
    a2 = sp.Rational(3*a*m - 6*b, 2)
    c = sp.Rational(64, 100)     # c = 0.64, from eq. (waveconstraint)
    a1 = c - 3*m**2
    a0 = 0
    res = check_applicability(a0, a1, a2, a3)
    report(res)

    print("\n### Example 3: KD system, corrected recovery (a1 = -0.04, a0 = C* = 0.544)")
    a3 = sp.Rational(1, 2)
    a2 = sp.Rational(-6, 5)
    a1 = sp.Rational(-4, 100)
    a0 = sp.Rational(544, 1000)
    res = check_applicability(a0, a1, a2, a3)
    report(res)

    print("\n### Example 4: fully symbolic check, R0 left free (trigonometric branch)")
    R0 = sp.symbols('R0')
    a3s = sp.symbols('a3', positive=True)
    a2s = -3 * a3s * R0
    a1s = 3 * a3s * R0**2 - 1
    a0s = -a3s * R0**3 + R0
    res = check_applicability(a0s, a1s, a2s, a3s)
    report(res)
