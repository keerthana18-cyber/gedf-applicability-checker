"""
examples/duffing.py

Undamped, unforced Duffing oscillator:  u'' + alpha u + beta u^3 = 0.

Already in master-ODE form directly (no traveling-wave reduction
needed): a3 = -beta, a2 = 0, a1 = -alpha, a0 = 0.

Checks both sign cases for beta (hardening/softening spring) and
reports which auxiliary-function branch(es) apply in each case.

Run:  python3 examples/duffing.py
"""

import sympy as sp
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gedf_applicability_checker import check_applicability, report

Lambda = sp.symbols("Lambda")

if __name__ == "__main__":
    print("Undamped Duffing oscillator: u'' + alpha u + beta u^3 = 0\n")

    for label, beta, alpha, phi_name, phi_func in [
        ("beta=1 (hardening), alpha=-1, hyperbolic branch", 1, -1, "cosh", sp.cosh),
        ("beta=-1 (softening), alpha=1, trigonometric branch", -1, 1, "cos", sp.cos),
    ]:
        a3 = -beta
        a1 = -alpha
        print(f"--- {label} ---")
        result = check_applicability(a0=0, a1=a1, a2=0, a3=a3)
        report(result)

        if result.applicable:
            H1 = result.H1
            U = H1 / phi_func(Lambda)
            residual = sp.simplify(sp.diff(U, Lambda, 2) + alpha * U + beta * U**3)
            print(f"\nCandidate solution: u(t) = {H1} / {phi_name}(t)  "
                  f"(with Lambda = t here, since this is an ODE, not a PDE reduction)")
            print(f"Symbolic check U'' + alpha U + beta U^3 = {residual}  "
                  f"({'verified' if residual == 0 else 'MISMATCH'})\n")
