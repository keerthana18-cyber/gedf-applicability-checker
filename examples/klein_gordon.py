"""
examples/klein_gordon.py

Nonlinear Klein-Gordon (phi^4-type) equation:  u_tt - u_xx - u + u^3 = 0.

Traveling-wave reduction (U(Lambda), Lambda = x - c t):
    (c^2 - 1) U'' = U - U^3
    =>  a3 = -1/(c^2-1), a1 = 1/(c^2-1), a2 = 0, a0 = 0.

Here a3 and a1 both depend on the wave speed c, so the required c is
solved for directly (condition C4: a1 = +1 for the hyperbolic branch)
rather than passed in as a fixed coefficient.

Run:  python3 examples/klein_gordon.py
"""

import sympy as sp
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gedf_applicability_checker import check_applicability, report

Lambda, c = sp.symbols("Lambda c", positive=True)

if __name__ == "__main__":
    print("Nonlinear Klein-Gordon (phi^4-type): u_tt - u_xx - u + u^3 = 0")

    a1_expr = 1 / (c**2 - 1)
    c_solutions = sp.solve(sp.Eq(a1_expr, 1), c)  # (C4) for hyperbolic branch: a1 = +1
    print(f"Solving a1 = 1 (hyperbolic branch) for c: c = {c_solutions}\n")

    for c_val in c_solutions:
        a3_val = sp.simplify((-1 / (c**2 - 1)).subs(c, c_val))
        a1_val = sp.simplify(a1_expr.subs(c, c_val))
        print(f"--- c = {c_val} ---")
        result = check_applicability(a0=0, a1=a1_val, a2=0, a3=a3_val)
        report(result)

        if result.applicable:
            H1 = result.H1
            U = H1 / sp.cosh(Lambda)
            residual = sp.simplify((c_val**2 - 1) * sp.diff(U, Lambda, 2) - (U - U**3))
            print(f"\nCandidate solution: u(x,t) = {H1} * sech(x - ({c_val}) t)")
            print(f"Symbolic check (c^2-1)U'' - (U - U^3) = {residual}  "
                  f"({'verified' if residual == 0 else 'MISMATCH'})\n")
