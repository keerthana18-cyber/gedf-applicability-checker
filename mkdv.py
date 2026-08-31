"""
examples/mkdv.py

Modified KdV equation:  u_t + 6 u^2 u_x + u_xxx = 0.

Traveling-wave reduction (U(Lambda), Lambda = x - c t), integrated once:
    U'' = c U - 2 U^3 - K        =>  a3 = -2, a2 = 0, a1 = c, a0 = -K.

This example checks applicability for the natural choice c = 1, K = 0,
which is admissible (hyperbolic branch, a3 < 0), and reports the
resulting exact soliton solution.

Run:  python3 examples/mkdv.py
"""

import sympy as sp
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gedf_applicability_checker import check_applicability, report

Lambda = sp.symbols("Lambda")

if __name__ == "__main__":
    print("Modified KdV: u_t + 6 u^2 u_x + u_xxx = 0")
    print("Reduced master ODE (c=1, K=0): U'' = U - 2U^3\n")

    result = check_applicability(a0=0, a1=1, a2=0, a3=-2)
    report(result)

    if result.applicable:
        H1 = result.H1
        U = H1 / sp.cosh(Lambda)
        residual = sp.simplify(sp.diff(U, Lambda, 2) - (U - 2 * U**3))
        print(f"\nCandidate solution: u(x,t) = {H1} * sech(x - t)")
        print(f"Symbolic check U'' - (U - 2U^3) = {residual}  "
              f"({'verified' if residual == 0 else 'MISMATCH'})")
