"""
examples/kdv_nonapplicable.py

Plain KdV equation:  u_t + 6 u u_x + u_xxx = 0.

Traveling-wave reduction (U(Lambda), Lambda = x - c t), integrated once:
    U'' = c U - 3 U^2 - K       =>  a3 = 0.

Since a3 = 0, this equation falls outside the scope of the n=1 GEDF
representation entirely (Theorem 3.1 / eq. 2 require a3 != 0): the
balancing principle gives a different value of n for a quadratic
nonlinearity, so the theory developed in this paper does not apply.

This example demonstrates that the implementation correctly identifies
and rejects an equation lying outside the admissible class, rather
than silently misapplying the trigonometric/hyperbolic case machinery.

Run:  python3 examples/kdv_nonapplicable.py
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gedf_applicability_checker import check_applicability

if __name__ == "__main__":
    print("Plain KdV: u_t + 6 u u_x + u_xxx = 0")
    print("Reduced master ODE: U'' = c U - 3U^2 - K   =>   a3 = 0\n")

    try:
        check_applicability(a0=0, a1=1, a2=-3, a3=0)
    except ValueError as e:
        print(f"Correctly rejected: {e}")
        print("\nThis is expected: the n=1 GEDF ansatz requires a cubic master")
        print("ODE (a3 != 0). Plain KdV reduces to a QUADRATIC master ODE, so")
        print("it falls outside this theory's scope by construction, not due")
        print("to a failure of conditions (C1)-(C4).")
