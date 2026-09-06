#!/usr/bin/env python3
# Authors: Dr. Davide Batic and Dr. Denys Dutykh
#          (Mathematics Department, Khalifa University of Science and
#           Technology, Abu Dhabi, UAE)
"""Numerical spot check of Eq. (117) for A_{beta,s}.

The defining Bessel-K integral is evaluated directly, without using its
closed-form value.  The result is then compared with the Gamma-product in the
manuscript.  Requires SciPy.
"""

from __future__ import annotations

from math import inf, isfinite, pi, sqrt

from scipy.integrate import quad
from scipy.special import gamma, kv


PAIRS = ((0.0, 0.0), (1.0, 0.0), (0.0, 2.0), (0.5, 1.0))
TOLERANCE = 1.0e-12


def integrand(y: float, q: float) -> float:
    """Return y^(q+2) K_1(y)^2, with its q >= 0 endpoint limit."""
    if y == 0.0:
        return 1.0 if q == 0.0 else 0.0
    return y ** (q + 2.0) * kv(1.0, y) ** 2


def direct_integral(beta: float, s: float) -> tuple[float, float]:
    """Evaluate the defining integral on [0,1] and [1,infinity)."""
    q = beta + s
    if q <= -1.0:
        raise ValueError("The defining integral requires beta + s > -1.")

    opts = {"epsabs": 2.0e-14, "epsrel": 2.0e-14, "limit": 1000}
    low, low_error = quad(integrand, 0.0, 1.0, args=(q,), **opts)
    high, high_error = quad(integrand, 1.0, inf, args=(q,), **opts)
    return low + high, low_error + high_error


def closed_form(beta: float, s: float) -> float:
    """Evaluate the Gamma-product displayed in Eq. (117)."""
    q = beta + s
    return (
        sqrt(pi)
        * gamma((q + 3.0) / 2.0)
        * gamma((q + 5.0) / 2.0)
        * gamma((q + 1.0) / 2.0)
        / (4.0 * gamma((q + 4.0) / 2.0))
    )


def main() -> None:
    print(" beta      s          direct integral          closed form        rel. diff.")
    print("------  ------  ----------------------  ----------------------  -----------")
    maximum_relative_difference = 0.0

    for beta, s in PAIRS:
        direct, error = direct_integral(beta, s)
        closed = closed_form(beta, s)
        relative_difference = abs(direct - closed) / abs(closed)
        maximum_relative_difference = max(
            maximum_relative_difference, relative_difference
        )
        assert isfinite(error)
        print(
            f"{beta:6.2f}  {s:6.2f}  {direct:22.16g}  {closed:22.16g}  "
            f"{relative_difference:11.3e}"
        )

    print(f"\nmaximum relative difference = {maximum_relative_difference:.3e}")
    if maximum_relative_difference >= TOLERANCE:
        raise SystemExit(
            f"Spot check failed: relative difference exceeds {TOLERANCE:.1e}."
        )


if __name__ == "__main__":
    main()
