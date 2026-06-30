#!/usr/bin/env python3
# Authors: Dr. Davide Batic and Dr. Denys Dutykh
#          (Mathematics Department, Khalifa University of Science and
#           Technology, Abu Dhabi, UAE)
"""Effective potentials for causal geodesics of the T-duality-inspired
Bardeen black hole.

With x = r/l0, mu = M/l0, lambda = L/l0, and
    f(x) = 1 - 2 mu x^2 / (x^2 + 1)^{3/2},
the radial first integral is x_dot^2 = E^2 - V(x), with

    timelike:  V_t(x) = f(x) (1 + lambda^2 / x^2),
    null:      V_n(x) = f(x) lambda^2 / x^2.

For lambda != 0 the centrifugal term ~ lambda^2/x^2 produces a divergent
barrier as x -> 0, so non-radial causal geodesics cannot reach the centre.
"""

import numpy as np

import ku_style as ks
from ku_style import plt

MU = 2.0  # black-hole branch (mu > 3 sqrt(3)/4 ~ 1.299)


def f(x):
    return 1.0 - 2.0 * MU * x**2 / (x**2 + 1.0) ** 1.5


def V_t(x, lam):
    return f(x) * (1.0 + lam**2 / x**2)


def V_n(x, lam):
    return f(x) * lam**2 / x**2


def main():
    x = np.linspace(0.06, 5.0, 2000)
    colors = {0: ks.NAVY, 1: ks.TEAL, 2: ks.CORAL}

    fig, (axT, axN) = plt.subplots(1, 2, figsize=(10.4, 4.5), sharex=True)

    # Timelike panel.
    for lam in (0, 1, 2):
        axT.plot(x, V_t(x, lam), color=colors[lam],
                 label=rf"$\lambda={lam}$")
    axT.axhline(1.0, color=ks.GREY, ls=":", lw=1.0)
    axT.text(3.3, 1.04, r"$E^2=1$", color=ks.GREY, fontsize=9)
    axT.axhline(0.0, color=ks.GREY, lw=0.8, alpha=0.6)
    axT.set_title(r"timelike: $V_t=f\,(1+\lambda^2/x^2)$")
    axT.set_xlabel(r"$x=r/\ell_0$")
    axT.set_ylabel(r"$V_t(x)$")
    axT.set_xlim(0, 5)
    axT.set_ylim(-2.4, 4.0)
    ks.grid(axT)
    ks.legend(axT, loc="upper right")

    # Null panel.
    for lam in (1, 2):
        axN.plot(x, V_n(x, lam), color=colors[lam], label=rf"$\lambda={lam}$")
    axN.axhline(0.0, color=ks.GREY, lw=0.8, alpha=0.6)
    axN.set_title(r"null: $V_n=f\,\lambda^2/x^2$")
    axN.set_xlabel(r"$x=r/\ell_0$")
    axN.set_ylabel(r"$V_n(x)$")
    axN.set_xlim(0, 5)
    axN.set_ylim(-2.4, 4.0)
    axN.annotate(r"$\sim\lambda^2/x^2$ barrier", (0.5, 3.3), fontsize=9,
                 color=ks.GREY)
    ks.grid(axN)
    ks.legend(axN, loc="upper right")

    fig.tight_layout()
    ks.savefig(fig, "geodesic_potentials")


if __name__ == "__main__":
    main()
