#!/usr/bin/env python3
# Authors: Dr. Davide Batic and Dr. Denys Dutykh
#          (Mathematics Department, Khalifa University of Science and
#           Technology, Abu Dhabi, UAE)
"""Exact proper time of the E=1 Painleve-Gullstrand infalling congruence
versus the endpoint expression obtained by replacing the integral of a
nonconstant function by an endpoint value.

For m(r) = M r^3 / (r^2 + l0^2)^{3/2}, the exact proper time of the marginal
radial congruence is

    tau_exact(r; r_i) = (l0^{3/2}/sqrt(2M)) [ I(u_i) - I(u) ],
    u = (1 + r^2/l0^2)^{1/4},
    I(u) = (2/3) u^3 + arctan(u) + (1/2) ln((u-1)/(u+1)),

which diverges only logarithmically as r -> 0, whereas the endpoint
prescription gives the spurious power law tau_JN ~ r^{-2}.

Parameters M = 5, l0 = 1, r_i = 4 l0, matching the original comparison.
"""

import numpy as np

import ku_style as ks
from ku_style import plt

M, L0, XI = 5.0, 1.0, 4.0


def Theta(r):
    return -3.0 * np.sqrt(2.0 * M) * (r**2 + 2.0 * L0**2) / (
        2.0 * (r**2 + L0**2) ** 1.75
    )


def Ffun(r):
    return r**2 * (3.0 * r**2 + 10.0 * L0**2) / (3.0 * (r**2 + 2.0 * L0**2) ** 2)


def Ical(u):
    return (2.0 / 3.0) * u**3 + np.arctan(u) + 0.5 * np.log((u - 1.0) / (u + 1.0))


def tau_exact(r):
    u = (1.0 + r**2 / L0**2) ** 0.25
    ui = (1.0 + XI**2 / L0**2) ** 0.25
    return (L0**1.5 / np.sqrt(2.0 * M)) * (Ical(ui) - Ical(u))


def tau_JN(r):
    return (1.0 / Theta(r) - 1.0 / Theta(XI)) / Ffun(r)


def main():
    slope = L0**1.5 / np.sqrt(2.0 * M)  # asymptotic slope of tau_exact vs ln(1/x)

    # Check the logarithmic asymptotics: d(tau_exact)/d ln(1/x) -> slope as x->0.
    slope_est = (tau_exact(1e-5) - tau_exact(1e-4)) / np.log(10.0)
    assert abs(slope_est - slope) / slope < 1e-2, "log slope mismatch"

    x = np.logspace(np.log10(3e-3), np.log10(0.985 * XI), 1400)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.4, 4.5))

    # Left: log-log comparison.
    axL.loglog(x, tau_exact(x), color=ks.NAVY, label=r"$\tau_{\rm exact}$")
    axL.loglog(x, tau_JN(x), color=ks.CORAL, ls="--", label=r"$\tau_{\rm JN}$")
    axL.set_xlabel(r"$x=r/\ell_0$")
    axL.set_ylabel(r"proper time")
    axL.annotate(r"$\tau_{\rm exact}\sim\ln(1/x)$", (5e-3, tau_exact(5e-3)),
                 xytext=(1.4e-2, 0.32 * tau_exact(5e-3)), color=ks.NAVY, fontsize=10)
    xJ = 2.0e-2
    axL.text(xJ, 1.8 * tau_JN(xJ), r"$\tau_{\rm JN}\sim x^{-2}$", color=ks.CORAL,
             fontsize=10, ha="left", va="bottom")
    ks.grid(axL)
    ks.legend(axL, loc="upper right")

    # Right: exact proper time is asymptotically linear in ln(1/x).
    xi = np.log(1.0 / x)
    axR.plot(xi, tau_exact(x), color=ks.NAVY, label=r"$\tau_{\rm exact}$")
    ref = slope * (xi + np.log(XI))  # = slope * ln(x_i/x)
    axR.plot(xi, ref, color=ks.GREY, ls=":", lw=1.2,
             label=r"$\frac{\ell_0^{3/2}}{\sqrt{2M}}\ln(x_i/x)$")
    axR.set_xlabel(r"$\ln(1/x)$")
    axR.set_ylabel(r"$\tau_{\rm exact}$")
    ks.grid(axR)
    ks.legend(axR, loc="upper left")

    fig.tight_layout()
    ks.savefig(fig, "proper_time_comparison")


if __name__ == "__main__":
    main()
