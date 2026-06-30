#!/usr/bin/env python3
# Authors: Dr. Davide Batic and Dr. Denys Dutykh
#          (Mathematics Department, Khalifa University of Science and
#           Technology, Abu Dhabi, UAE)
"""Horizon phase diagram and surface gravities of the T-duality-inspired
Bardeen black hole.

Left panel: the horizon function mu = H(x) = (x^2+1)^{3/2}/(2 x^2), with a
minimum mu_ext = 3 sqrt(3)/4 at x = sqrt(2). For mu > mu_ext a horizontal line
meets the curve at the inner (x_-) and outer (x_+) horizons.

Right panel: the dimensionless surface gravities
    l0 kappa_+ = (x_+^2 - 2)/(2 x_+ (x_+^2 + 1)),
    l0 kappa_- = (2 - x_-^2)/(2 x_- (x_-^2 + 1)),
both vanishing at extremality; kappa_- grows ~ sqrt(2 mu) at large mass while
kappa_+ ~ 1/(4 mu).
"""

import numpy as np
from scipy.optimize import brentq

import ku_style as ks
from ku_style import plt


def H(x):
    return (x**2 + 1.0) ** 1.5 / (2.0 * x**2)


def kappa(xh):
    return abs(xh**2 - 2.0) / (2.0 * xh * (xh**2 + 1.0))


def horizons(mu):
    xm = brentq(lambda x: H(x) - mu, 1e-6, np.sqrt(2.0))
    xp = brentq(lambda x: H(x) - mu, np.sqrt(2.0), 1e4)
    return xm, xp


def main():
    mu_ext = 3.0 * np.sqrt(3.0) / 4.0

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.4, 4.5))

    # Left: horizon phase diagram.
    x = np.linspace(0.35, 6.0, 2000)
    axL.plot(x, H(x), color=ks.NAVY)
    axL.axhline(mu_ext, color=ks.GREY, ls="--", lw=0.8)
    axL.plot([np.sqrt(2.0)], [mu_ext], "o", color=ks.CORAL, ms=5)
    mu_demo = 2.0
    xm, xp = horizons(mu_demo)
    axL.axhline(mu_demo, color=ks.TEAL, ls=":", lw=1.0)
    axL.plot([xm, xp], [mu_demo, mu_demo], "o", color=ks.TEAL, ms=5)
    axL.text(0.42, mu_ext - 0.06, r"$\mu_{\rm ext}=3\sqrt{3}/4$", color=ks.GREY,
             fontsize=10, va="top")
    axL.annotate(r"$x_-$", (xm, mu_demo), xytext=(xm + 0.08, mu_demo + 0.05),
                 fontsize=10, color=ks.TEAL, ha="left", va="bottom")
    axL.annotate(r"$x_+$", (xp, mu_demo), xytext=(xp - 0.08, mu_demo + 0.05),
                 fontsize=10, color=ks.TEAL, ha="right", va="bottom")
    axL.set_xlabel(r"$x=r/\ell_0$")
    axL.set_ylabel(r"$\mu=M/\ell_0$")
    axL.set_xlim(0.25, 6)
    axL.set_ylim(1.05, 3.2)
    axL.set_title(r"horizon phase diagram")
    ks.grid(axL)

    # Right: surface gravities vs mu.
    mus = np.linspace(mu_ext * 1.0008, 4.0, 600)
    kp = np.empty_like(mus)
    km = np.empty_like(mus)
    for i, mu in enumerate(mus):
        xm, xp = horizons(mu)
        kp[i] = kappa(xp)
        km[i] = kappa(xm)
    axR.plot(mus, kp, color=ks.NAVY, label=r"$\ell_0\kappa_+$")
    axR.plot(mus, km, color=ks.CORAL, label=r"$\ell_0\kappa_-$")
    axR.axvline(mu_ext, color=ks.GREY, ls="--", lw=0.8)
    axR.text(mu_ext + 0.04, 1.5, r"$\mu_{\rm ext}$", color=ks.GREY, fontsize=10,
             rotation=90, va="center", ha="left")
    axR.set_xlabel(r"$\mu=M/\ell_0$")
    axR.set_ylabel(r"$\ell_0\,\kappa_\pm$")
    axR.set_title(r"surface gravities")
    axR.set_xlim(1.2, 4)
    ks.grid(axR)
    ks.legend(axR, loc="upper left")

    fig.tight_layout()
    ks.savefig(fig, "horizon_phase_kappa")


if __name__ == "__main__":
    main()
