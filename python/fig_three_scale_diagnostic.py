#!/usr/bin/env python3
# Authors: Dr. Davide Batic and Dr. Denys Dutykh
#          (Mathematics Department, Khalifa University of Science and
#           Technology, Abu Dhabi, UAE)
"""Three local radial diagnostics of the T-duality-inspired Bardeen core.

In dimensionless radius x = r/l0, with l0 the zero-point length, the
Newtonian radial acceleration, the Ricci focusing term, and the Raychaudhuri
compression rate read

    a~(x)      = -x (x^2 - 2) / (x^2 + 1)^{5/2}                 (= l0^2/M a_r)
    R~(x)      =  3 (3 x^2 - 2) / (x^2 + 1)^{7/2}               (= l0^3/M R_uu)
    Theta_dot~ = -3 x^2 (3 x^2 + 10) / (2 (x^2 + 1)^{7/2})      (= l0^3/M dTheta/dtau)

The figure shows that the Newtonian force-reversal radius x_F = sqrt(2), the
strong-energy-condition radius x_SEC = sqrt(2/3), and the Raychaudhuri
compression extremum x_Theta = sqrt((-19 + sqrt(541))/9) are three distinct
scales.
"""

import numpy as np

import ku_style as ks
from ku_style import plt


def a_tilde(x):
    return -x * (x**2 - 2.0) / (x**2 + 1.0) ** 2.5


def R_tilde(x):
    return 3.0 * (3.0 * x**2 - 2.0) / (x**2 + 1.0) ** 3.5


def thetadot_tilde(x):
    return -3.0 * x**2 * (3.0 * x**2 + 10.0) / (2.0 * (x**2 + 1.0) ** 3.5)


def main():
    x_F = np.sqrt(2.0)
    x_SEC = np.sqrt(2.0 / 3.0)
    x_T = np.sqrt((-19.0 + np.sqrt(541.0)) / 9.0)

    # Sanity checks against the analytic special radii.
    assert abs(a_tilde(x_F)) < 1e-12, "a~ must vanish at x_F=sqrt(2)"
    assert abs(R_tilde(x_SEC)) < 1e-12, "R~ must vanish at x_SEC=sqrt(2/3)"
    dx = 1e-6
    dprime = (thetadot_tilde(x_T + dx) - thetadot_tilde(x_T - dx)) / (2 * dx)
    assert abs(dprime) < 1e-5, "Theta_dot~ must be extremal at x_T"

    x = np.linspace(1e-3, 4.0, 1600)
    fig, ax = plt.subplots(figsize=(7.2, 5.0))

    ax.axhline(0.0, color=ks.GREY, lw=0.8, alpha=0.6)
    ax.plot(x, a_tilde(x), color=ks.NAVY, label=r"$\tilde a(x)=\ell_0^2 a_r/M$")
    ax.plot(x, R_tilde(x), color=ks.CORAL,
            label=r"$\tilde R(x)=\ell_0^3\,R_{\mu\nu}u^\mu u^\nu/M$")
    ax.plot(x, thetadot_tilde(x), color=ks.TEAL,
            label=r"$\tilde{\dot\Theta}(x)=\ell_0^3\,\dot\Theta/M$")

    for xv in (x_T, x_SEC, x_F):
        ax.axvline(xv, color=ks.GREY, ls="--", lw=0.8, alpha=0.7)

    ax.set_xlim(0, 4)
    ax.set_ylim(1.15 * thetadot_tilde(x_T), 2.6)
    ytop = 2.32
    ax.annotate(r"$x_{\dot\Theta}\simeq0.688$", (x_T, ytop), rotation=90,
                va="top", ha="right", fontsize=9, color=ks.GREY)
    ax.annotate(r"$x_{\rm SEC}=\sqrt{2/3}$", (x_SEC, ytop), rotation=90,
                va="top", ha="right", fontsize=9, color=ks.GREY)
    ax.annotate(r"$x_F=\sqrt{2}$", (x_F, ytop), rotation=90,
                va="top", ha="right", fontsize=9, color=ks.GREY)

    ax.set_xlabel(r"$x=r/\ell_0$")
    ax.set_ylabel(r"dimensionless local diagnostics")
    ks.grid(ax)
    ks.legend(ax, loc="lower right")
    fig.tight_layout()
    ks.savefig(fig, "three_scale_diagnostic")


if __name__ == "__main__":
    main()
