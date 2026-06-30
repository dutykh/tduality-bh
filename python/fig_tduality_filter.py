#!/usr/bin/env python3
# Authors: Dr. Davide Batic and Dr. Denys Dutykh
#          (Mathematics Department, Khalifa University of Science and
#           Technology, Abu Dhabi, UAE)
"""Zero-point-length form factor and the classical-versus-filtered
Cauchy-horizon mass inflation.

Left: the path-integral-duality form factor Z(z) = z K_1(z) and its squared
flux factor W(z) = z^2 K_1(z)^2, with W(z) -> 1 (z << 1) and
W(z) ~ (pi z / 2) e^{-2 z} (z >> 1).

Right: with dimensionless advanced time V = v/l0, surface gravity l0 kappa = 1
and Price exponent p = 2, the classical mass-inflation factor V^{-p} e^{kappa V}
versus the filtered factor V^{-p} e^{-beta kappa V} for beta = 0, 1, 2.
"""

import numpy as np
from scipy.special import kv

import ku_style as ks
from ku_style import plt


def Z(z):
    return z * kv(1, z)


def W(z):
    return (z * kv(1, z)) ** 2


def main():
    # --- Left panel: form factor ---
    z = np.logspace(-1.3, 0.95, 800)
    figL_ok = np.allclose(W(1e-3), 1.0, atol=1e-3)
    assert figL_ok, "W(z)->1 as z->0 failed"
    # UV asymptote check
    zb = 6.0
    assert abs(W(zb) - 0.5 * np.pi * zb * np.exp(-2 * zb)) / W(zb) < 0.2

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10.4, 4.5))

    axL.semilogy(z, Z(z), color=ks.NAVY, label=r"$Z(z)=zK_1(z)$")
    axL.semilogy(z, W(z), color=ks.CORAL, label=r"$W(z)=z^2K_1(z)^2$")
    zU = np.linspace(2.0, 9.0, 200)
    axL.semilogy(zU, 0.5 * np.pi * zU * np.exp(-2 * zU), color=ks.GREY, ls=":",
                 lw=1.2, label=r"$\frac{\pi z}{2}e^{-2z}$")
    axL.axhline(1.0, color=ks.GREY, lw=0.6, alpha=0.5)
    axL.set_xlabel(r"$z=\ell_0\,\Omega$")
    axL.set_ylabel(r"form factor")
    axL.set_xlim(0, 9)
    axL.set_ylim(1e-6, 2.0)
    axL.set_title(r"zero-point-length filter")
    ks.grid(axL)
    ks.legend(axL, loc="lower left")

    # --- Right panel: classical vs filtered mass inflation ---
    V = np.linspace(0.6, 12.0, 1000)
    p, k = 2.0, 1.0
    classical = V ** (-p) * np.exp(k * V)
    axR.semilogy(V, classical, color=ks.CORAL,
                 label=r"classical $V^{-p}e^{\kappa V}$")
    for beta, c in [(0, ks.NAVY), (1, ks.TEAL), (2, ks.SKY)]:
        axR.semilogy(V, V ** (-p) * np.exp(-beta * k * V), color=c,
                     label=rf"filtered $\beta={beta}$")
    v_TD = np.log(1.0 / 0.1) / k  # onset for l0*omega = 0.1
    axR.axvline(v_TD, color=ks.GREY, ls="--", lw=0.8)
    axR.text(v_TD + 0.1, 1e3, r"$v_{\rm TD}$", color=ks.GREY, fontsize=10)
    axR.set_xlabel(r"$V=v/\ell_0$")
    axR.set_ylabel(r"$\delta\mathcal{M}(V)$")
    axR.set_xlim(0, 12)
    axR.set_ylim(1e-6, 1e5)
    axR.set_title(r"mass inflation ($\ell_0\kappa=1,\ p=2$)")
    ks.grid(axR)
    ks.legend(axR, loc="upper center", ncol=1)

    fig.tight_layout()
    ks.savefig(fig, "tduality_filter")


if __name__ == "__main__":
    main()
