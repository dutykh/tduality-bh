#!/usr/bin/env python3
# Authors: Dr. Davide Batic and Dr. Denys Dutykh
#          (Mathematics Department, Khalifa University of Science and
#           Technology, Abu Dhabi, UAE)
"""Light-ring and horizon thresholds of the T-duality-inspired Bardeen
black hole in the (x, mu) plane, x = r/l0, mu = M/l0.

Horizon condition f(x)=0:        mu = H(x) = (x^2 + 1)^{3/2} / (2 x^2),
Circular-null condition x f'-2f:  mu = P(x) = (x^2 + 1)^{5/2} / (3 x^4).

H has its minimum mu_ext = 3 sqrt(3)/4 at x = sqrt(2); P has its minimum
mu_LR = 25 sqrt(5)/48 at x = 2. Since mu_LR < mu_ext, a horizonless yet
ultracompact band mu_LR < mu < mu_ext carries a stable (inner) and an
unstable (outer) light ring.
"""

import numpy as np

import ku_style as ks
from ku_style import plt


def H(x):
    return (x**2 + 1.0) ** 1.5 / (2.0 * x**2)


def P(x):
    return (x**2 + 1.0) ** 2.5 / (3.0 * x**4)


def main():
    mu_ext = 3.0 * np.sqrt(3.0) / 4.0
    mu_LR = 25.0 * np.sqrt(5.0) / 48.0
    assert abs(H(np.sqrt(2.0)) - mu_ext) < 1e-12
    assert abs(P(2.0) - mu_LR) < 1e-12
    assert mu_LR < mu_ext

    x = np.linspace(0.55, 6.0, 2000)
    fig, ax = plt.subplots(figsize=(7.4, 5.2))

    ax.axhspan(mu_LR, mu_ext, color=ks.SKY, alpha=0.25, lw=0)
    ax.plot(x, H(x), color=ks.NAVY, label=r"horizon $\mu=H(x)$")
    ax.plot(x, P(x), color=ks.CORAL, label=r"light ring $\mu=P(x)$")

    ax.axhline(mu_ext, color=ks.NAVY, ls="--", lw=0.8, alpha=0.7)
    ax.axhline(mu_LR, color=ks.CORAL, ls="--", lw=0.8, alpha=0.7)
    ax.plot([np.sqrt(2.0)], [mu_ext], "o", color=ks.NAVY, ms=5)
    ax.plot([2.0], [mu_LR], "o", color=ks.CORAL, ms=5)

    ax.text(5.0, mu_ext + 0.02, r"$\mu_{\rm ext}=\frac{3\sqrt{3}}{4}$",
            color=ks.NAVY, fontsize=10, va="bottom")
    ax.text(5.0, mu_LR - 0.02, r"$\mu_{\rm LR}=\frac{25\sqrt{5}}{48}$",
            color=ks.CORAL, fontsize=10, va="top")
    ax.text(4.3, 0.5 * (mu_LR + mu_ext), "ultracompact branch", fontsize=9,
            color=ks.GREY, ha="center", va="center")

    ax.set_xlabel(r"$x=r/\ell_0$")
    ax.set_ylabel(r"$\mu=M/\ell_0$")
    ax.set_xlim(0.5, 6)
    ax.set_ylim(1.02, 3.2)
    ks.grid(ax)
    ks.legend(ax, loc="upper center")
    fig.tight_layout()
    ks.savefig(fig, "lightring_structure")


if __name__ == "__main__":
    main()
