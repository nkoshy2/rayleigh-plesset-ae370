"""Case C: small-amplitude forced response vs linear theory (rubric 3(i)).

At eps = 1e-4 the response is linear to within O(eps), so the steady amplitude
from the full nonlinear solver must match the damped driven oscillator of
Eq. 10 across the resonance.
"""

import math

import numpy as np

from common import AQUA, BLUE, save
import matplotlib.pyplot as plt

from rp.integrators import dp5
from rp.model import case_C, linear_amplitude, rhs

EPS = 1e-4
N_SETTLE = 700.0   # dimensionless time for the start-up transient to decay
N_MEASURE = 40     # periods over which the amplitude is measured


def steady_amplitude(prm):
    f = lambda t, y: rhs(t, y, prm)
    Om = prm.extra["Omega"]
    period = 2 * math.pi / Om
    T = (math.ceil(N_SETTLE / period) + N_MEASURE) * period
    sol = dp5(f, (0.0, T), [1.0, 0.0], rtol=1e-11, atol=1e-14,
              h_max=0.05 * period, positive=0)
    tail = sol.t > T - N_MEASURE * period
    x = sol.y[tail, 0] - 1.0
    return 0.5 * (x.max() - x.min()), sol.nfev


def main():
    ratios = np.array([0.5, 0.7, 0.85, 0.95, 1.0, 1.05, 1.15, 1.3, 1.6, 2.0])
    num, lin = [], []
    for r in ratios:
        prm = case_C(R0=10e-6, eps=EPS, omega_ratio=r)
        a, nfev = steady_amplitude(prm)
        num.append(a)
        lin.append(linear_amplitude(prm))
        print(f"w/w0={r:4.2f}: numerical {a:.6e}, linear {lin[-1]:.6e}, "
              f"rel diff {abs(a / lin[-1] - 1):.1e}, nfev {nfev}")
    prm = case_C(R0=10e-6, eps=EPS)
    w0 = prm.extra["w0"]
    print(f"w0 t_R = {w0:.4f}, f0 = {w0 / prm.t_R / (2 * math.pi) / 1e3:.1f} kHz, "
          f"Q = {w0 / (4 * prm.inv_Re):.1f}")

    dense = np.linspace(0.4, 2.1, 400)
    lin_dense = [linear_amplitude(case_C(R0=10e-6, eps=EPS, omega_ratio=r))
                 for r in dense]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
    ax = axes[0]
    ax.semilogy(dense, np.array(lin_dense) / EPS, color=BLUE, label="linear theory")
    ax.semilogy(ratios, np.array(num) / EPS, "o", color=AQUA,
                label="DP5, full equation")
    ax.set_xlabel(r"$\omega/\omega_0$")
    ax.set_ylabel(r"amplitude $|\tilde R - 1|/\varepsilon$")
    ax.set_title(r"(a) $R_0=10\ \mu$m, $\varepsilon=10^{-4}$")
    ax.legend(loc="upper right", fontsize=7)
    ax = axes[1]
    ax.semilogy(ratios, np.abs(np.array(num) / np.array(lin) - 1), "o-", color=AQUA)
    ax.set_xlabel(r"$\omega/\omega_0$")
    ax.set_ylabel("relative difference")
    ax.set_title("(b) Numerical vs linear amplitude")
    fig.tight_layout()
    save(fig, "fig_verification_caseC")


if __name__ == "__main__":
    main()
