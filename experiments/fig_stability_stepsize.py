"""Stability check along real trajectories, and adaptive step-size history.

fig_stability:
  (a) Boundaries |R(z)| = 1 of RK4 and DP5 with the scaled eigenvalues
      h*lambda of the Jacobian at every accepted DP5 step.
  (b) max |h*lambda| over time against each method's stability limits. If the
      curve stays far below the limits, the step was set by accuracy, not
      stability (method Section 2.2).
fig_stepsize: R(t) and the accepted step size h(t) for a violent collapse.
"""

import numpy as np

from common import AQUA, BLUE, INK_2, ORANGE, save
import matplotlib.pyplot as plt

from rp.integrators import dp5, rk4, stability_polynomial
from rp.model import case_B, case_C, jacobian, rhs

RTOL = 1e-8


def scaled_eigs(prm, T):
    f = lambda t, y: rhs(t, y, prm)
    sol = dp5(f, (0, T), [1, 0], rtol=RTOL, atol=RTOL * 1e-2, positive=0)
    h = np.diff(sol.t)[:-1]  # drop the final step clipped to land on T
    lam = np.array([np.linalg.eigvals(jacobian(sol.t[k], sol.y[k], prm))
                    for k in range(len(h))])
    return sol, h, h[:, None] * lam


def imag_limit(R):
    y = 1e-3
    while abs(R(1j * y)) <= 1.0:
        y += 1e-3
    return y


def real_limit(R):
    x = 1e-3
    while abs(R(-x)) <= 1.0:
        x += 1e-3
    return x


def main():
    prm_b = case_B(R0=100e-6, Pi=20.0)
    sol_b, h_b, z_b = scaled_eigs(prm_b, 3.0)
    prm_c = case_C(R0=10e-6, eps=0.3, omega_ratio=1.0)
    Om = prm_c.extra["Omega"]
    sol_c, h_c, z_c = scaled_eigs(prm_c, 20 * 2 * np.pi / Om)
    for name, z in [("Case B Pi=20", z_b), ("Case C eps=0.3", z_c)]:
        print(f"{name}: max |h lambda| = {np.max(np.abs(z)):.3f}, "
              f"min Re = {z.real.min():.3f}, max Re = {z.real.max():.3f}")
    limits = {}
    for name in ["rk4", "dp5"]:
        R = stability_polynomial(name)
        limits[name] = (real_limit(R), imag_limit(R))
        print(f"{name}: real-axis limit {limits[name][0]:.3f}, "
              f"imaginary-axis limit {limits[name][1]:.3f}")

    # How many fixed RK4 steps match the adaptive accuracy for Pi = 20 to t = 3?
    f = lambda t, y: rhs(t, y, prm_b)
    ref = dp5(f, (0, 3.0), [1, 0], rtol=1e-13, atol=1e-15, positive=0)
    err_dp = abs(sol_b.y[-1, 0] - ref.y[-1, 0])
    n = 100
    while True:
        # coarse RK4 runs overshoot to R < 0 and return NaN; they just fail the test
        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            s = rk4(f, (0, 3.0), [1, 0], n, store=False)
        if abs(s.y[-1, 0] - ref.y[-1, 0]) <= err_dp:
            break
        n = int(n * 1.25)
    print(f"Pi=20, rtol={RTOL:g}: DP5 error {err_dp:.1e} with {sol_b.nfev} f-evals "
          f"({sol_b.n_accept} accepted, {sol_b.n_reject} rejected); RK4 needs "
          f"~{n} steps = {4 * n} f-evals for the same error")
    print(f"step size range: {h_b.min():.2e} to {h_b.max():.2e} "
          f"(ratio {h_b.max() / h_b.min():.0f})")

    # ---------------------------------------------------------- stability --
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0),
                             gridspec_kw=dict(width_ratios=[1.0, 1.4]))
    ax = axes[0]
    xs = np.linspace(-4, 0.6, 500)
    ys = np.linspace(-3.6, 3.6, 700)
    Z = xs[None, :] + 1j * ys[:, None]
    for name, color in [("rk4", BLUE), ("dp5", ORANGE)]:
        R = np.vectorize(stability_polynomial(name))(Z)
        ax.contour(xs, ys, np.abs(R), levels=[1.0], colors=color, linewidths=1.6)
    ax.plot([], [], color=BLUE, label="RK4")
    ax.plot([], [], color=ORANGE, label="DP5")
    ax.scatter(z_b.real.ravel(), z_b.imag.ravel(), s=6, color=AQUA,
               linewidths=0, zorder=3, label=r"$h\lambda$ (B and C)")
    ax.scatter(z_c.real.ravel(), z_c.imag.ravel(), s=6, color=AQUA,
               linewidths=0, zorder=3)
    ax.axhline(0, color=INK_2, lw=0.6)
    ax.axvline(0, color=INK_2, lw=0.6)
    ax.set_xlim(-4, 0.6)
    ax.set_ylim(-3.6, 3.6)
    ax.set_aspect("equal")
    ax.set_xlabel(r"Re$(h\lambda)$")
    ax.set_ylabel(r"Im$(h\lambda)$")
    ax.set_title(r"(a) $|R(z)|=1$ and $h\lambda$")
    ax.legend(loc="lower left", fontsize=7)

    ax = axes[1]
    ax.semilogy(sol_b.t[:len(h_b)], np.abs(z_b).max(axis=1), color=AQUA,
                label=r"Case B, $\Pi=20$")
    tc = sol_c.t[:len(h_c)] * (3.0 / sol_c.t[len(h_c)])  # rescale to share axis
    ax.semilogy(tc, np.abs(z_c).max(axis=1), color=INK_2, lw=1.0,
                label=r"Case C, $\varepsilon=0.3$ (20 periods, rescaled)")
    for val, color, text in [(limits["rk4"][0], BLUE, "RK4 real limit 2.79"),
                             (limits["dp5"][1], ORANGE, "DP5 imag. limit 1.00")]:
        ax.axhline(val, color=color, ls="--", lw=1.0)
        ax.annotate(text, (0.02, val), xycoords=("axes fraction", "data"),
                    textcoords="offset points", xytext=(0, 3), color=INK_2,
                    fontsize=7)
    ax.set_ylim(1e-3, 6)
    ax.set_xlabel(r"$t/t_R$ (Case B)")
    ax.set_ylabel(r"$\max |h\lambda|$")
    ax.set_title(r"(b) Largest $|h\lambda|$ per accepted DP5 step")
    ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    save(fig, "fig_stability")

    # ---------------------------------------------------------- step size --
    fig, axes = plt.subplots(2, 1, figsize=(4.6, 3.2), sharex=True)
    axes[0].plot(sol_b.t, sol_b.y[:, 0], color=BLUE)
    axes[0].set_ylabel(r"$R/R_0$")
    axes[0].set_title(r"Case B, $\Pi=20$, DP5 rtol $=10^{-8}$")
    axes[1].semilogy(sol_b.t[:len(h_b)], h_b, ".", color=AQUA, ms=3)
    axes[1].set_ylabel(r"step $h$ ($t_R$)")
    axes[1].set_xlabel(r"$t/t_R$")
    fig.tight_layout()
    save(fig, "fig_stepsize")


if __name__ == "__main__":
    main()
