"""Case A: convergence to the exact Rayleigh solution (rubric 3(i), correctness).

Integrates the empty-cavity collapse to t* = t_exact(R = 0.5 R0) and compares
x(t*) with the exact value 0.5. Every method should converge to the exact
answer at its design rate. (Closer to R = 0 the asymptotic range starts at
smaller h because the derivatives of f grow; see the method section.) Also
checks the collapse time deep into the collapse, at R = 1e-3 R0.
"""

import numpy as np

from common import METHOD_COLOR, save, slope_guide
import matplotlib.pyplot as plt

from rp.integrators import Event, dp5, dp5_fixed, rk4
from rp.model import RAYLEIGH_CONSTANT, case_A, rayleigh_time, rhs

X_STAR = 0.5


def main():
    prm = case_A()
    f = lambda t, y: rhs(t, y, prm)
    t_star = rayleigh_time(X_STAR)
    print(f"t* = {t_star:.12f} (exact time to reach R = {X_STAR} R0)")

    rows = {"RK4": [], "DP5 fixed": [], "DP5 adaptive": []}
    for n in [16, 32, 64, 128, 256, 512, 1024]:
        s = rk4(f, (0, t_star), [1, 0], n, store=False)
        rows["RK4"].append((t_star / n, s.nfev, abs(s.y[-1, 0] - X_STAR)))
    for n in [8, 16, 32, 64, 128, 256, 512]:
        s = dp5_fixed(f, (0, t_star), [1, 0], n, store=False)
        rows["DP5 fixed"].append((t_star / n, s.nfev, abs(s.y[-1, 0] - X_STAR)))
    for tol in 10.0 ** -np.arange(4, 13):
        s = dp5(f, (0, t_star), [1, 0], rtol=tol, atol=tol * 1e-2, positive=0)
        rows["DP5 adaptive"].append((tol, s.nfev, abs(s.y[-1, 0] - X_STAR)))

    for k in ["RK4", "DP5 fixed"]:
        r = np.array(rows[k])
        p = np.log(r[-2, 2] / r[-1, 2]) / np.log(r[-2, 0] / r[-1, 0])
        print(f"{k:13s} observed order (finest pair) = {p:.2f}")

    ev = Event(lambda t, y: y[0] - 1e-3, direction=-1, terminal=True)
    s = dp5(f, (0, 2), [1, 0], rtol=1e-12, atol=1e-14, events=[ev], positive=0)
    t_ex = rayleigh_time(1e-3)
    print(f"t(x=1e-3): numerical {s.t[-1]:.12f}, exact {t_ex:.12f}, "
          f"diff {abs(s.t[-1] - t_ex):.1e}, t_c = {RAYLEIGH_CONSTANT:.6f}, "
          f"nfev = {s.nfev}")

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8))
    ax = axes[0]
    for k in ["RK4", "DP5 fixed"]:
        r = np.array(rows[k])
        ax.loglog(r[:, 0], r[:, 2], "o-", color=METHOD_COLOR[k], label=k)
    r4, r5 = np.array(rows["RK4"]), np.array(rows["DP5 fixed"])
    slope_guide(ax, r4[-1, 0], r4[-1, 2] / 8, r4[-4, 0], 4, "$h^4$")
    slope_guide(ax, r5[-1, 0], r5[-1, 2] / 8, r5[-4, 0], 5, "$h^5$")
    ax.set_xlabel(r"step size $h$ (units of $t_R$)")
    ax.set_ylabel(r"$|\tilde R(t^*) - 0.5|$")
    ax.set_title("(a) Error vs exact solution, fixed step")
    ax.legend(loc="upper left")

    ax = axes[1]
    for k in rows:
        r = np.array(rows[k])
        ax.loglog(r[:, 1], r[:, 2], "o-", color=METHOD_COLOR[k], label=k)
    ax.set_xlabel(r"evaluations of $\mathbf{f}$")
    ax.set_ylabel(r"$|\tilde R(t^*) - 0.5|$")
    ax.set_title("(b) Error vs cost")
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, "fig_verification_caseA")


if __name__ == "__main__":
    main()
