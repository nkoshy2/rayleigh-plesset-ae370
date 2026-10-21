"""Case B: convergence rate through collapse and rebound, and energy conservation.

(a) Fixed-step error at t = 3 t_R (past the first minimum and into the
    rebound) against a tight-tolerance reference; confirms orders 4 and 5 on
    the full physics (gas, viscosity, surface tension).
(b) Inviscid Case B conserves the energy of Eq. 8 exactly; the drift in E is
    an independent check that we solve the right equation.
"""

import numpy as np

from common import METHOD_COLOR, save, slope_guide
import matplotlib.pyplot as plt

from rp.integrators import dp5, dp5_fixed, rk4
from rp.model import case_B, energy, rhs, with_terms

T_END = 3.0


def main():
    prm = case_B(R0=100e-6, Pi=5.0)
    f = lambda t, y: rhs(t, y, prm)
    ref = dp5(f, (0, T_END), [1, 0], rtol=1e-13, atol=1e-15, positive=0)
    x_ref = ref.y[-1, 0]
    print(f"reference x(3) = {x_ref:.14f} with {ref.nfev} evaluations")

    rows = {"RK4": [], "DP5 fixed": []}
    for n in [50, 100, 200, 400, 800, 1600, 3200]:
        s = rk4(f, (0, T_END), [1, 0], n, store=False)
        rows["RK4"].append((T_END / n, s.nfev, abs(s.y[-1, 0] - x_ref)))
        s = dp5_fixed(f, (0, T_END), [1, 0], n // 2, store=False)
        rows["DP5 fixed"].append((2 * T_END / n, s.nfev, abs(s.y[-1, 0] - x_ref)))
    for k, r in rows.items():
        r = np.array(r)
        p = np.log(r[-2, 2] / r[-1, 2]) / np.log(r[-2, 0] / r[-1, 0])
        print(f"{k:10s} observed order (finest pair) = {p:.2f}")

    # Energy drift, inviscid Pi = 10 (stronger collapse), over several cycles
    prm_i = with_terms(case_B(R0=100e-6, Pi=10.0), viscosity=False)
    fi = lambda t, y: rhs(t, y, prm_i)
    T_E = 8.0
    E_rk4, E_dp = [], []
    for n in [400, 800, 1600, 3200, 6400, 12800]:
        s = rk4(fi, (0, T_E), [1, 0], n)
        E = energy(s.y, prm_i)
        E_rk4.append((s.nfev, np.max(np.abs(E - E[0])) / abs(E[0])))
    for tol in 10.0 ** -np.arange(5, 13):
        s = dp5(fi, (0, T_E), [1, 0], rtol=tol, atol=tol * 1e-2, positive=0)
        E = energy(s.y, prm_i)
        E_dp.append((s.nfev, np.max(np.abs(E - E[0])) / abs(E[0])))
    print("energy drift RK4:", [f"{e:.1e}" for _, e in E_rk4])
    print("energy drift DP5:", [f"{e:.1e}" for _, e in E_dp])

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8))
    ax = axes[0]
    for k, r in rows.items():
        r = np.array(r)
        ax.loglog(r[:, 0], r[:, 2], "o-", color=METHOD_COLOR[k], label=k)
    r4, r5 = np.array(rows["RK4"]), np.array(rows["DP5 fixed"])
    slope_guide(ax, r4[-1, 0], r4[-1, 2] / 8, r4[-4, 0], 4, "$h^4$")
    slope_guide(ax, r5[-1, 0], r5[-1, 2] / 8, r5[-4, 0], 5, "$h^5$")
    ax.set_xlabel(r"step size $h$ (units of $t_R$)")
    ax.set_ylabel(r"$|\tilde R(3) - \tilde R_{\mathrm{ref}}(3)|$")
    ax.set_title(r"(a) Case B, $\Pi=5$: convergence rate")
    ax.legend(loc="upper left")

    ax = axes[1]
    a = np.array(E_rk4)
    ax.loglog(a[:, 0], a[:, 1], "o-", color=METHOD_COLOR["RK4"], label="RK4")
    a = np.array(E_dp)
    ax.loglog(a[:, 0], a[:, 1], "o-", color=METHOD_COLOR["DP5 adaptive"],
              label="DP5 adaptive")
    ax.set_xlabel(r"evaluations of $\mathbf{f}$")
    ax.set_ylabel(r"$\max_t |E - E_0| / |E_0|$")
    ax.set_title(r"(b) Inviscid $\Pi=10$: energy drift, $0\leq\tilde t\leq 8$")
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, "fig_convergence_caseB")


if __name__ == "__main__":
    main()
