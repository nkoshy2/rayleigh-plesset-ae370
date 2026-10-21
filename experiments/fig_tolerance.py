"""Which tolerance do we need? (rubric 3(ii))

For Case B across the pressure ratios we will sweep, measure the relative error
of the Q1 quantities (R_min, peak inward wall speed, peak gas pressure)
against a reference at rtol = 1e-13, as a function of rtol. Also records cost.
"""

import numpy as np

from common import INK_2, SEQ_BLUES, save
import matplotlib.pyplot as plt

from rp.diagnostics import first_collapse
from rp.model import case_B

PIS = [2.0, 5.0, 10.0, 20.0]
RTOLS = 10.0 ** -np.arange(4, 12)
TARGET = 1e-6


def main():
    results = {}
    for Pi in PIS:
        prm = case_B(R0=100e-6, Pi=Pi)
        ref = first_collapse(prm, rtol=1e-13, atol=1e-15)
        rows = []
        for tol in RTOLS:
            c = first_collapse(prm, rtol=tol, atol=tol * 1e-2)
            rows.append((tol, c.nfev,
                         abs(c.x_min / ref.x_min - 1),
                         abs(c.v_min / ref.v_min - 1),
                         abs(c.p_g_max / ref.p_g_max - 1)))
        results[Pi] = np.array(rows)
        print(f"Pi={Pi:>4g}: R_min/R0={ref.x_min:.6f}, |v|max={abs(ref.v_min):.4f}, "
              f"p_g,max={ref.p_g_max:.4e} Pa, Mach={ref.mach:.3f}")

    print(f"\nLargest rtol giving all three relative errors < {TARGET:g}:")
    for Pi, r in results.items():
        ok = r[np.max(r[:, 2:], axis=1) < TARGET]
        tol = ok[0, 0] if len(ok) else float("nan")
        print(f"  Pi={Pi:>4g}: rtol = {tol:.0e}, nfev = {int(ok[0, 1])}")
    print("\nAt rtol = 1e-8:")
    for Pi, r in results.items():
        k = int(np.argmin(np.abs(np.log10(r[:, 0]) + 8)))
        print(f"  Pi={Pi:>4g}: errors x {r[k, 2]:.1e}, v {r[k, 3]:.1e}, "
              f"p {r[k, 4]:.1e}; nfev {int(r[k, 1])}")

    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.6), sharey=True)
    titles = [r"(a) $R_{\min}$", r"(b) peak inward $|\dot R|$",
              r"(c) peak gas pressure"]
    for j, ax in enumerate(axes):
        for color, Pi in zip(SEQ_BLUES, PIS):
            r = results[Pi]
            ax.loglog(r[:, 0], np.maximum(r[:, 2 + j], 1e-16), "o-", color=color,
                      ms=3.5, label=rf"$\Pi={Pi:g}$")
        ax.axhline(TARGET, color=INK_2, ls="--", lw=1.0)
        ax.set_xlabel("rtol")
        ax.set_title(titles[j])
        ax.invert_xaxis()
    axes[0].set_ylabel("relative error")
    axes[0].annotate(r"target $10^{-6}$", (RTOLS[-1], TARGET),
                     textcoords="offset points", xytext=(-2, 3), fontsize=7,
                     color=INK_2, ha="right")
    axes[2].legend(loc="lower left", fontsize=7)
    fig.tight_layout()
    save(fig, "fig_tolerance")


if __name__ == "__main__":
    main()
