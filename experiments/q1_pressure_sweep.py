"""Q1: how hard does a bubble collapse after a pressure jump? (rubric 4)

Case B first collapse for pressure ratios Pi = 2, 5, 10, 20 and initial radii
10, 100 and 1000 um: R_min / R0, peak inward wall speed, peak gas pressure and
wall Mach number. Prints the numbers in report Table 2 and Section 4.1.
"""

import matplotlib.pyplot as plt

from common import AQUA, BLUE, ORANGE, save

from rp.diagnostics import first_collapse
from rp.model import case_B

PIS = [2.0, 5.0, 10.0, 20.0]
R0S = [10e-6, 100e-6, 1000e-6]
COLORS = [BLUE, ORANGE, AQUA]


def main():
    results = {}
    for R0 in R0S:
        rows = []
        print(f"\nR0 = {R0 * 1e6:g} um")
        print(f"  {'Pi':>4} {'R_min/R0':>9} {'|Rdot|max m/s':>14} "
              f"{'p_g,max MPa':>12} {'Mach':>7}")
        for Pi in PIS:
            prm = case_B(R0=R0, Pi=Pi)
            c = first_collapse(prm)
            speed = abs(c.v_min) * prm.u_R
            rows.append((Pi, c.x_min, c.p_g_max / 1e6))
            print(f"  {Pi:>4g} {c.x_min:>9.4f} {speed:>14.2f} "
                  f"{c.p_g_max / 1e6:>12.3f} {c.mach:>7.4f}")
        results[R0] = rows

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.7))
    for color, R0 in zip(COLORS, R0S):
        Pi, x_min, p_max = zip(*results[R0])
        label = rf"$R_0={R0 * 1e6:g}\,\mu$m"
        ax1.plot(Pi, x_min, "o-", color=color, label=label)
        ax2.semilogy(Pi, p_max, "o-", color=color, label=label)
    ax1.set_xlabel(r"pressure ratio $\Pi$")
    ax1.set_ylabel(r"$R_{\min}/R_0$")
    ax1.set_title(r"(a) first-collapse minimum radius")
    ax2.set_xlabel(r"pressure ratio $\Pi$")
    ax2.set_ylabel("peak gas pressure (MPa)")
    ax2.set_title(r"(b) peak gas pressure")
    ax1.legend()
    fig.tight_layout()
    save(fig, "q1_pressure_sweep")


if __name__ == "__main__":
    main()
