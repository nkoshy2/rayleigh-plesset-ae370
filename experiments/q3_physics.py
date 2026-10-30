"""Q3: which physics resists the collapse, and when? (rubric 4)

Case B at Pi = 5 for R0 = 1, 10 and 100 um, comparing the full model with runs
where viscosity (mu = 0) or surface tension (sigma = 0) is removed. Removing a
term is done through case_B, so p_g0 is recomputed for the matching initial
equilibrium. Prints the numbers quoted in report Section 4.3.
"""

import matplotlib.pyplot as plt

from common import AQUA, BLUE, ORANGE, save

from rp.diagnostics import first_collapse
from rp.model import case_B

PI = 5.0
R0S = [1e-6, 10e-6, 100e-6]
VARIANTS = [("Full model", {}, BLUE),
            ("No viscosity", {"mu": 0.0}, ORANGE),
            ("No surface tension", {"sigma": 0.0}, AQUA)]


def main():
    results = {}
    for name, kw, _ in VARIANTS:
        results[name] = [first_collapse(case_B(R0=R0, Pi=PI, **kw), rtol=1e-8,
                                         atol=1e-10).x_min
                         for R0 in R0S]
    print(f"R_min / R0 at Pi = {PI:g}")
    for i, R0 in enumerate(R0S):
        vals = ", ".join(f"{name} {results[name][i]:.4f}" for name, _, _ in VARIANTS)
        print(f"  R0 = {R0 * 1e6:>5g} um: {vals}")

    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    for name, _, color in VARIANTS:
        ax.semilogx([R0 * 1e6 for R0 in R0S], results[name], "o-",
                    color=color, label=name)
    ax.set_xlabel(r"initial radius $R_0$ ($\mu$m)")
    ax.set_ylabel(r"$R_{\min}/R_0$")
    ax.legend()
    fig.tight_layout()
    save(fig, "q3_physics")


if __name__ == "__main__":
    main()
