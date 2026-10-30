"""Q2: response to oscillating pressure, and where it stops being linear (rubric 4)

Case C at R0 = 10 um, forcing amplitudes eps = 0.01, 0.1, 0.3 and frequencies
w / w0 = 0.7 to 1.3. Each run lasts 200 forcing periods; the response is half
the peak-to-peak variation of R / R0 over the last 20 periods, by which time
the start-up transient (quality factor about 54) has decayed. Prints the
numbers quoted in report Section 4.2.
"""

import math

import numpy as np
import matplotlib.pyplot as plt

from common import AQUA, BLUE, ORANGE, save

from rp.integrators import dp5
from rp.model import case_C, rhs

EPSS = [0.01, 0.1, 0.3]
RATIOS = np.round(np.arange(0.7, 1.31, 0.1), 2)
COLORS = [BLUE, ORANGE, AQUA]
N_PERIODS, N_MEASURE = 200, 20


def steady_response(eps, ratio):
    """Return (half peak-to-peak, max, min) of R / R0 over the last periods."""
    prm = case_C(R0=10e-6, eps=eps, omega_ratio=ratio)
    T = 2.0 * math.pi / prm.extra["Omega"]
    # h_max keeps at least 50 samples per period so the max/min are resolved.
    sol = dp5(lambda t, y: rhs(t, y, prm), (0.0, N_PERIODS * T), [1.0, 0.0],
              rtol=1e-8, atol=1e-10, positive=0, h_max=T / 50)
    x = sol.y[sol.t >= (N_PERIODS - N_MEASURE) * T, 0]
    return 0.5 * (x.max() - x.min()), x.max(), x.min()


def main():
    results = {}
    for eps in EPSS:
        print(f"\neps = {eps:g}")
        rows = []
        for ratio in RATIOS:
            amp, x_max, x_min = steady_response(eps, ratio)
            rows.append(amp)
            print(f"  w/w0 = {ratio:.1f}: response {amp:.4f}, "
                  f"R_max/R0 = {x_max:.3f}, R_min/R0 = {x_min:.3f}")
        results[eps] = rows

    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    for color, eps in zip(COLORS, EPSS):
        ax.plot(RATIOS, results[eps], "o-", color=color,
                label=rf"$\varepsilon={eps:g}$")
    ax.set_xlabel(r"forcing frequency $\omega/\omega_0$")
    ax.set_ylabel(r"$(R_{\max}-R_{\min})/(2R_0)$")
    ax.legend()
    fig.tight_layout()
    save(fig, "q2_resonance")


if __name__ == "__main__":
    main()
