"""Shared plotting style and helpers for the experiment scripts."""

from __future__ import annotations

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

FIG_DIR = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG_DIR.mkdir(exist_ok=True)

# Reference categorical palette (fixed order) and a one-hue sequential ramp.
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
METHOD_COLOR = {"RK4": BLUE, "DP5 fixed": ORANGE, "DP5 adaptive": AQUA}
SEQ_BLUES = ["#9ec5f4", "#5598e7", "#256abf", "#0d366b"]
INK, INK_2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.edgecolor": INK_2,
    "axes.labelcolor": INK,
    "xtick.color": INK_2,
    "ytick.color": INK_2,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "lines.linewidth": 1.6,
    "lines.markersize": 4.5,
    "legend.frameon": False,
    "mathtext.fontset": "cm",
})


def save(fig, name: str) -> pathlib.Path:
    path = FIG_DIR / f"{name}.pdf"
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{name}.png", bbox_inches="tight")
    print(f"wrote {path}")
    return path


def slope_guide(ax, x0, y0, x1, slope, label):
    """Draw a dashed reference line y ~ x^slope through (x0, y0) on log axes."""
    y1 = y0 * (x1 / x0) ** slope
    ax.plot([x0, x1], [y0, y1], ls="--", lw=1.0, color=INK_2)
    ax.annotate(label, (x1, y1), textcoords="offset points", xytext=(4, 0),
                color=INK_2, fontsize=8, va="center")
