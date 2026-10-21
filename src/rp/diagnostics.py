"""Collapse diagnostics built on the adaptive integrator (Q1 quantities)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .integrators import Event, dp5
from .model import C_SOUND, Params, rhs


@dataclass
class Collapse:
    x_min: float        # R_min / R0
    t_min: float        # time of R_min, in units of t_R
    v_min: float        # most negative x' before R_min (peak inward wall speed)
    p_g_max: float      # peak gas pressure [Pa]
    mach: float         # |R'|_max / c
    nfev: int
    status: str


def first_collapse(prm: Params, rtol: float = 1e-10, atol: float = 1e-12,
                   t_end: float = 10.0, **kw) -> Collapse:
    """Integrate to the first minimum of R and report the collapse quantities.

    The minimum of R is the upward zero crossing of x'; the peak inward speed
    is the upward zero crossing of x'' (the minimum of x'). Both are located
    by the integrator's fifth-order event search.
    """
    f = lambda t, y: rhs(t, y, prm)
    ev_min = Event(lambda t, y: y[1], direction=+1, terminal=True)
    ev_vmin = Event(lambda t, y: f(t, y)[1], direction=+1, terminal=False)
    sol = dp5(f, (0.0, t_end), [1.0, 0.0], rtol=rtol, atol=atol,
              events=[ev_min, ev_vmin], store=False, positive=0, **kw)
    if sol.status != "event":
        raise RuntimeError(f"no minimum found before t = {t_end}: {sol.message}")
    t_min = next(te for j, te in sol.t_events if j == 0)
    x_min = float(sol.y[-1, 0])
    v_events = [ye[1] for j, ye in sol.y_events if j == 1]
    v_min = float(min(v_events)) if v_events else float("nan")
    p_g0 = prm.extra.get("p_g0", prm.Pg * prm.dp)
    p_g_max = p_g0 * x_min ** (-3.0 * prm.kappa)
    mach = abs(v_min) * prm.u_R / C_SOUND
    return Collapse(x_min, t_min, v_min, p_g_max, mach, sol.nfev, sol.status)
