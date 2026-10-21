"""From-scratch explicit Runge-Kutta integrators (report Section 2).

* ``rk4``        classical fixed-step fourth-order method (Algorithm 1)
* ``dp5_fixed``  Dormand-Prince 5(4) with a fixed step, used only to show order 5
* ``dp5``        adaptive Dormand-Prince 5(4) with FSAL and events (Algorithm 2)

All take f(t, y) -> ndarray and count evaluations of f, the machine-independent
cost measure used in the report. No library time integrator is used anywhere.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np

RHS = Callable[[float, np.ndarray], np.ndarray]

# --------------------------------------------------------------------------- #
# Butcher tableaus
# --------------------------------------------------------------------------- #

RK4_C = np.array([0.0, 0.5, 0.5, 1.0])
RK4_A = np.array([[0.0, 0.0, 0.0, 0.0],
                  [0.5, 0.0, 0.0, 0.0],
                  [0.0, 0.5, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0]])
RK4_B = np.array([1 / 6, 1 / 3, 1 / 3, 1 / 6])

DP_C = np.array([0.0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0])
DP_A = np.array([
    [0, 0, 0, 0, 0, 0, 0],
    [1 / 5, 0, 0, 0, 0, 0, 0],
    [3 / 40, 9 / 40, 0, 0, 0, 0, 0],
    [44 / 45, -56 / 15, 32 / 9, 0, 0, 0, 0],
    [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729, 0, 0, 0],
    [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656, 0, 0],
    [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0],
], dtype=float)
DP_B = np.array([35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0])
DP_BSTAR = np.array([5179 / 57600, 0, 7571 / 16695, 393 / 640,
                     -92097 / 339200, 187 / 2100, 1 / 40])
DP_E = DP_B - DP_BSTAR  # error weights: e = h * sum(E_i k_i)


# --------------------------------------------------------------------------- #
# Single steps
# --------------------------------------------------------------------------- #

def rk4_step(f: RHS, t: float, y: np.ndarray, h: float) -> np.ndarray:
    """One classical RK4 step (Algorithm 1). Costs 4 evaluations of f."""
    k1 = f(t, y)
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1)
    k3 = f(t + 0.5 * h, y + 0.5 * h * k2)
    k4 = f(t + h, y + h * k3)
    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def dp5_step(f: RHS, t: float, y: np.ndarray, h: float, k1: np.ndarray):
    """One Dormand-Prince step given k1 = f(t, y).

    Returns (y_new, err, k7) where y_new is the fifth-order solution,
    err = y_new - y_new* is the embedded error estimate and k7 = f(t+h, y_new)
    is the next step's k1 (FSAL). Costs 6 evaluations of f.
    """
    K = np.empty((7, y.size))
    K[0] = k1
    for i in range(1, 7):
        K[i] = f(t + DP_C[i] * h, y + h * (DP_A[i, :i] @ K[:i]))
    y_new = y + h * (DP_B @ K)
    err = h * (DP_E @ K)
    return y_new, err, K[6]


# --------------------------------------------------------------------------- #
# Fixed-step drivers
# --------------------------------------------------------------------------- #

@dataclass
class Solution:
    t: np.ndarray
    y: np.ndarray
    nfev: int
    n_accept: int = 0
    n_reject: int = 0
    status: str = "success"
    message: str = ""
    t_events: list = field(default_factory=list)
    y_events: list = field(default_factory=list)

    @property
    def h(self) -> np.ndarray:
        return np.diff(self.t)


def rk4(f: RHS, t_span: Sequence[float], y0: Sequence[float], n_steps: int,
        store: bool = True) -> Solution:
    """Integrate with n_steps equal RK4 steps over t_span."""
    t0, tf = map(float, t_span)
    h = (tf - t0) / n_steps
    y = np.array(y0, dtype=float)
    ts, ys = [t0], [y.copy()]
    for k in range(n_steps):
        y = rk4_step(f, t0 + k * h, y, h)
        if store:
            ts.append(t0 + (k + 1) * h)
            ys.append(y.copy())
    if not store:
        ts.append(tf)
        ys.append(y.copy())
    return Solution(np.array(ts), np.array(ys), nfev=4 * n_steps,
                    n_accept=n_steps)


def dp5_fixed(f: RHS, t_span: Sequence[float], y0: Sequence[float],
              n_steps: int, store: bool = True) -> Solution:
    """Dormand-Prince fifth-order solution with a fixed step (no control)."""
    t0, tf = map(float, t_span)
    h = (tf - t0) / n_steps
    y = np.array(y0, dtype=float)
    k1 = f(t0, y)
    ts, ys = [t0], [y.copy()]
    for k in range(n_steps):
        y, _, k1 = dp5_step(f, t0 + k * h, y, h, k1)
        if store:
            ts.append(t0 + (k + 1) * h)
            ys.append(y.copy())
    if not store:
        ts.append(tf)
        ys.append(y.copy())
    return Solution(np.array(ts), np.array(ys), nfev=1 + 6 * n_steps,
                    n_accept=n_steps)


# --------------------------------------------------------------------------- #
# Adaptive Dormand-Prince 5(4)
# --------------------------------------------------------------------------- #

@dataclass
class Event:
    """Root of g(t, y) = 0 to locate during integration.

    direction: +1 for g crossing upward, -1 downward, 0 either way.
    terminal: stop the integration at the first root.
    """
    g: Callable[[float, np.ndarray], float]
    direction: int = 0
    terminal: bool = False


def _error_norm(err, y_old, y_new, atol, rtol) -> float:
    sc = atol + rtol * np.maximum(np.abs(y_old), np.abs(y_new))
    return float(np.sqrt(np.mean((err / sc) ** 2)))


def _initial_step(f, t0, y0, f0, tf, atol, rtol, order=5) -> float:
    """Starting step size (Hairer, Norsett & Wanner, Sec. II.4)."""
    sc = atol + rtol * np.abs(y0)
    d0 = float(np.sqrt(np.mean((y0 / sc) ** 2)))
    d1 = float(np.sqrt(np.mean((f0 / sc) ** 2)))
    h0 = 1e-6 if (d0 < 1e-5 or d1 < 1e-5) else 0.01 * d0 / d1
    h0 = min(h0, abs(tf - t0))
    f1 = f(t0 + h0, y0 + h0 * f0)
    d2 = float(np.sqrt(np.mean(((f1 - f0) / sc) ** 2))) / h0
    if max(d1, d2) <= 1e-15:
        h1 = max(1e-6, h0 * 1e-3)
    else:
        h1 = (0.01 / max(d1, d2)) ** (1.0 / (order + 1))
    return min(100.0 * h0, h1, abs(tf - t0))


def dp5(f: RHS, t_span: Sequence[float], y0: Sequence[float], *,
        rtol: float = 1e-8, atol: float = 1e-10, h0: Optional[float] = None,
        h_min: float = 1e-14, h_max: float = math.inf,
        max_steps: int = 2_000_000, safety: float = 0.9,
        fac_min: float = 0.2, fac_max: float = 5.0,
        positive: Optional[int] = None, events: Sequence[Event] = (),
        store: bool = True) -> Solution:
    """Adaptive Dormand-Prince 5(4) integration (report Algorithm 2).

    Steps are accepted when the scaled RMS error norm E <= 1. The next step is
    h * min(fac_max, max(fac_min, safety * E^(-1/5))), capped at 1 after a
    rejection. A non-finite step, or (if ``positive`` is an index) one that
    makes y[positive] <= 0, is rejected and retried with h/4. Pass
    positive=0 for the bubble radius. Event roots are found by re-stepping from the last
    accepted point (bracketed secant / bisection), so they inherit the
    method's fifth-order accuracy.
    """
    t0, tf = map(float, t_span)
    y = np.array(y0, dtype=float)
    k1 = f(t0, y)
    nfev = 1
    if h0 is None:
        h = _initial_step(f, t0, y, k1, tf, atol, rtol)
        nfev += 1
    else:
        h = float(h0)
    h = min(h, h_max)

    t = t0
    ts, ys = [t0], [y.copy()]
    sol = Solution(np.empty(0), np.empty(0), nfev=0)
    g_old = [ev.g(t, y) for ev in events]
    n_acc = n_rej = 0
    status, message = "success", ""
    rejected_last = False

    while t < tf:
        if n_acc + n_rej >= max_steps:
            status, message = "failed", "max_steps reached"
            break
        h = min(h, tf - t)
        if tf - (t + h) < 1e-12 * max(1.0, abs(tf)):
            h = tf - t  # land exactly on tf, avoid a sliver of a step
        # Trial stages may overshoot to R < 0 and produce NaN; such a step is
        # rejected below, so the floating-point warning is expected noise.
        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            y_new, err, k7 = dp5_step(f, t, y, h, k1)
        nfev += 6

        bad = not np.all(np.isfinite(y_new)) or (
            positive is not None and y_new[positive] <= 0.0)
        if bad:
            E = math.inf
        else:
            E = _error_norm(err, y, y_new, atol, rtol)

        if E <= 1.0:
            t_new = t + h
            # ---- events on the accepted step [t, t_new] ----
            stop = False
            for j, ev in enumerate(events):
                g_new = ev.g(t_new, y_new)
                crossed = (g_old[j] < 0 <= g_new and ev.direction >= 0) or (
                    g_old[j] > 0 >= g_new and ev.direction <= 0)
                if crossed and g_old[j] != 0.0:
                    te, ye, extra = _locate(f, ev.g, t, y, k1, h, g_old[j], g_new)
                    nfev += extra
                    sol.t_events.append((j, te))
                    sol.y_events.append((j, ye))
                    if ev.terminal:
                        t_new, y_new, stop = te, ye, True
                        break
                g_old[j] = g_new
            t, y, k1 = t_new, y_new, k7
            n_acc += 1
            if store or stop:
                ts.append(t)
                ys.append(y.copy())
            if stop:
                status, message = "event", "terminal event"
                break
            fac = fac_max if E == 0.0 else safety * E ** (-0.2)
            fac = min(fac_max, max(fac_min, fac))
            if rejected_last:
                fac = min(fac, 1.0)
            h = min(h * fac, h_max)
            rejected_last = False
        else:
            n_rej += 1
            rejected_last = True
            if bad:
                h *= 0.25
            else:
                h *= min(1.0, max(fac_min, safety * E ** (-0.2)))
            if h < h_min:
                status = "failed"
                message = f"step size fell below h_min at t = {t:.6g}"
                break

    if not store and (not ts or ts[-1] != t):
        ts.append(t)
        ys.append(y.copy())
    sol.t, sol.y = np.array(ts), np.array(ys)
    sol.nfev, sol.n_accept, sol.n_reject = nfev, n_acc, n_rej
    sol.status, sol.message = status, message
    return sol


def _locate(f, g, t, y, k1, h, g0, g1, tol=1e-13, max_iter=60):
    """Find s in (0, h] with g(t+s, step(y, s)) = 0 by the Illinois method.

    Each trial re-steps the Dormand-Prince formula from (t, y) with step s,
    which is as accurate as an ordinary step. Returns (t_event, y_event, nfev).
    """
    a, b, ga, gb = 0.0, h, g0, g1
    nfev = 0
    side = 0
    ys = None
    s = b
    for _ in range(max_iter):
        s = (a * gb - b * ga) / (gb - ga)
        if not (a < s < b):
            s = 0.5 * (a + b)
        ys, _, _ = dp5_step(f, t, y, s, k1)
        nfev += 6
        gs = g(t + s, ys)
        if gs == 0.0 or (b - a) < tol * max(1.0, abs(t)):
            break
        if (gs > 0) == (gb > 0):
            b, gb = s, gs
            if side == -1:
                ga *= 0.5
            side = -1
        else:
            a, ga = s, gs
            if side == 1:
                gb *= 0.5
            side = 1
        if abs(gs) < 1e-15:
            break
    return t + s, ys, nfev


def stability_polynomial(method: str) -> Callable[[complex], complex]:
    """R(z) for y' = lambda y, z = h lambda (report Section 2.5)."""
    if method == "rk4":
        A, b = RK4_A, RK4_B
    elif method == "dp5":
        A, b = DP_A, DP_B
    else:
        raise ValueError(method)
    s = len(b)
    eye = np.eye(s)
    ones = np.ones(s)

    def R(z: complex) -> complex:
        # R(z) = 1 + z b^T (I - z A)^{-1} 1, exact for explicit methods
        return 1.0 + z * (b @ np.linalg.solve(eye - z * A, ones))

    return R
