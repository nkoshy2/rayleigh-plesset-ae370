"""Rayleigh-Plesset model in dimensionless form, plus the three test cases.

Dimensional equation (report Eq. 1):

    R R'' + (3/2) R'^2 = (1/rho) [ p_v + p_g0 (R0/R)^(3 kappa) - p_inf(t)
                                   - 2 sigma / R - 4 mu R' / R ]

Scaling R = R0 * x and t = t_R * s with t_R = R0 sqrt(rho / dp) gives
(report Eq. 12), with primes now d/ds,

    x x'' + (3/2) x'^2 = Pg x^(-3 kappa) - P(s) - S / x - (4 / Re) x' / x

where Pg = p_g0 / dp, P(s) = (p_inf - p_v) / dp, S = 2 sigma / (R0 dp) and
Re = R0 sqrt(rho dp) / mu. We store inv_Re = 1 / Re so that mu = 0 is exact.

The integrators advance y = [x, x'] with ``rhs(t, y, prm)``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

# Water at 20 C, 1 atm (NIST WebBook, checked 2026-10-21). See report Table 1.
RHO = 998.2          # kg/m^3
MU = 1.0016e-3       # Pa s
SIGMA = 0.07274      # N/m
P_V = 2.339e3        # Pa
P_ATM = 101325.0     # Pa
C_SOUND = 1482.0     # m/s, only used for the wall Mach number diagnostic


@dataclass(frozen=True)
class Params:
    """Dimensionless parameters of Eq. 12 plus the scales to convert back."""

    Pg: float
    S: float
    inv_Re: float
    kappa: float
    P: Callable[[float], float]
    # Scales (for converting results back to SI); not used by rhs().
    R0: float = 1.0
    t_R: float = 1.0
    dp: float = 1.0
    label: str = ""
    extra: dict = field(default_factory=dict)

    @property
    def u_R(self) -> float:
        """Velocity scale R0 / t_R in m/s."""
        return self.R0 / self.t_R


def rhs(t: float, y: np.ndarray, prm: Params) -> np.ndarray:
    """Dimensionless right-hand side f(t, y) for y = [x, x'] (report Eq. rhs)."""
    x, v = y[0], y[1]
    g = (prm.Pg * x ** (-3.0 * prm.kappa) - prm.P(t) - prm.S / x
         - 4.0 * prm.inv_Re * v / x - 1.5 * v * v)
    return np.array([v, g / x])


def jacobian(t: float, y: np.ndarray, prm: Params) -> np.ndarray:
    """Analytic Jacobian df/dy, used for the stability (eigenvalue) check."""
    x, v = y[0], y[1]
    k3 = 3.0 * prm.kappa
    g = (prm.Pg * x ** (-k3) - prm.P(t) - prm.S / x
         - 4.0 * prm.inv_Re * v / x - 1.5 * v * v)
    dg_dx = -k3 * prm.Pg * x ** (-k3 - 1.0) + prm.S / x**2 + 4.0 * prm.inv_Re * v / x**2
    dg_dv = -4.0 * prm.inv_Re / x - 3.0 * v
    return np.array([[0.0, 1.0],
                     [(dg_dx * x - g) / x**2, dg_dv / x]])


def rhs_dimensional(t: float, y: np.ndarray, *, rho, mu, sigma, p_v, p_g0,
                    kappa, R0, p_inf) -> np.ndarray:
    """SI right-hand side of Eq. 1, used only to test the nondimensionalization."""
    R, Rd = y[0], y[1]
    p_wall = (p_v + p_g0 * (R0 / R) ** (3.0 * kappa) - p_inf(t)
              - 2.0 * sigma / R - 4.0 * mu * Rd / R)
    return np.array([Rd, (p_wall / rho - 1.5 * Rd * Rd) / R])


# --------------------------------------------------------------------------- #
# Test cases (report Section 1.4)
# --------------------------------------------------------------------------- #

def _const(value: float) -> Callable[[float], float]:
    return lambda t: value


def case_A(R0: float = 1e-3, p_inf: float = P_ATM, rho: float = RHO) -> Params:
    """Empty-cavity Rayleigh collapse: p_g0 = p_v = sigma = mu = 0."""
    t_R = R0 * math.sqrt(rho / p_inf)
    return Params(Pg=0.0, S=0.0, inv_Re=0.0, kappa=1.4, P=_const(1.0),
                  R0=R0, t_R=t_R, dp=p_inf, label="Case A")


def case_B(R0: float = 100e-6, Pi: float = 5.0, kappa: float = 1.4,
           p0: float = P_ATM, rho: float = RHO, mu: float = MU,
           sigma: float = SIGMA, p_v: float = P_V) -> Params:
    """Gas bubble in equilibrium at p0, ambient jumps to p1 = Pi p0 at t = 0."""
    p1 = Pi * p0
    dp = p1 - p_v
    p_g0 = p0 - p_v + 2.0 * sigma / R0
    t_R = R0 * math.sqrt(rho / dp)
    inv_Re = mu / (R0 * math.sqrt(rho * dp))
    return Params(Pg=p_g0 / dp, S=2.0 * sigma / (R0 * dp), inv_Re=inv_Re,
                  kappa=kappa, P=_const(1.0), R0=R0, t_R=t_R, dp=dp,
                  label=f"Case B (Pi={Pi:g})",
                  extra=dict(Pi=Pi, p0=p0, p1=p1, p_g0=p_g0, p_v=p_v))


def linear_frequency(prm: Params) -> tuple[float, float]:
    """Dimensionless natural frequency w0 t_R and damping beta t_R (Eq. 10)."""
    w0 = math.sqrt(3.0 * prm.kappa * prm.Pg - prm.S)
    beta = 2.0 * prm.inv_Re
    return w0, beta


def case_C(R0: float = 10e-6, eps: float = 0.01, omega_ratio: float = 1.0,
           kappa: float = 1.4, p0: float = P_ATM, rho: float = RHO,
           mu: float = MU, sigma: float = SIGMA, p_v: float = P_V) -> Params:
    """Equilibrium at p0, then p_inf = p0 [1 - eps sin(w t)] with w = omega_ratio * w0."""
    dp = p0 - p_v
    p_g0 = p0 - p_v + 2.0 * sigma / R0
    t_R = R0 * math.sqrt(rho / dp)
    inv_Re = mu / (R0 * math.sqrt(rho * dp))
    Pg, S = p_g0 / dp, 2.0 * sigma / (R0 * dp)
    w0 = math.sqrt(3.0 * kappa * Pg - S)
    Om = omega_ratio * w0
    amp = eps * p0 / dp

    def P(t: float) -> float:
        # (p_inf - p_v)/dp = 1 - eps p0 sin(Om t) / dp
        return 1.0 - amp * math.sin(Om * t)

    return Params(Pg=Pg, S=S, inv_Re=inv_Re, kappa=kappa, P=P, R0=R0, t_R=t_R,
                  dp=dp, label=f"Case C (eps={eps:g}, w/w0={omega_ratio:g})",
                  extra=dict(eps=eps, Omega=Om, w0=w0, forcing=amp, p0=p0,
                             p_g0=p_g0))


def with_terms(prm: Params, *, gas: bool = True, viscosity: bool = True,
               surface_tension: bool = True) -> Params:
    """Copy of prm with selected physics switched off (for Q3)."""
    return Params(Pg=prm.Pg if gas else 0.0,
                  S=prm.S if surface_tension else 0.0,
                  inv_Re=prm.inv_Re if viscosity else 0.0,
                  kappa=prm.kappa, P=prm.P, R0=prm.R0, t_R=prm.t_R, dp=prm.dp,
                  label=prm.label, extra=dict(prm.extra))


# --------------------------------------------------------------------------- #
# Reference quantities used for verification (report Section 3)
# --------------------------------------------------------------------------- #

RAYLEIGH_CONSTANT = math.sqrt(1.5) * math.gamma(5 / 6) * math.gamma(0.5) / (
    3.0 * math.gamma(5 / 6 + 0.5))  # 0.914681...


def rayleigh_time(x: float, n: int = 512) -> float:
    """Exact dimensionless time for an empty cavity to reach x = R/R0 (Eq. 5).

    t(x) = sqrt(3/2) * int_x^1 u^(3/2) / sqrt(1 - u^3) du. The substitution
    u = 1 - s^2 removes the endpoint singularity, leaving a smooth integrand
    that Gauss-Legendre quadrature handles to machine precision.
    """
    if not 0.0 <= x <= 1.0:
        raise ValueError("x must be in [0, 1]")
    smax = math.sqrt(1.0 - x)
    nodes, weights = np.polynomial.legendre.leggauss(n)
    s = 0.5 * smax * (nodes + 1.0)
    u = 1.0 - s * s
    integrand = u**1.5 * 2.0 / np.sqrt(3.0 - 3.0 * s * s + s**4)
    return math.sqrt(1.5) * 0.5 * smax * float(np.dot(weights, integrand))


def energy(y: np.ndarray, prm: Params) -> np.ndarray:
    """Dimensionless energy of Eq. 8 divided by (4 pi / 3) R0^3 dp.

    Conserved when inv_Re = 0, P is constant and kappa != 1. Accepts y with
    shape (2,) or (N, 2).
    """
    y = np.asarray(y)
    x, v = y[..., 0], y[..., 1]
    P = prm.P(0.0)
    return (1.5 * x**3 * v**2 + P * x**3 + 1.5 * prm.S * x**2
            + prm.Pg * x ** (3.0 * (1.0 - prm.kappa)) / (prm.kappa - 1.0))


def linear_amplitude(prm: Params) -> float:
    """Steady amplitude of x = R/R0 - 1 from the linearized Eq. 10 (Case C)."""
    w0, beta = linear_frequency(prm)
    Om, F = prm.extra["Omega"], prm.extra["forcing"]
    return F / math.sqrt((w0**2 - Om**2) ** 2 + (2.0 * beta * Om) ** 2)
