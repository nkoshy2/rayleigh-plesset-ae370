"""Tests of the Rayleigh-Plesset model and its verification references."""

import math

import numpy as np
import pytest

from rp.diagnostics import first_collapse
from rp.integrators import Event, dp5
from rp.model import (MU, P_ATM, P_V, RAYLEIGH_CONSTANT, RHO, SIGMA, case_A,
                      case_B, case_C, energy, jacobian, linear_amplitude,
                      rayleigh_time, rhs, rhs_dimensional, with_terms)


def test_rhs_exists():
    assert callable(rhs)


def test_rayleigh_constant():
    assert RAYLEIGH_CONSTANT == pytest.approx(0.914681356, abs=1e-9)
    assert rayleigh_time(0.0) == pytest.approx(RAYLEIGH_CONSTANT, abs=1e-13)
    assert rayleigh_time(1.0) == 0.0


def test_dimensionless_rhs_matches_dimensional():
    """f_SI = (R0 / t_R^2) * f_dimless at matching states (checks Eq. 12)."""
    prm = case_B(R0=50e-6, Pi=3.0, kappa=1.3)
    p0, p1 = prm.extra["p0"], prm.extra["p1"]
    for x, v in [(1.0, 0.0), (0.6, -0.8), (1.4, 0.3)]:
        f_nd = rhs(0.1, np.array([x, v]), prm)
        R, Rd = x * prm.R0, v * prm.u_R
        f_si = rhs_dimensional(0.1, np.array([R, Rd]), rho=RHO, mu=MU,
                               sigma=SIGMA, p_v=P_V, p_g0=prm.extra["p_g0"],
                               kappa=1.3, R0=prm.R0, p_inf=lambda t: p1)
        assert f_si[1] == pytest.approx(f_nd[1] * prm.R0 / prm.t_R**2, rel=1e-12)


def test_case_C_starts_in_equilibrium():
    prm = case_C()
    assert rhs(0.0, np.array([1.0, 0.0]), prm)[1] == pytest.approx(0.0, abs=1e-14)


def test_case_B_starts_collapsing():
    prm = case_B(Pi=5.0)
    acc = rhs(0.0, np.array([1.0, 0.0]), prm)[1]
    # x'' = (p0 - p1) / dp  in dimensionless form
    assert acc == pytest.approx((P_ATM - 5 * P_ATM) / prm.dp, rel=1e-12)


def test_jacobian_matches_finite_differences():
    prm = case_B(Pi=8.0)
    y = np.array([0.4, -1.7])
    J = jacobian(0.0, y, prm)
    eps = 1e-7
    for j in range(2):
        dy = np.zeros(2)
        dy[j] = eps
        col = (rhs(0.0, y + dy, prm) - rhs(0.0, y - dy, prm)) / (2 * eps)
        assert np.allclose(J[:, j], col, rtol=1e-6, atol=1e-6)


# ---------------------------------------------- verification against exact --

def test_case_A_matches_exact_rayleigh_solution():
    prm = case_A()
    f = lambda t, y: rhs(t, y, prm)
    for x_end in [0.5, 0.1, 0.02]:
        ev = Event(lambda t, y, xe=x_end: y[0] - xe, direction=-1, terminal=True)
        sol = dp5(f, (0.0, 2.0), [1.0, 0.0], rtol=1e-12, atol=1e-14, events=[ev],
                   positive=0)
        assert sol.status == "event"
        assert sol.t[-1] == pytest.approx(rayleigh_time(x_end), abs=1e-9)
        # first integral, Eq. 4 in dimensionless form: v^2 = (2/3)(1/x^3 - 1)
        x, v = sol.y[-1]
        assert v**2 == pytest.approx((2 / 3) * (x**-3 - 1), rel=1e-8)


def test_case_B_inviscid_energy_conserved():
    prm = with_terms(case_B(Pi=10.0), viscosity=False)
    f = lambda t, y: rhs(t, y, prm)
    sol = dp5(f, (0.0, 8.0), [1.0, 0.0], rtol=1e-11, atol=1e-13, positive=0)
    E = energy(sol.y, prm)
    assert np.max(np.abs(E - E[0])) / abs(E[0]) < 1e-8


def test_case_C_small_amplitude_matches_linear_theory():
    for ratio in [0.7, 1.0, 1.3]:
        prm = case_C(R0=10e-6, eps=1e-4, omega_ratio=ratio)
        f = lambda t, y: rhs(t, y, prm)
        Om = prm.extra["Omega"]
        n_per = 40
        T = (math.ceil(700 * Om / (2 * math.pi)) + n_per) * 2 * math.pi / Om
        sol = dp5(f, (0.0, T), [1.0, 0.0], rtol=1e-11, atol=1e-14,
                  h_max=0.05 * 2 * math.pi / Om, positive=0)
        tail = sol.t > T - n_per * 2 * math.pi / Om
        x = sol.y[tail, 0] - 1.0
        amp = 0.5 * (x.max() - x.min())
        assert amp == pytest.approx(linear_amplitude(prm), rel=5e-3)


def test_collapse_diagnostics_reasonable():
    c = first_collapse(case_B(Pi=5.0))
    assert 0.40 < c.x_min < 0.43
    assert c.v_min < 0
    assert c.mach < 0.1
    # stronger step -> deeper collapse, approaching Rayleigh timing
    c20 = first_collapse(case_B(Pi=20.0))
    assert c20.x_min < c.x_min
    assert c20.t_min < c.t_min
