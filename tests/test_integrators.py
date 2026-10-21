"""Tests of the integrators on problems with known answers (no bubble physics)."""

import math
from fractions import Fraction as Fr

import numpy as np
import pytest

from rp.integrators import (DP_A, DP_B, DP_BSTAR, DP_C, RK4_A, RK4_B, RK4_C,
                            Event, dp5, dp5_fixed, rk4, stability_polynomial)


# ---------------------------------------------------------------- tableaus --

def order_conditions(A, b, c):
    """Residuals of the eight order conditions through order 4 (report Eq. order)."""
    return np.array([
        b.sum() - 1,
        b @ c - 1 / 2,
        b @ c**2 - 1 / 3,
        b @ A @ c - 1 / 6,
        b @ c**3 - 1 / 4,
        (b * c) @ A @ c - 1 / 8,
        b @ A @ c**2 - 1 / 12,
        b @ A @ A @ c - 1 / 24,
    ])


@pytest.mark.parametrize("A,b,c", [(RK4_A, RK4_B, RK4_C), (DP_A, DP_B, DP_C),
                                   (DP_A, DP_BSTAR, DP_C)])
def test_order_conditions_through_4(A, b, c):
    assert np.allclose(A.sum(axis=1), c, atol=1e-14)
    assert np.max(np.abs(order_conditions(A, b, c))) < 1e-14


def test_dp5_fifth_order_conditions():
    """The b weights also satisfy the nine order-5 conditions; b* does not."""
    A, c = DP_A, DP_C

    def res(b):
        return np.array([
            b @ c**4 - 1 / 5,
            (b * c**2) @ A @ c - 1 / 10,
            b @ (A @ c) ** 2 - 1 / 20,
            (b * c) @ A @ c**2 - 1 / 15,
            b @ A @ c**3 - 1 / 20,
            (b * c) @ A @ A @ c - 1 / 30,
            b @ A @ (c * (A @ c)) - 1 / 40,
            b @ A @ A @ c**2 - 1 / 60,
            b @ A @ A @ A @ c - 1 / 120,
        ])

    assert np.max(np.abs(res(DP_B))) < 1e-14
    assert np.max(np.abs(res(DP_BSTAR))) > 1e-6


def test_fsal_last_row_equals_b():
    assert np.allclose(DP_A[6], DP_B)


# ---------------------------------------------------------- stability poly --

def test_stability_polynomials_match_report():
    R4 = stability_polynomial("rk4")
    R5 = stability_polynomial("dp5")
    for z in [-1.0, -2.5 + 0.3j, 1.2j, 0.4 - 0.7j]:
        assert abs(R4(z) - (1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24)) < 1e-13
        expect = (1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24 + z**5 / 120 + z**6 / 600)
        assert abs(R5(z) - expect) < 1e-13


def test_real_axis_stability_limits():
    R4 = stability_polynomial("rk4")
    R5 = stability_polynomial("dp5")
    assert abs(R4(-2.78)) <= 1 < abs(R4(-2.79))
    assert abs(R5(-3.30)) <= 1 < abs(R5(-3.31))


# ------------------------------------------------------- convergence orders --

def f_test(t, y):
    # y' = -2 t y^2, y(0) = 1  ->  y = 1 / (1 + t^2). Nonlinear and non-autonomous.
    return np.array([-2.0 * t * y[0] ** 2])


def exact_test(t):
    return 1.0 / (1.0 + t * t)


def observed_order(solver, ns):
    """Order from the two finest runs (coarser ones can be pre-asymptotic)."""
    errs = [abs(solver(f_test, (0.0, 2.0), [1.0], n).y[-1, 0] - exact_test(2.0))
            for n in ns]
    return math.log(errs[-1] / errs[-2]) / math.log(ns[-1] / ns[-2])


def test_rk4_is_fourth_order():
    assert observed_order(rk4, [80, 160]) == pytest.approx(-4.0, abs=0.15)


def test_dp5_fixed_is_fifth_order():
    assert observed_order(dp5_fixed, [80, 160]) == pytest.approx(-5.0, abs=0.3)


# ---------------------------------------------------------------- adaptive --

def test_dp5_adaptive_meets_tolerance_and_lands_on_tf():
    for rtol in [1e-6, 1e-9, 1e-12]:
        sol = dp5(f_test, (0.0, 2.0), [1.0], rtol=rtol, atol=rtol * 1e-2)
        assert sol.status == "success"
        assert sol.t[-1] == 2.0
        err = abs(sol.y[-1, 0] - exact_test(2.0))
        assert err < 50 * rtol


def test_dp5_tolerance_proportionality():
    errs, tols = [], [1e-5, 1e-7, 1e-9, 1e-11]
    for tol in tols:
        sol = dp5(f_test, (0.0, 2.0), [1.0], rtol=tol, atol=tol)
        errs.append(abs(sol.y[-1, 0] - exact_test(2.0)))
    slope = np.polyfit(np.log(tols), np.log(errs), 1)[0]
    assert 0.7 < slope < 1.3


def test_dp5_harmonic_oscillator_long_run():
    f = lambda t, y: np.array([y[1], -y[0]])
    T = 20 * 2 * math.pi
    sol = dp5(f, (0.0, T), [1.0, 0.0], rtol=1e-10, atol=1e-12)
    assert np.allclose(sol.y[-1], [1.0, 0.0], atol=1e-7)


def test_terminal_event_location():
    # y' = -1, y(0) = 1: crosses y = 0.3 at t = 0.7 exactly
    f = lambda t, y: np.array([-1.0])
    ev = Event(lambda t, y: y[0] - 0.3, direction=-1, terminal=True)
    sol = dp5(f, (0.0, 5.0), [1.0], rtol=1e-8, events=[ev], positive=None)
    assert sol.status == "event"
    assert sol.t[-1] == pytest.approx(0.7, abs=1e-12)


def test_rejects_steps_that_cross_zero():
    # y' = -sqrt(y) style collapse; y must stay positive: y' = -1/y, y(0)=1,
    # exact y = sqrt(1 - 2t) hits 0 at t = 0.5.
    f = lambda t, y: np.array([-1.0 / y[0]])
    sol = dp5(f, (0.0, 1.0), [1.0], rtol=1e-8, h_min=1e-12, positive=0)
    assert sol.status == "failed"
    assert np.all(sol.y[:, 0] > 0)
    assert sol.t[-1] == pytest.approx(0.5, abs=1e-6)
