"""Rayleigh-Plesset model: right-hand side of the first-order system.

State vector y = [R, Rdot]. The governing equation is

    R R'' + (3/2) R'^2 = (1/rho) [ p_v + p_g0 (R0/R)^(3 kappa) - p_inf(t)
                                   - 2 sigma / R - 4 mu R' / R ]

with p_g0 = p_inf(0) - p_v + 2 sigma / R0 unless a case overrides it
(Case A sets p_g0 = 0). See report/sections/premise.tex for symbols and units.
"""

from __future__ import annotations

import numpy as np


def rhs(t: float, y: np.ndarray, params) -> np.ndarray:
    """Return dy/dt = [R', R''] for the Rayleigh-Plesset system.

    Parameters
    ----------
    t : float
        Time [s] (or dimensionless time if params are dimensionless).
    y : ndarray, shape (2,)
        State [R, R'] in [m, m/s].
    params : object
        Holds rho, mu, sigma, p_v, p_g0, kappa, R0 and a callable p_inf(t).

    Returns
    -------
    ndarray, shape (2,)
        Time derivative of the state.
    """
    raise NotImplementedError("Model RHS not implemented yet.")
