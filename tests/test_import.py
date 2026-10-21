"""Smoke test: the package imports and exposes the RHS entry point."""

import rp.model


def test_rhs_exists():
    assert callable(rp.model.rhs)
