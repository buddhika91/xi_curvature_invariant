from __future__ import annotations

import mpmath as mp

from xi_curvature.core import xi, xi_constants, k0_analytic, curvature_fd, set_precision


def test_functional_equation_center_symmetry():
    set_precision(50)
    z = mp.mpc("0.2", "1.3")
    left = xi(mp.mpf("0.5") + z)
    right = xi(mp.mpf("0.5") - z)
    assert abs(left - right) < mp.mpf("1e-40")


def test_k0_matches_reference_digits():
    c = xi_constants(100)
    ref = mp.mpf("0.046209986230837941577867")
    assert abs(c.k0 - ref) < mp.mpf("1e-24")


def test_n_xi_matches_reference_digits():
    c = xi_constants(100)
    ref = mp.mpf("135.970291698259424392746")
    assert abs(c.n_xi - ref) < mp.mpf("1e-21")


def test_finite_difference_close_to_k0():
    set_precision(80)
    k_fd = curvature_fd(mp.mpf("0"), mp.mpf("1e-5"))
    k0 = k0_analytic()
    # Finite difference has truncation error near 1e-12 at h=1e-5.
    assert abs(k_fd - k0) < mp.mpf("1e-10")
