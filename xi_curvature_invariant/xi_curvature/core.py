"""Core high-precision functions for the completed Riemann Xi-function.

The convention used here is
    Xi(s) = 1/2 * s * (s - 1) * pi^(-s/2) * Gamma(s/2) * zeta(s)

The transverse coordinate is centered at the critical line:
    s = 1/2 + sigma + i gamma.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import mpmath as mp


@dataclass(frozen=True)
class XiConstants:
    """Container for the central invariant values."""

    dps: int
    xi_half: mp.mpf
    xi_second_half: mp.mpf
    k0: mp.mpf
    n_xi: mp.mpf


def set_precision(dps: int = 80) -> None:
    """Set mpmath decimal precision."""
    if dps < 30:
        raise ValueError("dps should be at least 30 for stable Xi calculations.")
    mp.mp.dps = dps


def xi(s: complex | mp.mpc | mp.mpf) -> mp.mpc:
    """Completed Riemann Xi-function."""
    s = mp.mpc(s)
    return mp.mpf("0.5") * s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def centered_xi(z: complex | mp.mpc | mp.mpf) -> mp.mpc:
    """Xi(1/2 + z)."""
    return xi(mp.mpf("0.5") + z)


def log_landscape(sigma: mp.mpf, gamma: mp.mpf) -> mp.mpf:
    """L(sigma, gamma) = log |Xi(1/2 + sigma + i gamma)|^2."""
    z = mp.mpf("0.5") + mp.mpf(sigma) + 1j * mp.mpf(gamma)
    value = xi(z)
    return mp.log(abs(value) ** 2)


def curvature_fd(gamma: mp.mpf, h: mp.mpf = mp.mpf("1e-5")) -> mp.mpf:
    """Finite-difference transverse curvature K(gamma).

    K(gamma) = 1/2 * d^2/dsigma^2 log|Xi(1/2+sigma+i gamma)|^2 at sigma=0.
    The centered finite-difference formula already includes the 1/2 factor.
    """
    gamma = mp.mpf(gamma)
    h = mp.mpf(h)
    return (log_landscape(h, gamma) - 2 * log_landscape(0, gamma) + log_landscape(-h, gamma)) / (2 * h * h)


def curvature_log_derivative(gamma: mp.mpf) -> mp.mpf:
    """Curvature from the logarithmic derivative of F(z)=Xi(1/2+z).

    K(gamma) = Re[F''(i gamma)/F(i gamma) - (F'(i gamma)/F(i gamma))^2].
    """
    gamma = mp.mpf(gamma)
    z = 1j * gamma
    f = centered_xi(z)
    fp = mp.diff(centered_xi, z, 1)
    fpp = mp.diff(centered_xi, z, 2)
    return mp.re(fpp / f - (fp / f) ** 2)


def k0_analytic() -> mp.mpf:
    """K0 candidate at gamma=0: Xi''(1/2)/Xi(1/2)."""
    second = mp.diff(xi, mp.mpf("0.5"), 2)
    return mp.re(second / xi(mp.mpf("0.5")))


def xi_constants(dps: int = 80) -> XiConstants:
    """Compute Xi(1/2), Xi''(1/2), K0, and N_Xi at chosen precision."""
    set_precision(dps)
    half = mp.mpf("0.5")
    xi_half = xi(half)
    xi_second = mp.diff(xi, half, 2)
    k0 = mp.re(xi_second / xi_half)
    n_xi = 2 * mp.pi / k0
    return XiConstants(dps=dps, xi_half=xi_half, xi_second_half=xi_second, k0=k0, n_xi=n_xi)


def transformed_curvature_fd(
    gamma: mp.mpf,
    transform: Callable[[mp.mpf], mp.mpf],
    fprime0: mp.mpf,
    h: mp.mpf = mp.mpf("1e-5"),
) -> tuple[mp.mpf, mp.mpf]:
    """Raw and chain-normalized curvature under sigma = f(x).

    The finite difference is taken in x, while the landscape receives sigma=f(x).
    Raw curvature transforms by [f'(0)]^2. Normalized curvature divides by [f'(0)]^2.
    """
    gamma = mp.mpf(gamma)
    h = mp.mpf(h)
    fp0 = mp.mpf(fprime0)

    def lx(x: mp.mpf) -> mp.mpf:
        return log_landscape(transform(x), gamma)

    raw = (lx(h) - 2 * lx(mp.mpf("0")) + lx(-h)) / (2 * h * h)
    normalized = raw / (fp0 ** 2)
    return raw, normalized
