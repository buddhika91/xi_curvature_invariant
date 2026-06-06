"""Scan utilities for curvature, h-stability, and canonicality tests."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import mpmath as mp
import pandas as pd

from .core import curvature_fd, set_precision, transformed_curvature_fd


@dataclass(frozen=True)
class TransformSpec:
    name: str
    func: Callable[[mp.mpf], mp.mpf]
    fprime0: mp.mpf


def default_transforms() -> list[TransformSpec]:
    """Coordinate transforms used in the canonicality attack."""
    return [
        TransformSpec("identity sigma=x", lambda x: x, mp.mpf("1")),
        TransformSpec("scale sigma=2x", lambda x: 2 * x, mp.mpf("2")),
        TransformSpec("scale sigma=0.5x", lambda x: mp.mpf("0.5") * x, mp.mpf("0.5")),
        TransformSpec("tanh sigma=tanh(x)", lambda x: mp.tanh(x), mp.mpf("1")),
        TransformSpec("sinh sigma=sinh(x)", lambda x: mp.sinh(x), mp.mpf("1")),
        TransformSpec("cubic sigma=x+0.1x^3", lambda x: x + mp.mpf("0.1") * x**3, mp.mpf("1")),
        TransformSpec("cubic sigma=x-0.1x^3", lambda x: x - mp.mpf("0.1") * x**3, mp.mpf("1")),
    ]


def h_stability(gamma: str = "0", hs: Iterable[str] | None = None, dps: int = 80) -> pd.DataFrame:
    set_precision(dps)
    if hs is None:
        hs = ["1e-3", "5e-4", "1e-4", "5e-5", "1e-5", "5e-6"]
    rows = []
    for h in hs:
        k = curvature_fd(mp.mpf(gamma), mp.mpf(h))
        rows.append({"h": h, "K": mp.nstr(k, 50), "2pi_over_K": mp.nstr(2 * mp.pi / k, 50)})
    return pd.DataFrame(rows)


def canonicality_scan(gamma: str = "0.0001", h: str = "1e-5", dps: int = 80) -> pd.DataFrame:
    set_precision(dps)
    base = curvature_fd(mp.mpf(gamma), mp.mpf(h))
    rows = []
    for spec in default_transforms():
        raw, normalized = transformed_curvature_fd(mp.mpf(gamma), spec.func, spec.fprime0, mp.mpf(h))
        rows.append(
            {
                "transform": spec.name,
                "fprime0": mp.nstr(spec.fprime0, 30),
                "raw_K": mp.nstr(raw, 50),
                "raw_2pi_over_K": mp.nstr(2 * mp.pi / raw, 50),
                "normalized_K": mp.nstr(normalized, 50),
                "normalized_2pi_over_K": mp.nstr(2 * mp.pi / normalized, 50),
                "normalized_minus_baseline": mp.nstr(normalized - base, 30),
            }
        )
    return pd.DataFrame(rows)


def gamma_window_scan(
    gamma_min: str = "0.0001",
    windows: Iterable[str] | None = None,
    points: int = 300,
    h: str = "1e-5",
    dps: int = 70,
) -> pd.DataFrame:
    """Coarse gamma-window scan for the minimum of K(gamma).

    This is intentionally simple and reproducible. For dense scans use smaller windows
    or increase points cautiously, since high-precision Xi evaluations are expensive.
    """
    set_precision(dps)
    if windows is None:
        windows = ["1", "5", "10", "30", "80"]
    rows = []
    g0 = mp.mpf(gamma_min)
    for w in windows:
        g1 = mp.mpf(w)
        best_g = None
        best_k = None
        for i in range(points):
            g = g0 + (g1 - g0) * i / max(points - 1, 1)
            k = curvature_fd(g, mp.mpf(h))
            if best_k is None or k < best_k:
                best_k = k
                best_g = g
        rows.append({"window_end": w, "gamma_min": mp.nstr(best_g, 40), "K_min": mp.nstr(best_k, 50), "2pi_over_K": mp.nstr(2 * mp.pi / best_k, 50)})
    return pd.DataFrame(rows)
