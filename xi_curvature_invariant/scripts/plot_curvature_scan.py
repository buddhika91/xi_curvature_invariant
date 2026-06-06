#!/usr/bin/env python3
"""Create a simple K(gamma) plot over a user-defined range."""
from __future__ import annotations

import argparse
from pathlib import Path

import sys
from pathlib import Path as _Path
sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))


import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt

from xi_curvature.core import curvature_fd, k0_analytic, set_precision


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--gamma-max", type=float, default=5.0)
    parser.add_argument("--points", type=int, default=200)
    parser.add_argument("--h", type=str, default="1e-5")
    parser.add_argument("--out", type=str, default="results/curvature_scan.png")
    args = parser.parse_args()

    set_precision(args.dps)
    gammas = np.linspace(0.0001, args.gamma_max, args.points)
    ks = [float(curvature_fd(mp.mpf(str(g)), mp.mpf(args.h))) for g in gammas]
    k0 = float(k0_analytic())

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(8, 5))
    plt.plot(gammas, ks, label="K(gamma)")
    plt.axhline(k0, linestyle="--", label="K0 analytic")
    plt.xlabel("gamma")
    plt.ylabel("K(gamma)")
    plt.title("Transverse Curvature of the Riemann Xi Log-Landscape")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=160)
    print(f"Saved plot: {out.resolve()}")


if __name__ == "__main__":
    main()
