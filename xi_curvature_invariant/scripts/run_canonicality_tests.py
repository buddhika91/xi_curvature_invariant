#!/usr/bin/env python3
"""Run chain-normalized reparameterization tests for K0."""
from __future__ import annotations

import argparse
from pathlib import Path

import sys
from pathlib import Path as _Path
sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))


from xi_curvature.scans import canonicality_scan, h_stability, gamma_window_scan


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--gamma", type=str, default="0.0001")
    parser.add_argument("--h", type=str, default="1e-5")
    parser.add_argument("--outdir", type=str, default="results")
    parser.add_argument("--window-points", type=int, default=120, help="points per window for gamma scan")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    can = canonicality_scan(gamma=args.gamma, h=args.h, dps=args.dps)
    hstab = h_stability(gamma="0", dps=args.dps)
    windows = gamma_window_scan(points=args.window_points, h=args.h, dps=max(50, args.dps - 10))

    can.to_csv(outdir / "canonicality_tests.csv", index=False)
    hstab.to_csv(outdir / "h_stability.csv", index=False)
    windows.to_csv(outdir / "gamma_window_scan.csv", index=False)

    print("CANONICALITY TESTS")
    print(can.to_string(index=False))
    print("\nH-STABILITY")
    print(hstab.to_string(index=False))
    print("\nGAMMA WINDOW SCAN")
    print(windows.to_string(index=False))
    print(f"\nWrote CSV files to: {outdir.resolve()}")


if __name__ == "__main__":
    main()
