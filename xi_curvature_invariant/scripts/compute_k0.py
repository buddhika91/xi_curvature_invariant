#!/usr/bin/env python3
"""Compute K0 = Xi''(1/2)/Xi(1/2) and N_Xi = 2pi/K0."""
from __future__ import annotations

import argparse
from pathlib import Path

import sys
from pathlib import Path as _Path
sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))


import mpmath as mp

from xi_curvature.core import xi_constants


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=100, help="mpmath decimal precision")
    parser.add_argument("--out", type=str, default="results/k0_values.txt", help="output text file")
    args = parser.parse_args()

    c = xi_constants(args.dps)
    text = "\n".join(
        [
            "RIEMANN XI CURVATURE INVARIANT",
            "=" * 80,
            f"dps                 : {c.dps}",
            f"Xi(1/2)             : {mp.nstr(c.xi_half, 80)}",
            f"Xi''(1/2)           : {mp.nstr(c.xi_second_half, 80)}",
            f"K0 = Xi''/Xi        : {mp.nstr(c.k0, 80)}",
            f"N_Xi = 2pi/K0       : {mp.nstr(c.n_xi, 80)}",
            "",
        ]
    )
    print(text)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
