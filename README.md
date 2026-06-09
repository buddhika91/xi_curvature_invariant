# xi-curvature-invariant

Numerical support package for the paper:

**A Canonical Curvature Invariant of the Riemann Ξ-Function and Its Numerical Proximity to the Fine-Structure Constant**

This repository computes and tests the candidate curvature invariant

```math
K_0 = \frac{\Xi''(1/2)}{\Xi(1/2)}
```

for the completed Riemann function

```math
\Xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s).
```

It also computes the associated dimensionless spectral count

```math
N_\Xi=\frac{2\pi}{K_0}.
```

Current reference values:

```text
K0   ≈ 0.04620998623...
N_Xi ≈ 135.97029169...
```

## What this code does

The package provides reproducible scripts for:

1. Computing `Xi(1/2)`, `Xi''(1/2)`, `K0`, and `N_Xi`.
2. Evaluating the transverse curvature function

```math
K(\gamma)=\frac12\left.\partial_\sigma^2
\log\left|\Xi\left(\frac12+\sigma+i\gamma\right)\right|^2\right|_{\sigma=0}.
```

3. Testing finite-difference stability in the transverse step size `h`.
4. Testing coordinate canonicality under linear and nonlinear reparameterizations.
5. Running coarse gamma-window scans for the curvature-floor conjecture

```math
K(\gamma) \ge K(0).
```

## Installation

```bash
git clone xi_curvature_invariant
cd xi-curvature-invariant
python -m venv .venv
.venv\Scripts\activate   # Windows
# source .venv/bin/activate  # macOS/Linux
pip install -e .[test]
```

## Quick start

Compute the central invariant:

```bash
python scripts/compute_k0.py --dps 100
```

Run canonicality and stability tests:

```bash
python scripts/run_canonicality_tests.py --dps 80 --gamma 0.0001 --h 1e-5
```

Create a simple plot of `K(gamma)`:

```bash
python scripts/plot_curvature_scan.py --gamma-max 5 --points 200
```

Run tests:

```bash
pytest -q
```

## Expected output

`compute_k0.py` should report values close to:

```text
K0 = Xi''/Xi  : approximately 0.04620998623...
N_Xi = 2pi/K0 : approximately 135.97029169...
```

`run_canonicality_tests.py` writes CSV files into `results/`:

```text
results/canonicality_tests.csv
results/h_stability.csv
results/gamma_window_scan.csv
```

## Repository structure

```text
xi-curvature-invariant/
├── xi_curvature/
│   ├── __init__.py
│   ├── core.py              # Xi, landscape, curvature, K0
│   └── scans.py             # h-stability, canonicality, gamma scans
├── scripts/
│   ├── compute_k0.py
│   ├── run_canonicality_tests.py
│   └── plot_curvature_scan.py
├── tests/
│   └── test_core.py
├── results/
├── pyproject.toml
├── requirements.txt
├── LICENSE
└── README.md
```

## Scientific status

This repository supports a focused mathematical note. It does **not** claim to prove a physical derivation of the fine-structure constant.

Established computationally here:

```math
K_0 = \Xi''(1/2)/\Xi(1/2)
```

and

```math
N_\Xi = 2\pi/K_0 = 135.97029169\ldots.
```

Open mathematical problem:

```math
K(\gamma)\ge K(0)\quad\text{for all real }\gamma.
```

The numerical proximity between `N_Xi` and `alpha^{-1}` is recorded as an observation only.

## Citation

If you use this code, cite the accompanying paper/note when available.


## Numerical note

The finite-difference curvature scan near `gamma=0` and the direct numerical derivative `Xi''(1/2)/Xi(1/2)` agree at the scale relevant for the paper. Very small differences in the final displayed digits can occur depending on whether the value is taken from a finite-difference scan at `gamma=1e-4`, exact `gamma=0`, or direct high-order differentiation. The scripts report the method used so these values are reproducible.
