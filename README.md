# SATFM: Spatially Adaptive Terrain-Kinematic Constraints

This repository contains a clean, public-facing implementation of the SATFM
synthetic benchmark workflow used in:

> Cross-looking InSAR reconstruction of three-dimensional kinematics and
> moving-layer geometry of a reservoir-bank landslide

Availability of the observational and site-specific datasets supporting the
paper is described in its Data Availability statement.

## What Is Included

- SATFM core solver with spatially adaptive terrain-kinematic constraints.
- Baseline solvers: unconstrained weighted least squares, APFM, and SPFM.
- Synthetic landslide scene generator with known 3-D velocity truth.
- Multi-geometry LOS observation simulator.
- Manuscript Table 2 reproduction script for the 15 m, 150 x 150, 30-realization synthetic benchmark.

## Repository Layout

```text
configs/                  Parameter files for benchmark runs
data/synthetic/           Expected Table 2 summary and data notes
docs/                     Method notes and data policy
scripts/                  Command-line entry points
src/satfm/                SATFM Python package
tests/                    Smoke tests for installation and solver health
```

## Installation

Create a clean Python environment, then install the listed dependencies:

```bash
pip install -r requirements.txt
```

Python 3.10 or newer is recommended. No package installation step is required:
the command-line scripts add `src/` to the Python path at runtime.

## Reproduce Manuscript Table 2

Run the 15 m synthetic benchmark reported in manuscript Table 2:

```bash
python scripts/run_manuscript_table2.py --output outputs/manuscript_table2
```

This runs 30 Monte Carlo realizations on a `150 x 150` grid with a `15 m`
pixel size and the manuscript Eq. (12) edge kernel (`legacy_linear` with
direction-confidence gating). The expected North-component RMSEs are:

```text
WLS        8.902 +/- 0.036 mm/yr
APFM      12.817 +/- 0.034 mm/yr
SPFM      8.637 +/- 0.024 mm/yr
Full SATFM 2.119 +/- 0.027 mm/yr
```

The expected summary is also provided in
`data/synthetic/manuscript_table2_expected.csv`.

## Method Notes

SATFM solves the 3-D velocity field from multi-geometry LOS observations by
combining:

1. local weighted observation equations,
2. a terrain-direction weak constraint,
3. a terrain-parallel weak constraint,
4. spatially adaptive edge-aware smoothing.

The manuscript Table 2 reproduction uses the Eq. (12) edge kernel
(`legacy_linear`) with direction-confidence gating, matching the reported
experiments.

## Citation

If you use this code or synthetic benchmark, please cite:

> Han, M., Dou, J., Yang, Y., et al. *Cross-looking InSAR reconstruction of
> three-dimensional kinematics and moving-layer geometry of a reservoir-bank
> landslide*.
