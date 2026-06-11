# SATFM: Spatially Adaptive Terrain-Kinematic Constraints

This repository contains a clean, public-facing implementation of the SATFM
synthetic benchmark workflow used in:

> Spatially Adaptive Terrain-Kinematic Constraints for Resolving North-South
> Landslide Motion from Multi-Geometry InSAR in the Three Gorges Reservoir Area

The release is intentionally limited to code and synthetic data. It does not
include proprietary or access-restricted SAR products, GNSS observations, UAV
data, borehole data, or intermediate real-case outputs.

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
WLS        8.903 +/- 0.030 mm/yr
APFM      12.825 +/- 0.034 mm/yr
SPFM      8.641 +/- 0.026 mm/yr
Full SATFM 2.122 +/- 0.027 mm/yr
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

Please cite the associated manuscript when using this code or synthetic dataset.
The manuscript has not yet been assigned a DOI; citation metadata can be added
after acceptance or publication.
