# Data Availability Notes

This public repository is designed for synthetic reproducibility. The manuscript
Table 2 synthetic benchmark is reproduced by `scripts/run_manuscript_table2.py`.
The synthetic scene, multi-geometry LOS observations, and noise are generated
deterministically in code with fixed random seeds; no separate synthetic data
file is distributed.

## Publicly Included

- Source code for SATFM and baseline inversions.
- Synthetic scene generation code for the manuscript Table 2 benchmark.
- Synthetic LOS observation simulator.
- Manuscript Table 2 reproduction script and expected summary table.

## Not Included

- Sentinel-1 or other SAR products.
- GNSS observations.
- UAV-derived products.
- Borehole logs or interpreted slip-surface points.
- Real-case DEM derivatives, masks, shapefiles, or intermediate rasters.
- Any files copied from commercial, institutional, or restricted datasets.

## Manuscript Wording

Recommended wording for the manuscript:

> The SATFM source code required to reproduce the manuscript Table 2 synthetic
> benchmark is available in a public GitHub repository. The synthetic scene,
> multi-geometry LOS observations, and noise are generated deterministically by
> the released scripts with fixed random seeds. Real-case SAR, GNSS, UAV,
> borehole, and derived geospatial products are not redistributed because of
> data-provider and access restrictions; access details are provided in the
> manuscript data availability statement.
