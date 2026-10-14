# Data Availability Notes

This public repository provides the source code and synthetic benchmark associated
with the paper. The synthetic benchmark reported in Table 2 is reproduced by
`scripts/run_manuscript_table2.py`. The synthetic scene, multi-geometry LOS
observations, and noise are generated deterministically in code using the fixed
seed range `2026100801`-`2026100830`; no separate synthetic data file is
distributed.

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

> The Sentinel-1 SAR data are freely available from the European Space Agency
> (https://dataspace.copernicus.eu). The SATFM source code and synthetic benchmark
> associated with this paper are openly available in this public GitHub
> repository. The Fucheng-1 SAR datasets, GNSS observation records, UAV-generated
> DEMs and DOMs, and borehole survey data can be obtained from the corresponding
> author upon reasonable request.
