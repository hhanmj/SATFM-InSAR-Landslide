# Synthetic Benchmark Outputs

This folder holds the expected summary for the manuscript Table 2 synthetic
benchmark.

- `manuscript_table2_expected.csv`: expected ENU RMSE summary (mean +/- standard
  deviation, mm/yr) produced by `scripts/run_manuscript_table2.py`.

The synthetic scene, multi-geometry LOS observations, and noise used for Table 2
are generated deterministically in code (`src/satfm/manuscript_table2.py`) with
fixed random seeds; no separate synthetic data file is distributed. Running the
reproduction script regenerates all reported values on a `150 x 150` grid with a
`15 m` pixel size.

Do not place real SAR, GNSS, UAV, borehole, DEM, or shapefile data in this
folder unless the release rights have been checked.
