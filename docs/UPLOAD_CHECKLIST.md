# Public Repository Checklist

Before updating the public repository:

1. Confirm that the MIT license terms and copyright holder are correct.
2. Confirm that no real SAR, GNSS, UAV, borehole, DEM, shapefile, raster, or
   institution-restricted product is committed.
3. From a clean environment, install dependencies:

```bash
pip install -r requirements.txt
```

4. Run the reproducibility checks:

```bash
python scripts/run_manuscript_table2.py --output outputs/manuscript_table2
python -m pytest -q
```

5. Confirm that `outputs/manuscript_table2/Table_2_maintext_ENU_RMSE_mean_pm_sd_mm_yr.csv`
   matches `data/synthetic/manuscript_table2_expected.csv`.
6. Check that `docs/METHOD_NOTES.md` matches the exact equations in the
   submitted manuscript.
7. Add the final GitHub URL to the manuscript Data Availability statement.
8. Add `CITATION.cff` only after the manuscript has stable citation metadata
   such as a DOI.
