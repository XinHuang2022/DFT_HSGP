# ERA5 North Sea reanalysis

This folder packages the North Sea ERA5 experiment used in the manuscript.
It contains the extracted six-hourly eastward 10-m wind component on the fixed
17-by-17 ocean grid, the established five spatial folds, the analytic and
fixed-rank Equation-9 implementations, and the frequency-adaptive HSGP
implementation.

## Entry points

- Interactive: `main_era5_reanalysis.m`
- Batch preflight: `run_era5_experiment("preflight", 1)`
- Analytic/fixed-rank paper workflow:
  `run_era5_experiment("fixed_paper", 5)`
- Frequency-adaptive paper workflow:
  `run_era5_experiment("adaptive_paper", 5)`

The interactive script defaults to `"retained"`, which reads the preserved
tables without starting an optimization. New runs are written beneath
`results/` in timestamped directories.

## Data contract

The local extracted dataset contains 360 six-hourly observations at 289
marine sites. The experiment models `u10_time_by_site.csv`. Site columns must
match `site_manifest_17x17.csv`, and all five folds must match the retained
fold manifest. Raw ERA5 NetCDF files and download software are intentionally
excluded from this analysis folder.

## Scientific workflows

The fixed workflow reproduces the spatially analytic reference and the
fixed-rank HSGP comparisons. The adaptive workflow uses the extended rank
ladder and relaxed tolerances recorded for the manuscript result. It starts
from new pilot fits because large intermediate optimization checkpoints are
not distributed in this compact package.

See `DATA_PROVENANCE.md` and `PROVENANCE.md` before publishing the folder.
