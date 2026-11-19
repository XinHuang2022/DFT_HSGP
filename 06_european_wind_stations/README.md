# European wind stations

This folder packages the European station experiment based on the extracted
GHCNh panel used in the manuscript. It contains 360 six-hourly wind-speed
records at 839 stations, the deterministic five whole-site folds, the
fixed-rank comparator, Algorithm 1, and Algorithm 2.

## Entry points

- Interactive: `main_european_wind_stations.m`
- Batch preflight: `run_european_wind_experiment("preflight", 1)`
- Fixed-rank comparator: `run_european_wind_experiment("fixed_paper", 5)`
- Algorithm 1: `run_european_wind_experiment("algorithm1", 5)`
- Algorithm 2: `run_european_wind_experiment("algorithm2", 5)`
- Complete comparison: `run_european_wind_experiment("paper", 5)`

The full modes use five outer workers. Each worker processes one spatial fold,
while frequency-level calculations inside a fold remain serial. This avoids
nested pools and matches the computational structure of the retained run.
Because the adaptive procedures reuse fold-specific pilot fits, the complete
workflow runs the fixed-rank folds before Algorithms 1 and 2.

The interactive script defaults to `"retained"`. New fits are written under
`results/` and do not overwrite the preserved summaries.

## Data contract

The included analysis-ready panel contains 839 stations with complete records
at 360 target times from June through August 2024. The station manifest fixes
the column order. Five deterministic, spatially distributed whole-site folds
contain 167 or 168 held-out stations each. Their exact assignment is stored in
`input/data/cv_fold_manifest.csv` and checked before fitting.
