# Pacific wind application

This folder packages the 480-time-point Pacific wind experiment on the fixed
17-by-17 grid. It includes the extracted data used by the analysis, the
deterministic five whole-site folds, the fixed-rank HSGP workflow, and the
frequency-adaptive HSGP workflow.

## Entry points

- Interactive: `main_pacific_wind.m`
- Batch preflight: `run_pacific_wind_experiment("preflight", 1)`
- Fixed-rank paper workflow: `run_pacific_wind_experiment("fixed_paper", 5)`
- Frequency-adaptive workflow:
  `run_pacific_wind_experiment("adaptive_paper", 5)`

The interactive script defaults to `"retained"`, so opening it does not start
an optimization. New runs are written to timestamped folders under `results/`.

## Data and fold contract

The wind matrix has 480 six-hourly records at 289 sites. The coordinate file
defines a complete 17-by-17 grid. Both workflows reuse the preserved five-fold
whole-site partition. The preflight check verifies the local matrix, station
order, fold assignment, and model interfaces before a long run begins.

The adaptive paper profile uses the rank ladder and tolerances retained for
the manuscript experiment. Large fit checkpoints and prediction matrices are
not duplicated in this compact package.
