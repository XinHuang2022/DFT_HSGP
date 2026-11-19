# Fuentes analytic-reference simulation

This folder contains the 20-replicate experiment used to validate
frequency-specific HSGP ranks against the analytic Fuentes covariance.  The
model has an analytic spatial half-spectrum, so this experiment primarily
tests approximation fidelity and rank allocation rather than avoidance of
numerical spatial inversion.

The held-out sites are excluded from parameter fitting and rank selection.
The adaptive profiles are selected at the replicate-specific analytic-model
estimate and compared with the analytic covariance and a common HSGP rank.

## Entry points

- Interactive: `main_fuentes_simulation.m`
- Batch preflight: `run_fuentes_experiment("preflight", 1)`
- Full paper experiment: `run_fuentes_experiment("paper", 5)`
- Structural test: `run('tests/run_smoke_test.m')`

The sectioned main script defaults to retained outputs.  Change
`executionMode` to `"paper"` only when a complete rerun is intended.

## Folder contents

- `src/` contains the experiment driver and reporting functions.
- `src/model/` contains the minimal validated Fuentes/HSGP dependency set.
- `input/` contains the fixed parameter estimates, prediction regression
  targets, compact reference metadata, and retained result summaries.
- `results/` receives timestamped reruns.

The compact MAT file in `input/` contains only the geometry, configuration,
fold manifest, source configuration, and oracle rank profile required by this
workflow.  It replaces a much larger historical MAT file without changing the
experiment.

