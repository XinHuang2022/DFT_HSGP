# Direct spectral-multiplier simulation

This folder contains the simulation in which the spatial spectral multiplier
is directly evaluable but the corresponding spatial half-spectrum is computed
by numerical Hankel inversion.  The simulator and numerical-inversion
reference do not use the HSGP basis.

The principal comparison uses the same 20 simulated records, site split, and
optimizer seeds for:

1. the numerical-inversion reference;
2. fixed HSGP at the stabilized internal-reference rank; and
3. frequency-adaptive HSGP selected from a nested basis.

The initial screened rank remains part of the adaptive workflow, but it is not
included as a principal manuscript comparator by the retained-summary helper.

## Entry points

- Interactive: `main_direct_multiplier_simulation.m`
- Batch preflight: `run_direct_multiplier_experiment("preflight", 1)`
- End-to-end smoke run: `run_direct_multiplier_experiment("smoke", 1)`
- Full paper experiment: `run_direct_multiplier_experiment("paper", 5)`

The main script defaults to retained outputs.  The full run is computationally
substantial and writes checkpoints and tables beneath `results/`.

## Folder contents

- `src/` contains the matched-run driver and configuration.
- `src/model/` contains the minimal direct-multiplier, quadrature, HSGP,
  likelihood, simulation, and prediction dependency set.
- `input/` contains the retained 20-replicate summaries and rank profiles used
  for regression checks and immediate inspection.
- `results/` receives new smoke and paper runs.

