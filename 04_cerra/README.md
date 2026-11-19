# CERRA North Sea reanalysis

This folder packages the dense-grid CERRA experiment used in the manuscript.
It contains the extracted native marine grid for the North Sea, the fixed
300-site evaluation set, and the fixed-rank HSGP workflows for both temporal
windows reported in the paper.

## Entry points

- Interactive: `main_cerra_reanalysis.m`
- 31-day preflight: `run_cerra_experiment("preflight_31day", 1)`
- 90-day preflight: `run_cerra_experiment("preflight_90day", 1)`
- 31-day paper workflow: `run_cerra_experiment("paper_31day", 6)`
- 90-day paper workflow: `run_cerra_experiment("paper_90day", 6)`

The interactive script defaults to `"retained"`. New runs use timestamped
directories beneath `results/` and do not overwrite the preserved summaries.

## Data contract

The local dataset contains 360 six-hourly observations at 3,701 marine grid
locations. The 31-day configuration selects the 124 January 2021 observations,
whereas the 90-day configuration uses the complete December 2020--February
2021 record. Both configurations reconstruct the same deterministic 300-site
evaluation set and fit the remaining 3,401 sites.

The rank ladder is `M = 250, 500, 750, 1000`, with radial ordering of the
two-dimensional Laplacian modes. Raw CERRA downloads and the acquisition
environment are excluded because the extracted matrix and manifests are
sufficient for the MATLAB experiment.

See `DATA_PROVENANCE.md` and `PROVENANCE.md` before publishing the folder.
