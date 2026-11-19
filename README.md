# Temporal Fourier likelihoods with spatial HSGP approximations

This repository contains the MATLAB implementations and analysis-ready data
used in the manuscript *Temporal Fourier Likelihoods with Spatial
Hilbert-Space Gaussian Process Approximations*. The code studies spatial
reconstruction for stationary space–time Gaussian processes observed at many
locations on a common regular time grid.

## Scientific motivation

A temporal discrete Fourier transform (DFT) converts a regularly sampled
space--time record into frequency-indexed vectors of spatial coefficients.
Under a Whittle likelihood, these vectors are treated as independent across
temporal frequencies, leaving one spatial covariance calculation at each
frequency.

For a model specified through a joint spectral density, the required spatial
covariance is the inverse spatial Fourier transform of the sampled spectral
multiplier. This transform may be unavailable in closed form and expensive to
evaluate repeatedly by numerical quadrature. The implementation in this
repository instead constructs a Hilbert-space Gaussian process (HSGP)
approximation by evaluating the spatial spectral multiplier at the
eigenfrequencies of a Laplacian basis on an extended domain.

The resulting workflow combines:

1. a temporal DFT of the observations;
2. a frequency-wise Whittle likelihood;
3. a shared nested Laplacian basis for the spatial calculations;
4. HSGP weights obtained directly from the sampled spectral multiplier; and
5. inverse-DFT reconstruction at spatial locations excluded from fitting.

The HSGP step replaces repeated spatial Fourier inversion with finite-basis
calculations. Its accuracy still depends on the domain extension and basis
truncation, which are assessed separately in the numerical experiments.

## Frequency-adaptive spatial rank

The spatial resolution needed to approximate the covariance can vary across
temporal frequencies. Retaining the same large basis at every frequency may
therefore perform unnecessary calculations. The manuscript develops two
procedures that select a rank separately for each frequency from a common
nested basis.

- **Algorithm 1** compares candidate covariances with a stabilized high-rank
  HSGP reference using a normalized Kullback–Leibler criterion.
- **Algorithm 2** uses a covariance-scaled trace certificate to screen
  candidate ranks. This provides a more conservative but potentially cheaper
  rank assessment when the number of sites or the candidate ladder is large.

In both procedures, the rank profile is selected from the fitted covariance
model and reassessed after covariance-parameter estimation. Held-out spatial
responses are not used to select the ranks. Analytic or numerical-inversion
covariances are used only in experiments where an external assessment is
available.

## Numerical studies

The repository is organized by dataset. Each folder is self-contained and has
its own data contract, MATLAB entry point, configuration, retained summaries,
and provenance record.

| Folder | Study | Role in the manuscript |
|---|---|---|
| `01_fuentes_simulation` | Fuentes spectral model | Evaluates HSGP approximation and adaptive rank selection against an analytic full-rank covariance. |
| `02_direct_multiplier_simulation` | Directly specified spectral multiplier | Compares numerical spatial inversion, stabilized fixed-rank HSGP, and frequency-adaptive HSGP when the half-spectrum is evaluated numerically. |
| `03_era5` | North Sea ERA5 reanalysis | Studies held-site reconstruction on a regular 17-by-17 marine grid. |
| `04_cerra` | North Sea CERRA reanalysis | Examines fixed-rank reconstruction on a denser 3,701-site grid using 31-day and 90-day records. |
| `05_pacific_wind` | Pacific wind field | Compares fixed and frequency-adaptive HSGP under five whole-site folds. |
| `06_european_wind_stations` | European GHCNh stations | Compares fixed rank, Algorithm 1, and Algorithm 2 at 839 irregularly located stations. |

The simulation studies distinguish approximation error from parameter
estimation and reconstruction error. The environmental applications assess
the methods on locations excluded for their complete temporal records, so the
reported errors measure spatial reconstruction rather than temporal
forecasting.

## Repository structure

Each dataset folder follows the same pattern:

```text
dataset/
|-- main_*.m                 sectioned interactive entry point
|-- run_*_experiment.m      batch entry point
|-- input/
|   |-- data/               analysis-ready data and protected manifests
|   `-- retained/           compact reference summaries
|-- src/                    model and experiment functions
|-- tests/                  lightweight preflight checks
|-- results/                timestamped outputs from new runs
|-- README.md               dataset-specific instructions
|-- DATA_PROVENANCE.md      data origin and redistribution notes
|-- PROVENANCE.md           source-code and retained-result origin
`-- SOURCE_MANIFEST.csv     file sizes and SHA-256 hashes
```

The retained outputs permit immediate inspection and regression checking.
They do not replace a complete rerun of an experiment.

## Requirements

The code has been prepared and preflight-tested with:

- MATLAB R2023a Update 2;
- Parallel Computing Toolbox for multi-worker experiments; and
- sufficient memory for the selected dataset and rank ladder.

The complete fitting workflows can be computationally substantial. The
default interactive mode reads retained summaries and never starts a long
optimization automatically.

## Quick start

Clone or download the repository, start MATLAB in the repository root, and
open one dataset folder. For example:

```matlab
cd('01_fuentes_simulation')
main_fuentes_simulation
```

The first section of every `main_*.m` file defines the execution mode. Common
modes are:

- `"retained"`: inspect preserved results without fitting;
- `"preflight"`: validate data, protected splits, and model interfaces;
- `"smoke"`: run a small end-to-end calculation where available; and
- `"paper"`: run the complete experiment.

For a batch preflight, call the dataset function directly. For example:

```matlab
cd('03_era5')
report = run_era5_experiment("preflight", 1);
```

The European station experiment uses fold-level parallelism:

```matlab
cd('06_european_wind_stations')
result = run_european_wind_experiment("paper", 5);
```

With five workers, the five spatial cross-validation folds are processed in
parallel. Frequency-level calculations remain serial within each worker to
avoid nested parallel pools.

## Reproducibility safeguards

The packaged workflows preserve the statistical and computational choices
used in the corresponding experiments:

- fixed fitting and reconstruction locations;
- protected cross-validation manifests;
- coordinate transformations and temporal sampling conventions;
- parameter bounds and optimizer settings;
- HSGP domain extensions and nested rank ladders;
- adaptive-rank tolerances and stopping rules; and
- retained summaries for regression checking.

Preflight functions check dimensions, finite values, site ordering, fold
assignments, and core model interfaces before a long run begins. New results
are written to timestamped directories and do not overwrite the retained
reference outputs.

## Data and provenance

Only analysis-ready datasets needed by the MATLAB experiments are included.
Large raw ERA5, CERRA, and station-download payloads are not duplicated. Each
dataset folder documents the origin of its data and code in
`DATA_PROVENANCE.md` and `PROVENANCE.md`.

Before a public release, the applicable data-source attribution and
redistribution conditions should be confirmed. A software licence and final
citation metadata should also be added when the manuscript and archived
repository release are available.

## Citation

If this code contributes to published work, please cite the manuscript
*Temporal Fourier Likelihoods with Spatial Hilbert-Space Gaussian Process
Approximations*. Formal BibTeX metadata will be added after the manuscript or
archived repository release receives a persistent identifier.

