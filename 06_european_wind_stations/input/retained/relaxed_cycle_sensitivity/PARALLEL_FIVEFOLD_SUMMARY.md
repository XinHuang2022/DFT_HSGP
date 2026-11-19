# Relaxed-cycle European GHCNh experiment

## Protected design

- 839 stations and 360 time points;
- deterministic five-fold whole-site partition;
- folds 2--5 executed concurrently with four MATLAB workers;
- unchanged Fuentes multiplier, mean model, rank ladder, `J_ref=48`,
  optimizer settings, seeds, and shared pilot checkpoints;
- held-site responses used only for final evaluation;
- 5,000-resample spatial-block bootstrap with 70 blocks of width 400 km.

## Relaxed procedures

- Algorithm 1: at most four adaptive fits and no more than 1% (one of 181)
  residual rank requests;
- Algorithm 2: at most three adaptive fits and no more than 5% (nine of 181)
  trace-certificate requests.

## Five-fold results

| Procedure | RMSE | MAE | NLPD | folds 2--5 parallel wall time |
|---|---:|---:|---:|---:|
| Fixed `J=48` | 1.9473 [1.8124, 2.1138] | 1.3925 [1.3127, 1.4920] | 2.0720 [2.0036, 2.1565] | 88 min |
| Algorithm 1, four-fit cap | 2.0445 [1.8869, 2.2437] | 1.4552 [1.3653, 1.5685] | 2.1138 [2.0402, 2.2073] | 60.7 min |
| Algorithm 2, three-fit cap | 1.9946 [1.8414, 2.1879] | 1.4180 [1.3346, 1.5232] | 2.0869 [2.0177, 2.1733] | 65.5 min |

Algorithm 2 stopped after two fits in folds 2--3 and three fits in folds 1,
4, and 5.  The validation-only exact audit found no tolerance violation; the
largest normalized KL discrepancy was 0.004493.

Algorithm 1 used four fits in every fold.  Folds 3--5 reached the cap before
the one-frequency stopping rule.  Their reconstructions were therefore
recomputed under the ranks actually used in the fourth parameter fit.  The
reported Algorithm-1 row is a capped sensitivity, not a fully stabilized
execution of the exact-KL procedure.

## Anomaly

The large `complete_task.mat` checkpoint failed after CSV outputs were written
in the relaxed parallel runs.  Metrics, predictions, rank profiles, parameters,
and timing CSV files are complete.  The zero-byte MAT files are invalid and
must not be used.
