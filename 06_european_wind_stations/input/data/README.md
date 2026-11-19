# Analysis data

`ghcnh_complete_panel_expanded.mat` contains the 360-by-839 analysis matrix,
coordinates, ordered site table, and target times. The CSV station manifest
records the extraction diagnostics, and `cv_fold_manifest.csv` protects the
five whole-site fold assignment.

This release uses only stations with complete matches at all 360 target times.
Incomplete records were screened during data preparation rather than imputed
inside the statistical model. Raw downloads and the screening pipeline are not
required to rerun the packaged inference experiment.
