# Fuentes source provenance

The release files were copied without mathematical changes from the following
validated implementations:

- experiment driver and configuration:
  `deliverables/fuentes_fitted_adaptive_hsgp_v0_3_0`;
- Gaussian KL calculation:
  `deliverables/fuentes_frequency_adaptive_hsgp_v0_1_0/src`;
- Fuentes model, simulation, analytic covariance, and HSGP routines:
  `wind_field_implement/Synthetic/eq9_half_spectral_syn_suite_v0_3_7_hsgp_recovery_program`;
- retained fitted-parameter results:
  `deliverables/fuentes_fitted_adaptive_hsgp_v0_3_0/results/fitted_adaptive_20261022_181153`.

Only path construction was changed so the experiment is local to this folder.
The compact reference MAT file was created by retaining the five fields read by
the experiment: `sourceCfg`, `sites`, `configurations`, `foldManifest`, and
`klRankProfiles`.

