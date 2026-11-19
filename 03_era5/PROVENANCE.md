# ERA5 code and result provenance

The fixed-rank and analytic source was copied from
`deliverables/era5_northsea_eq9_suite_v0_1_0`, which is byte-identical to the
trusted implementation under `wind_field_implement/Comp/ERA5`. The adaptive
source was copied from `deliverables/realdata_adaptive_hsgp_v0_1_0`.

The preserved fixed-rank tables come from the trusted
`stage1_fixed_theta_J_audit`, `stage2_hsgp_refit_screen`, and
`stage3_paired_cv` outputs. The preserved adaptive tables come from the
completed `Results_extended_relaxed/era5` run, which is the source of the
reported mean held-site RMSE near 0.0714.

Only path handling and release entry points were changed. The covariance
model, temporal DFT, Whittle objective, coordinate scaling, mean design,
folds, parameter bounds, rank ladders, and tolerances were not redefined.
