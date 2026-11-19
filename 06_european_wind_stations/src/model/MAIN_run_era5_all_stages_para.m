clear;clc;cfg=era5_northsea_config();cfg.parallel.useParallel=true;
run_era5_preflight();
run_era5_stage1_fixed_theta_audit(cfg);
run_era5_stage2_hsgp_refit_screen(cfg);
run_era5_stage3_paired_cv(cfg);
