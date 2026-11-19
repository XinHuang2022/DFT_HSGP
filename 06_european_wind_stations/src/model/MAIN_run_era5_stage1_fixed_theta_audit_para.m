clear; clc;
cfg = era5_northsea_config();
cfg.parallel.useParallel = true;
run_era5_stage1_fixed_theta_audit(cfg);
