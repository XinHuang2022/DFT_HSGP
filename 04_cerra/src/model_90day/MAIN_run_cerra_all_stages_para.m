clear;clc;cfg=cerra_northsea_config();cfg.parallel.useParallel=true;
run_cerra_dense_preflight(cfg);
run_cerra_dense_stage1_fixed_theta_screen(cfg);
run_cerra_dense_stage2_selective_hsgp_refit(cfg);
run_cerra_dense_stage3_exact_baselines(cfg);

