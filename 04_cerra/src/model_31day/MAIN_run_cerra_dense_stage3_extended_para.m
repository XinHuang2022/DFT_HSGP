clear;clc;cfg=cerra_northsea_config();cfg.parallel.useParallel=true;
cfg.stage3.exactSizes=cfg.stage3.extendedExactSizes;
run_cerra_dense_stage3_exact_baselines(cfg);
