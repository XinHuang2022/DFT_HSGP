clear; clc;
suiteRoot = fileparts(mfilename('fullpath'));
addpath(suiteRoot);
cfg = pacific_hsgp_config('full17_cv5');
cfg.parallel.useParallel = true;
cfg.parallel.numWorkers = 5;
run_pacific_hsgp_pipeline('cv', cfg, true);
