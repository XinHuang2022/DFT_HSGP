clear; clc;
suiteRoot = fileparts(mfilename('fullpath'));
addpath(suiteRoot);
run_adaptive_realdata('era5');

