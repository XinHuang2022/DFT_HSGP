experimentRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(experimentRoot);
report31 = run_cerra_experiment("preflight_31day", 1);
report90 = run_cerra_experiment("preflight_90day", 1);
assert(all(report31.passed) && all(report90.passed), ...
    'CERRA package preflight did not pass.');
fprintf('CERRA 31-day and 90-day package preflight tests passed.\n');
