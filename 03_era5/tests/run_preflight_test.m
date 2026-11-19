experimentRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(experimentRoot);
result = run_era5_experiment("preflight", 1);
assert(all(result.passed), 'ERA5 package preflight did not pass.');
fprintf('ERA5 package preflight test passed.\n');
