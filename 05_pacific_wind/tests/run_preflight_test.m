experimentRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(experimentRoot, fullfile(experimentRoot, 'src'), ...
    fullfile(experimentRoot, 'src', 'model'), ...
    fullfile(experimentRoot, 'src', 'adaptive'));
report = preflight_pacific_package(experimentRoot);
assert(all(report.passed));
retained = summarize_retained_pacific(experimentRoot);
assert(~isempty(retained.fixedMethodSummary));
assert(~isempty(retained.adaptiveMethodSummary));
fprintf('Pacific package preflight and retained-result checks passed.\n');
