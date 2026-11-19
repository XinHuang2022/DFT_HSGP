experimentRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(experimentRoot, fullfile(experimentRoot, 'src'), ...
    fullfile(experimentRoot, 'src', 'experiment'), ...
    fullfile(experimentRoot, 'src', 'model'), ...
    fullfile(experimentRoot, 'src', 'adaptive'));
report = preflight_european_package(experimentRoot);
assert(all(report.passed));
retained = summarize_retained_european(experimentRoot);
assert(~isempty(retained.fivefoldSummary));
assert(~isempty(retained.relaxedSummary));
fprintf('European station preflight and retained-result checks passed.\n');
