function run_preflight_test()
%RUN_PREFLIGHT_TEST Check the direct-multiplier release structure.

    experimentRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(experimentRoot, 'src'));
    addpath(fullfile(experimentRoot, 'src', 'model'));
    report = preflight_direct_multiplier(experimentRoot);
    assert(all(report.passed));
    fprintf('Direct-multiplier preflight passed.\n');
end

