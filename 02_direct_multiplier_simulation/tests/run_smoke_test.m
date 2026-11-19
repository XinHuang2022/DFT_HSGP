function result = run_smoke_test()
%RUN_SMOKE_TEST Execute one small end-to-end direct-multiplier run.

    experimentRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(experimentRoot, 'src'));
    addpath(fullfile(experimentRoot, 'src', 'model'));
    preflight_direct_multiplier(experimentRoot);
    cfg = matched_config('smoke');
    cfg.outputRoot = fullfile(experimentRoot, 'results', 'smoke');
    result = run_matched_replicate(cfg);
    fprintf('Direct-multiplier smoke run passed: %s\n', result.outputDir);
end

