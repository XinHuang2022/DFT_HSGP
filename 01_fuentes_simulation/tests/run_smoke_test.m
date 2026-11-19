function run_smoke_test()
%RUN_SMOKE_TEST Run the Fuentes structural and numerical preflight checks.

    experimentRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(experimentRoot, 'src'));
    addpath(fullfile(experimentRoot, 'src', 'model'));
    report = preflight_fuentes_simulation(fitted_adaptive_config());
    assert(all(report.passed));
    fprintf('Fuentes simulation preflight passed.\n');
end

