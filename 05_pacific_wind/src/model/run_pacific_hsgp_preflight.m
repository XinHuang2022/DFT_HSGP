function result = run_pacific_hsgp_preflight(preset)
%RUN_PACIFIC_HSGP_PREFLIGHT Base and HSGP checks before long jobs.

    if nargin < 1 || isempty(preset), preset = 'local7_cv5'; end
    result.base = run_pacific_unit_checks(preset);
    result.hsgp = run_pacific_hsgp_unit_checks(preset);
    fprintf('Pacific HSGP preflight passed: %d\n', ...
        all(result.base.passed) && all(result.hsgp.passed));
end
