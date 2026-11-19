function result = run_all(preset, useParallel)
%RUN_ALL Main entry point for the Pacific wind pilot suite.
%
%   result = run_all('audit');
%   result = run_all('smoke3');
%   result = run_all('pilot5', true);
%   result = run_all('local5_cv5', true);
%   result = run_all('local7_cv5', true);

    if nargin < 1 || isempty(preset)
        preset = 'audit';
    end
    if strcmpi(preset, 'audit')
        result = run_pacific_audit('audit');
        return;
    end
    run_pacific_unit_checks(preset);
    if nargin < 2
        result = run_pacific_cv(preset);
    else
        result = run_pacific_cv(preset, useParallel);
    end
end
