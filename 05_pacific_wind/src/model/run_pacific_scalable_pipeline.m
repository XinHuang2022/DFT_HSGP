function result = run_pacific_scalable_pipeline(stage, ...
        configOrPreset, useParallel)
%RUN_PACIFIC_SCALABLE_PIPELINE Dispatch the staged scalable workflow.
%
% Examples:
%   run_pacific_scalable_pipeline('checks', 'local7_cv5', false)
%   run_pacific_scalable_pipeline('nugget', 'full17_cv5', true)
%   run_pacific_scalable_pipeline('subset', 'full17_cv5', true)
%   run_pacific_scalable_pipeline('rff0', 'full17_cv5', true)
%   run_pacific_scalable_pipeline('rfffit', 'full17_cv5', true)

    if nargin < 1 || isempty(stage)
        stage = 'checks';
    end
    if nargin < 2 || isempty(configOrPreset)
        configOrPreset = 'full17_cv5';
    end
    if isstruct(configOrPreset)
        cfg = configOrPreset;
    else
        cfg = pacific_scalable_config(configOrPreset);
    end
    if nargin >= 3 && ~isempty(useParallel)
        cfg.parallel.useParallel = logical(useParallel);
    end
    stage = lower(string(stage));
    switch stage
        case "checks"
            result.base = run_pacific_unit_checks(cfg.preset);
            result.scalable = run_pacific_scalable_unit_checks(cfg.preset);
        case "nugget"
            result = run_pacific_nugget_study(cfg, ...
                cfg.parallel.useParallel);
        case "subset"
            result = run_pacific_subset_plugin_cv(cfg, ...
                cfg.parallel.useParallel);
        case "rff0"
            result = run_pacific_rff_stage0(cfg, ...
                cfg.parallel.useParallel);
        case "rfffit"
            result = run_pacific_rff_fit_cv(cfg, ...
                cfg.parallel.useParallel);
        otherwise
            error(['Unknown scalable stage %s. Use checks, nugget, ', ...
                'subset, rff0, or rfffit.'], char(stage));
    end
end
