function result = run_pacific_hsgp_pipeline(stage, configOrPreset, ...
        useParallel)
%RUN_PACIFIC_HSGP_PIPELINE Dispatch preflight, basis audit, or CV.

    if nargin < 1 || isempty(stage), stage = 'checks'; end
    if nargin < 2 || isempty(configOrPreset)
        configOrPreset = 'local7_cv5';
    end
    if isstruct(configOrPreset)
        cfg = configOrPreset;
    else
        cfg = pacific_hsgp_config(configOrPreset);
    end
    if nargin >= 3 && ~isempty(useParallel)
        cfg.parallel.useParallel = logical(useParallel);
    end
    switch lower(string(stage))
        case "checks"
            result = run_pacific_hsgp_preflight(cfg.preset);
        case "audit"
            result = run_pacific_hsgp_basis_audit(cfg);
        case "cv"
            result = run_pacific_hsgp_cv(cfg, cfg.parallel.useParallel);
        otherwise
            error('Use stage checks, audit, or cv.');
    end
end
