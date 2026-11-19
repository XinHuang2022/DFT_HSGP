function result = run_direct_multiplier_experiment(mode, numberOfWorkers, ...
        auditRule, klTolerance)
%RUN_DIRECT_MULTIPLIER_EXPERIMENT Batch entry point for the second simulation.

    if nargin < 1 || isempty(mode), mode = "preflight"; end
    if nargin < 2 || isempty(numberOfWorkers), numberOfWorkers = 1; end
    if nargin < 3 || isempty(auditRule), auditRule = "exact_kl"; end
    if nargin < 4 || isempty(klTolerance), klTolerance = 0.01; end

    experimentRoot = fileparts(mfilename('fullpath'));
    addpath(fullfile(experimentRoot, 'src'));
    addpath(fullfile(experimentRoot, 'src', 'model'));

    switch lower(string(mode))
        case "retained"
            result = summarize_retained_direct_multiplier(experimentRoot);
        case "preflight"
            result = preflight_direct_multiplier(experimentRoot);
        case "smoke"
            preflight_direct_multiplier(experimentRoot);
            cfg = matched_config('smoke');
            cfg.auditRule = char(auditRule);
            cfg.klTolerance = klTolerance;
            result = run_matched_replicate(cfg);
        case "paper"
            preflight_direct_multiplier(experimentRoot);
            result.outputDir = run_matched_20(numberOfWorkers, ...
                auditRule, klTolerance);
        otherwise
            error('Unknown direct-multiplier execution mode: %s', mode);
    end
end

