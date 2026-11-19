function result = run_fuentes_experiment(mode, numberOfWorkers)
%RUN_FUENTES_EXPERIMENT Batch entry point for the Fuentes simulation.

    if nargin < 1 || isempty(mode)
        mode = "preflight";
    end
    if nargin < 2 || isempty(numberOfWorkers)
        numberOfWorkers = 1;
    end

    experimentRoot = fileparts(mfilename('fullpath'));
    addpath(fullfile(experimentRoot, 'src'));
    addpath(fullfile(experimentRoot, 'src', 'model'));
    cfg = fitted_adaptive_config();
    cfg.nWorkers = numberOfWorkers;
    cfg.useParallel = numberOfWorkers > 1;

    switch lower(string(mode))
        case "retained"
            result = summarize_retained_fuentes(experimentRoot);
        case "preflight"
            result = preflight_fuentes_simulation(cfg);
        case "paper"
            preflight_fuentes_simulation(cfg);
            result = run_fitted_adaptive_fuentes(cfg);
        otherwise
            error('Unknown Fuentes execution mode: %s', mode);
    end
end

