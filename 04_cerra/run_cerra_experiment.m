function result = run_cerra_experiment(mode, numberOfWorkers)
%RUN_CERRA_EXPERIMENT Batch entry point for the dense CERRA application.

    if nargin < 1 || isempty(mode), mode = "preflight_90day"; end
    if nargin < 2 || isempty(numberOfWorkers), numberOfWorkers = 1; end

    experimentRoot = fileparts(mfilename('fullpath'));
    sourceRoot = fullfile(experimentRoot, 'src');
    addpath(sourceRoot);
    sourceCleanup = onCleanup(@() rmpath(sourceRoot));
    mode = lower(string(mode));

    if mode == "retained"
        result = summarize_retained_cerra(experimentRoot);
        return;
    end

    if contains(mode, "31day")
        modelRoot = fullfile(sourceRoot, 'model_31day');
        windowLabel = "31day";
    elseif contains(mode, "90day")
        modelRoot = fullfile(sourceRoot, 'model_90day');
        windowLabel = "90day";
    else
        error('CERRA mode must identify either the 31-day or 90-day window.');
    end
    addpath(modelRoot);
    modelCleanup = onCleanup(@() rmpath(modelRoot));

    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    outputDir = fullfile(experimentRoot, 'results', ...
        sprintf('%s_%s', char(mode), timestamp));
    if ~isfolder(outputDir), mkdir(outputDir); end

    previousOutput = getenv('CERRA_REPRO_OUTPUT_DIR');
    previousWorkers = getenv('CERRA_NUM_WORKERS');
    environmentCleanup = onCleanup(@() restore_environment( ...
        previousOutput, previousWorkers));
    setenv('CERRA_REPRO_OUTPUT_DIR', outputDir);
    setenv('CERRA_NUM_WORKERS', sprintf('%d', numberOfWorkers));

    cfg = cerra_northsea_config();
    cfg.parallel.numWorkers = numberOfWorkers;
    if startsWith(mode, "preflight")
        result = preflight_cerra_package(experimentRoot, cfg, windowLabel);
    elseif startsWith(mode, "paper")
        preflight_cerra_package(experimentRoot, cfg, windowLabel);
        run_cerra_dense_stage1_fixed_theta_screen(cfg);
        run_cerra_dense_stage2_selective_hsgp_refit(cfg);
        run_cerra_dense_stage3_exact_baselines(cfg);
        result = struct('outputDir', outputDir, ...
            'workflow', "CERRA " + windowLabel);
    else
        error('Unknown CERRA execution mode: %s', mode);
    end
end

function restore_environment(outputValue, workerValue)
    setenv('CERRA_REPRO_OUTPUT_DIR', outputValue);
    setenv('CERRA_NUM_WORKERS', workerValue);
end
