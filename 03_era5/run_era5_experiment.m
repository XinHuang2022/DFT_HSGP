function result = run_era5_experiment(mode, numberOfWorkers)
%RUN_ERA5_EXPERIMENT Batch entry point for the ERA5 North Sea application.

    if nargin < 1 || isempty(mode), mode = "preflight"; end
    if nargin < 2 || isempty(numberOfWorkers), numberOfWorkers = 1; end

    experimentRoot = fileparts(mfilename('fullpath'));
    sourceRoot = fullfile(experimentRoot, 'src');
    modelRoot = fullfile(sourceRoot, 'model');
    adaptiveRoot = fullfile(sourceRoot, 'adaptive');
    addpath(sourceRoot, modelRoot, adaptiveRoot);
    pathCleanup = onCleanup(@() rmpath(sourceRoot, modelRoot, adaptiveRoot));

    mode = lower(string(mode));
    if mode == "retained"
        result = summarize_retained_era5(experimentRoot);
        return;
    end

    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    outputDir = fullfile(experimentRoot, 'results', ...
        sprintf('%s_%s', char(mode), timestamp));
    if ~isfolder(outputDir), mkdir(outputDir); end

    previousOutput = getenv('ERA5_REPRO_OUTPUT_DIR');
    previousWorkers = getenv('ERA5_NUM_WORKERS');
    environmentCleanup = onCleanup(@() restore_environment( ...
        previousOutput, previousWorkers));
    setenv('ERA5_REPRO_OUTPUT_DIR', outputDir);
    setenv('ERA5_NUM_WORKERS', sprintf('%d', numberOfWorkers));

    switch mode
        case "preflight"
            result = preflight_era5_package(experimentRoot);
        case "fixed_paper"
            preflight_era5_package(experimentRoot);
            cfg = era5_northsea_config();
            cfg.parallel.numWorkers = numberOfWorkers;
            run_era5_stage1_fixed_theta_audit(cfg);
            run_era5_stage2_hsgp_refit_screen(cfg);
            run_era5_stage3_paired_cv(cfg);
            result = struct('outputDir', outputDir, ...
                'workflow', "analytic and fixed-rank ERA5");
        case "adaptive_paper"
            preflight_era5_package(experimentRoot);
            override = struct();
            override.nWorkers = numberOfWorkers;
            override.outputRoot = outputDir;
            override.referenceJ = [64 72 80 88 96];
            override.candidateJ = 16:4:92;
            override.epsilonReferencePerSite = 0.03;
            override.epsilonKlPerSite = 0.05;
            override.reuseTrustedPilot = false;
            result = run_adaptive_realdata('era5', override);
            result.outputDir = outputDir;
        otherwise
            error('Unknown ERA5 execution mode: %s', mode);
    end
end

function restore_environment(outputValue, workerValue)
    setenv('ERA5_REPRO_OUTPUT_DIR', outputValue);
    setenv('ERA5_NUM_WORKERS', workerValue);
end
