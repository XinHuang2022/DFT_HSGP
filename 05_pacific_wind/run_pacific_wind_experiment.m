function result = run_pacific_wind_experiment(mode, numberOfWorkers)
%RUN_PACIFIC_WIND_EXPERIMENT Batch entry point for the Pacific application.

    if nargin < 1 || isempty(mode), mode = "preflight"; end
    if nargin < 2 || isempty(numberOfWorkers), numberOfWorkers = 5; end

    experimentRoot = fileparts(mfilename('fullpath'));
    sourceRoot = fullfile(experimentRoot, 'src');
    modelRoot = fullfile(sourceRoot, 'model');
    adaptiveRoot = fullfile(sourceRoot, 'adaptive');
    addpath(sourceRoot, modelRoot, adaptiveRoot);
    pathCleanup = onCleanup(@() rmpath(sourceRoot, modelRoot, adaptiveRoot));

    mode = lower(string(mode));
    if mode == "retained"
        result = summarize_retained_pacific(experimentRoot);
        return;
    end

    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    outputDir = fullfile(experimentRoot, 'results', ...
        sprintf('%s_%s', char(mode), timestamp));
    if ~isfolder(outputDir), mkdir(outputDir); end
    previousOutput = getenv('PACIFIC_REPRO_OUTPUT_DIR');
    environmentCleanup = onCleanup(@() ...
        setenv('PACIFIC_REPRO_OUTPUT_DIR', previousOutput));
    setenv('PACIFIC_REPRO_OUTPUT_DIR', outputDir);

    switch mode
        case "preflight"
            result = preflight_pacific_package(experimentRoot);
        case "fixed_paper"
            preflight_pacific_package(experimentRoot);
            cfg = pacific_hsgp_config('full17_cv5');
            cfg.parallel.numWorkers = numberOfWorkers;
            run_pacific_hsgp_basis_audit(cfg);
            result = run_pacific_hsgp_cv(cfg, true);
            result.outputDir = outputDir;
        case "adaptive_paper"
            preflight_pacific_package(experimentRoot);
            override = struct();
            override.nWorkers = numberOfWorkers;
            override.outputRoot = fullfile(outputDir, 'adaptive_rank');
            override.referenceJ = [64 72 80 88 96];
            override.candidateJ = 16:4:92;
            override.epsilonReferencePerSite = 0.03;
            override.epsilonKlPerSite = 0.05;
            override.reuseTrustedPilot = false;
            result = run_adaptive_realdata('pacific', override);
            result.outputDir = outputDir;
        otherwise
            error('Unknown Pacific execution mode: %s', mode);
    end
end
