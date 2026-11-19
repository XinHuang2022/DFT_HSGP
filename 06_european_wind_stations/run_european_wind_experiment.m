function result = run_european_wind_experiment(mode, numberOfWorkers)
%RUN_EUROPEAN_WIND_EXPERIMENT Batch entry point for European stations.

    if nargin < 1 || isempty(mode), mode = "preflight"; end
    if nargin < 2 || isempty(numberOfWorkers), numberOfWorkers = 5; end

    experimentRoot = fileparts(mfilename('fullpath'));
    sourceRoot = fullfile(experimentRoot, 'src');
    experimentSource = fullfile(sourceRoot, 'experiment');
    modelRoot = fullfile(sourceRoot, 'model');
    adaptiveRoot = fullfile(sourceRoot, 'adaptive');
    addpath(sourceRoot, experimentSource, modelRoot, adaptiveRoot);
    pathCleanup = onCleanup(@() rmpath(sourceRoot, experimentSource, ...
        modelRoot, adaptiveRoot));

    mode = lower(string(mode));
    if mode == "retained"
        result = summarize_retained_european(experimentRoot);
        return;
    end

    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    outputDir = fullfile(experimentRoot, 'results', ...
        sprintf('%s_%s', char(mode), timestamp));
    if ~isfolder(outputDir), mkdir(outputDir); end
    previousOutput = getenv('EUROPEAN_WIND_REPRO_OUTPUT_DIR');
    environmentCleanup = onCleanup(@() ...
        setenv('EUROPEAN_WIND_REPRO_OUTPUT_DIR', previousOutput));
    setenv('EUROPEAN_WIND_REPRO_OUTPUT_DIR', outputDir);

    switch mode
        case "preflight"
            result = preflight_european_package(experimentRoot);
        case {"fixed_paper", "algorithm1", "algorithm2", "paper"}
            preflight_european_package(experimentRoot);
            result = run_european_fivefold(mode, numberOfWorkers);
            result.outputDir = outputDir;
        otherwise
            error('Unknown European station execution mode: %s', mode);
    end
end
