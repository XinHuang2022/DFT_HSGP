function result = run_pacific_subset_plugin_cv(configOrPreset, useParallel)
%RUN_PACIFIC_SUBSET_PLUGIN_CV Fit theta on subsets, condition on all sites.

    if nargin < 1 || isempty(configOrPreset)
        cfg = pacific_scalable_config('full17_cv5');
    elseif isstruct(configOrPreset)
        cfg = configOrPreset;
    else
        cfg = pacific_scalable_config(configOrPreset);
    end
    if nargin >= 2 && ~isempty(useParallel)
        cfg.parallel.useParallel = logical(useParallel);
    end
    ensure_directory(cfg.paths.subsetPluginDir);
    nuggetRule = pacific_resolve_nugget_rule(cfg);
    contract = load_pacific_wind_data(cfg);
    if any(~isfinite(contract.U), 'all')
        error('Subset plug-in CV requires a complete Pacific panel.');
    end
    subset = make_pacific_subset(contract, cfg.subset.name);
    [folds, foldManifest, foldSummary] = make_pacific_cv_folds( ...
        subset, cfg.cv.holdoutPattern);
    nFold = numel(folds);
    tasks = cell(nFold, 1);
    timerValue = tic;
    if cfg.parallel.useParallel && parallel_available()
        pool = gcp('nocreate');
        if isempty(pool) && cfg.parallel.autoStartPool
            pool = parpool(min(cfg.parallel.numWorkers, nFold));
        end
        parfor f = 1:nFold
            tasks{f} = run_pacific_subset_plugin_fold(cfg, subset, ...
                folds(f).targetIndex, folds(f).foldID, nuggetRule);
        end
    else
        for f = 1:nFold
            fprintf('Subset plug-in fold %d/%d.\n', f, nFold);
            tasks{f} = run_pacific_subset_plugin_fold(cfg, subset, ...
                folds(f).targetIndex, folds(f).foldID, nuggetRule);
        end
    end
    elapsedSeconds = toc(timerValue);
    metrics = concatenate(tasks, 'metrics');
    predictions = concatenate(tasks, 'predictions');
    fitDiagnostics = concatenate(tasks, 'fitDiagnostics');
    subsetSensitivity = concatenate(tasks, 'subsetSensitivity');
    parameterManifest = concatenate(tasks, 'parameterManifest');
    checks = concatenate(tasks, 'checks');
    timing = concatenate(tasks, 'timing');
    rffStage0Inputs = make_rff_stage0_inputs(tasks, subset);
    methodSummary = summarize_metrics(metrics);
    runSummary = table(string(cfg.version), string(cfg.preset), ...
        height(subset.siteTable), nFold, ...
        cfg.scalable.parameterSubsetSize, ...
        cfg.scalable.nSubsetSeedsForSensitivity, nuggetRule.mode, ...
        nuggetRule.residualVarianceFraction, cfg.parallel.useParallel, ...
        elapsedSeconds, all(checks.passed), ...
        mean(fitDiagnostics.convergenceCertified), ...
        'VariableNames', {'version', 'preset', 'nSites', 'nOuterFolds', ...
        'parameterSubsetSize', 'nSubsetCandidatesPerFold', ...
        'nuggetRule', ...
        'nuggetResidualVarianceFraction', 'parallelRequested', ...
        'elapsedSeconds', 'allImplementationChecksPassed', ...
        'optimizerConvergenceFraction'});

    writetable(foldManifest, fullfile(cfg.paths.subsetPluginDir, ...
        'outer_fold_manifest.csv'));
    writetable(foldSummary, fullfile(cfg.paths.subsetPluginDir, ...
        'outer_fold_summary.csv'));
    writetable(metrics, fullfile(cfg.paths.subsetPluginDir, ...
        'fold_metrics.csv'));
    writetable(methodSummary, fullfile(cfg.paths.subsetPluginDir, ...
        'method_summary.csv'));
    writetable(predictions, fullfile(cfg.paths.subsetPluginDir, ...
        'predictions.csv'));
    writetable(fitDiagnostics, fullfile(cfg.paths.subsetPluginDir, ...
        'fit_diagnostics.csv'));
    writetable(subsetSensitivity, fullfile(cfg.paths.subsetPluginDir, ...
        'subset_seed_sensitivity.csv'));
    writetable(parameterManifest, fullfile(cfg.paths.subsetPluginDir, ...
        'parameter_subset_manifest.csv'));
    writetable(checks, fullfile(cfg.paths.subsetPluginDir, ...
        'implementation_checks.csv'));
    writetable(timing, fullfile(cfg.paths.subsetPluginDir, ...
        'task_timing.csv'));
    writetable(runSummary, fullfile(cfg.paths.subsetPluginDir, ...
        'run_summary.csv'));
    save(fullfile(cfg.paths.subsetPluginDir, ...
        'subset_plugin_fits.mat'), 'tasks', 'cfg', 'folds', ...
        'nuggetRule', 'subset', 'rffStage0Inputs', '-v7.3');

    fprintf('\nSubset plug-in CV complete in %.1f minutes.\n', ...
        elapsedSeconds / 60);
    fprintf('  Nugget rule: %s.\n', char(nuggetRule.mode));
    fprintf('  Results: %s\n\n', cfg.paths.subsetPluginDir);
    result.cfg = cfg;
    result.metrics = metrics;
    result.methodSummary = methodSummary;
    result.fitDiagnostics = fitDiagnostics;
    result.subsetSensitivity = subsetSensitivity;
    result.checks = checks;
    result.runSummary = runSummary;
end

function inputs = make_rff_stage0_inputs(tasks, subset)
% Store only information that can be used without target responses.

    inputs = cell(numel(tasks), 1);
    for j = 1:numel(tasks)
        task = tasks{j};
        input.foldID = task.timing.foldID;
        input.Yobserved = subset.U(:, task.observedIndex);
        input.nTargetSites = numel(task.targetIndex);
        input.outerGeometry = task.outerGeometry;
        input.targetGeometry = task.targetGeometry;
        input.exactFit = task.fullFit;
        input.exactPrediction = task.fullPrediction;
        inputs{j} = input;
    end
end

function value = concatenate(tasks, name)
    pieces = cellfun(@(x) x.(name), tasks, 'UniformOutput', false);
    value = vertcat(pieces{:});
end

function summary = summarize_metrics(metrics)
    methods = unique(metrics.method, 'stable');
    rows = cell(numel(methods), 1);
    for j = 1:numel(methods)
        use = metrics.method == methods(j);
        rows{j} = table(methods(j), sum(use), ...
            finite_mean(metrics.rmse(use)), finite_se(metrics.rmse(use)), ...
            finite_mean(metrics.mae(use)), ...
            finite_mean(metrics.correlation(use)), ...
            finite_mean(metrics.nlpd(use)), ...
            finite_mean(metrics.crps(use)), ...
            finite_mean(metrics.coverage95(use)), ...
            finite_mean(metrics.meanIntervalWidth95(use)), ...
            'VariableNames', {'method', 'nTargetSites', 'meanRmse', ...
            'seRmseAcrossSites', 'meanMae', 'meanCorrelation', ...
            'meanNlpd', 'meanCrps', 'meanCoverage95', ...
            'meanIntervalWidth95'});
    end
    summary = vertcat(rows{:});
end

function value = finite_mean(x)
    x = x(isfinite(x));
    value = mean(x);
end

function value = finite_se(x)
    x = x(isfinite(x));
    if numel(x) < 2
        value = NaN;
    else
        value = std(x, 0) / sqrt(numel(x));
    end
end

function available = parallel_available()
    available = false;
    try
        available = license('test', 'Distrib_Computing_Toolbox');
    catch
        available = false;
    end
end
