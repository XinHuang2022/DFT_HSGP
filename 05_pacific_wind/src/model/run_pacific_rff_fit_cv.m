function result = run_pacific_rff_fit_cv(configOrPreset, useParallel)
%RUN_PACIFIC_RFF_FIT_CV Gated all-training-site RFF likelihood fit.

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
    decisionFile = fullfile(cfg.paths.rffStage0Dir, ...
        'rff_stage0_decision.csv');
    if ~isfile(decisionFile)
        error('Run run_pacific_rff_stage0 before RFF optimization.');
    end
    stage0 = readtable(decisionFile, 'TextType', 'string');
    if height(stage0) ~= 1
        error('RFF stage-0 decision file must contain exactly one row.');
    end
    if ~ismember('targetResponseFree', stage0.Properties.VariableNames) || ...
            ~stage0.targetResponseFree
        error(['RFF stage 0 must certify that outer-target responses ', ...
            'were not used to select the feature count.']);
    end
    if cfg.rff.fitEnabledOnlyAfterStage0Pass && ~stage0.stage0Passed
        error(['RFF stage 0 did not pass. Increase/repair the feature ', ...
            'bank before optimizing an approximate likelihood.']);
    end
    if cfg.rff.requireAccelerationBeforeFit && ...
            ~stage0.accelerationSupportedAtRecommendedP
        error(['The fixed-theta RFF approximation passed its accuracy ', ...
            'gate but did not demonstrate acceleration at the accepted ', ...
            'feature count. Do not optimize the RFF likelihood unless ', ...
            'you deliberately set cfg.rff.requireAccelerationBeforeFit ', ...
            'to false after reviewing rff_stage0_decision.csv.']);
    end
    nFeatures = stage0.recommendedFeatureCount;
    if ~(isfinite(nFeatures) && nFeatures >= 4)
        error('RFF stage 0 did not provide a valid feature count.');
    end
    fitFile = fullfile(cfg.paths.subsetPluginDir, ...
        'subset_plugin_fits.mat');
    if ~isfile(fitFile)
        error('Subset plug-in fits are unavailable: %s', fitFile);
    end
    loaded = load(fitFile, 'tasks', 'subset');
    pluginTasks = loaded.tasks;
    subset = loaded.subset;
    nFold = numel(pluginTasks);
    ensure_directory(cfg.paths.rffFitDir);
    outputs = cell(nFold, 1);
    timerValue = tic;
    if cfg.parallel.useParallel && parallel_available()
        pool = gcp('nocreate');
        if isempty(pool) && cfg.parallel.autoStartPool
            pool = parpool(min(cfg.parallel.numWorkers, nFold));
        end
        parfor f = 1:nFold
            outputs{f} = run_pacific_rff_fit_fold( ...
                cfg, subset, pluginTasks{f}, nFeatures);
        end
    else
        for f = 1:nFold
            fprintf('RFF fit fold %d/%d.\n', f, nFold);
            outputs{f} = run_pacific_rff_fit_fold( ...
                cfg, subset, pluginTasks{f}, nFeatures);
        end
    end
    elapsedSeconds = toc(timerValue);
    diagnostics = concatenate(outputs, 'diagnostics');
    metrics = concatenate(outputs, 'metrics');
    targetMetrics = concatenate(outputs, 'targetMetrics');
    summary = summarize_metrics(metrics);
    targetSummary = summarize_metrics(targetMetrics);
    pairedMetrics = make_paired_differences(metrics, false);
    pairedTargetMetrics = make_paired_differences(targetMetrics, true);
    pairedSummary = summarize_paired_differences(pairedMetrics);
    pairedTargetSummary = ...
        summarize_paired_differences(pairedTargetMetrics);
    runSummary = table(string(cfg.version), string(cfg.preset), nFold, ...
        nFeatures, numel(cfg.rff.fitFeatureSeeds), ...
        cfg.rff.fixNuggetAtExactSubsetEstimateDuringFirstFit, ...
        cfg.parallel.useParallel, elapsedSeconds, ...
        mean(diagnostics.convergenceCertified), ...
        mean(diagnostics.fitPassedNumericalGates), ...
        'VariableNames', {'version', 'preset', 'nOuterFolds', ...
        'nFeatures', 'nFeatureSeeds', 'nuggetFixedDuringFirstRffFit', ...
        'parallelRequested', 'elapsedSeconds', ...
        'convergenceFraction', 'numericalGatePassFraction'});
    writetable(diagnostics, fullfile(cfg.paths.rffFitDir, ...
        'rff_fit_diagnostics.csv'));
    writetable(metrics, fullfile(cfg.paths.rffFitDir, ...
        'rff_fit_metrics.csv'));
    writetable(targetMetrics, fullfile(cfg.paths.rffFitDir, ...
        'rff_fit_target_metrics.csv'));
    writetable(summary, fullfile(cfg.paths.rffFitDir, ...
        'rff_fit_summary.csv'));
    writetable(targetSummary, fullfile(cfg.paths.rffFitDir, ...
        'rff_fit_target_summary.csv'));
    writetable(pairedMetrics, fullfile(cfg.paths.rffFitDir, ...
        'rff_fit_paired_metrics.csv'));
    writetable(pairedTargetMetrics, fullfile(cfg.paths.rffFitDir, ...
        'rff_fit_paired_target_metrics.csv'));
    writetable(pairedSummary, fullfile(cfg.paths.rffFitDir, ...
        'rff_fit_paired_summary.csv'));
    writetable(pairedTargetSummary, fullfile(cfg.paths.rffFitDir, ...
        'rff_fit_paired_target_summary.csv'));
    writetable(runSummary, fullfile(cfg.paths.rffFitDir, ...
        'run_summary.csv'));
    save(fullfile(cfg.paths.rffFitDir, 'rff_fits.mat'), ...
        'outputs', 'cfg', 'stage0', '-v7.3');
    make_fit_figure(targetMetrics, cfg);
    fprintf('\nRFF likelihood CV complete in %.1f minutes.\n', ...
        elapsedSeconds / 60);
    fprintf('  Feature count: %d; results: %s\n\n', ...
        nFeatures, cfg.paths.rffFitDir);
    result.cfg = cfg;
    result.diagnostics = diagnostics;
    result.metrics = metrics;
    result.targetMetrics = targetMetrics;
    result.summary = summary;
    result.targetSummary = targetSummary;
    result.pairedMetrics = pairedMetrics;
    result.pairedTargetMetrics = pairedTargetMetrics;
    result.pairedSummary = pairedSummary;
    result.pairedTargetSummary = pairedTargetSummary;
    result.runSummary = runSummary;
end

function make_fit_figure(metrics, cfg)
    baseline = metrics(metrics.method == "exact_subset_baseline", :);
    rff = metrics(metrics.method == "rff_theta_rff_prediction", :);
    keys = {'foldID', 'featureSeed', 'targetStation'};
    [common, ia, ib] = intersect(baseline(:, keys), rff(:, keys), ...
        'rows', 'stable'); %#ok<ASGLU>
    difference = rff.rmse(ib) - baseline.rmse(ia);
    h = figure('Visible', cfg.figures.visible, 'Color', 'white', ...
        'Position', [100, 100, 950, 650]);
    scatter(baseline.rmse(ia), difference, 35, ...
        baseline.xKm(ia), 'filled');
    yline(0, '--k');
    xlabel('Exact subset plug-in RMSE (m/s)');
    ylabel('RFF fit/prediction minus exact RMSE (m/s)');
    title('RFF approximation cost by held-out site');
    colorbar; grid on;
    figureDir = fullfile(cfg.paths.rffFitDir, 'figures');
    ensure_directory(figureDir);
    pacific_save_figure(h, figureDir, 'rff_fit_rmse_difference');
end

function value = concatenate(outputs, name)
    parts = cellfun(@(x) x.(name), outputs, 'UniformOutput', false);
    value = vertcat(parts{:});
end

function summary = summarize_metrics(metrics)
    methods = unique(metrics.method, 'stable');
    rows = cell(numel(methods), 1);
    for j = 1:numel(methods)
        use = metrics.method == methods(j);
        rows{j} = table(methods(j), sum(use), ...
            mean(metrics.rmse(use)), mean(metrics.nlpd(use)), ...
            mean(metrics.crps(use)), mean(metrics.coverage95(use)), ...
            mean(metrics.meanIntervalWidth95(use)), ...
            'VariableNames', {'method', 'nResults', ...
            'meanRmse', 'meanNlpd', 'meanCrps', 'meanCoverage95', ...
            'meanIntervalWidth95'});
    end
    summary = vertcat(rows{:});
end

function differences = make_paired_differences(metrics, includeTarget)
% Three contrasts identify parameter, feature, and total RFF costs.

    if includeTarget
        keys = {'foldID', 'featureSeed', 'targetStation'};
    else
        keys = {'foldID', 'featureSeed'};
    end
    baseline = metrics(metrics.method == "exact_subset_baseline", :);
    exactTheta = metrics( ...
        metrics.method == "rff_theta_exact_prediction", :);
    rff = metrics(metrics.method == "rff_theta_rff_prediction", :);
    comparisons = [ ...
        "parameter_displacement", ...
        "feature_approximation", ...
        "total_rff_cost"];
    methodA = [ ...
        "rff_theta_exact_prediction", ...
        "rff_theta_rff_prediction", ...
        "rff_theta_rff_prediction"];
    methodB = [ ...
        "exact_subset_baseline", ...
        "rff_theta_exact_prediction", ...
        "exact_subset_baseline"];
    sourcesA = {exactTheta, rff, rff};
    sourcesB = {baseline, exactTheta, baseline};
    rows = cell(numel(comparisons), 1);
    for j = 1:numel(comparisons)
        left = sourcesA{j};
        right = sourcesB{j};
        [common, ia, ib] = intersect(left(:, keys), right(:, keys), ...
            'rows', 'stable');
        nPair = height(common);
        comparison = repmat(comparisons(j), nPair, 1);
        a = repmat(methodA(j), nPair, 1);
        b = repmat(methodB(j), nPair, 1);
        part = [common, table(comparison, a, b, ...
            left.rmse(ia) - right.rmse(ib), ...
            left.nlpd(ia) - right.nlpd(ib), ...
            left.crps(ia) - right.crps(ib), ...
            left.coverage95(ia) - right.coverage95(ib), ...
            left.meanIntervalWidth95(ia) - ...
                right.meanIntervalWidth95(ib), ...
            'VariableNames', {'comparison', 'methodA', 'methodB', ...
            'rmseDifference', 'nlpdDifference', 'crpsDifference', ...
            'coverageDifference', 'intervalWidthDifference'})];
        rows{j} = part;
    end
    differences = vertcat(rows{:});
end

function summary = summarize_paired_differences(differences)
    comparisons = unique(differences.comparison, 'stable');
    rows = cell(numel(comparisons), 1);
    for j = 1:numel(comparisons)
        use = differences.comparison == comparisons(j);
        n = sum(use);
        rmse = differences.rmseDifference(use);
        nlpd = differences.nlpdDifference(use);
        crps = differences.crpsDifference(use);
        coverage = differences.coverageDifference(use);
        width = differences.intervalWidthDifference(use);
        rows{j} = table(comparisons(j), n, ...
            mean(rmse), standard_error(rmse), ...
            mean(nlpd), standard_error(nlpd), ...
            mean(crps), standard_error(crps), ...
            mean(coverage), standard_error(coverage), ...
            mean(width), standard_error(width), ...
            'VariableNames', {'comparison', 'nPairs', ...
            'meanRmseDifference', 'seRmseDifference', ...
            'meanNlpdDifference', 'seNlpdDifference', ...
            'meanCrpsDifference', 'seCrpsDifference', ...
            'meanCoverageDifference', 'seCoverageDifference', ...
            'meanIntervalWidthDifference', ...
            'seIntervalWidthDifference'});
    end
    summary = vertcat(rows{:});
end

function value = standard_error(x)
    if numel(x) > 1
        value = std(x, 0) / sqrt(numel(x));
    else
        value = NaN;
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
