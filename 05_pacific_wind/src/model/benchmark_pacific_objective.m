function benchmark = benchmark_pacific_objective(preset, repetitions)
%BENCHMARK_PACIFIC_OBJECTIVE Time likelihood evaluations without fitting.
%
%   benchmark_pacific_objective('local7_cv5', 2)
%   benchmark_pacific_objective('full17_cv5', 1)
%
% This function deliberately bypasses the full-grid optimization gate.  It
% evaluates each profiled Whittle-REML objective at a reproducible interior
% parameter vector and writes objective_benchmark.csv.  Use its per-call
% timings together with the function-evaluation counts from a completed
% smaller run to estimate an HPC wall-time request.

    if nargin < 1 || isempty(preset)
        preset = 'local7_cv5';
    end
    if nargin < 2 || isempty(repetitions)
        repetitions = 1;
    end
    validateattributes(repetitions, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    cfg = pacific_eq9_config(preset);
    if ~cfg.runModel
        error('Choose a modeling preset for the objective benchmark.');
    end

    audit = run_pacific_audit(cfg.preset);
    if ~all(audit.dataChecks.passed)
        error('Pacific data-contract checks failed; inspect audit outputs.');
    end
    subset = make_pacific_subset(audit.contract, cfg.subset.name);
    [folds, ~, ~] = make_pacific_cv_folds( ...
        subset, cfg.cv.holdoutPattern);
    if isempty(folds)
        error('The selected preset has no CV folds.');
    end
    target = folds(1).targetIndex;
    observed = setdiff((1:height(subset.siteTable))', target, 'stable');
    Y = subset.U(:, observed);
    xyKm = subset.xyKm(observed, :);
    center = mean(xyKm, 1);
    scale = std(xyKm, 0, 1);
    scale(scale <= sqrt(eps)) = 1;
    design = make_pacific_mean_design(size(Y, 1), ...
        (xyKm - center) ./ scale);
    prepared = prepare_pacific_whittle_data(Y, ...
        (xyKm - center) ./ cfg.spatial.coordinateScaleKm, design);

    trainingVariance = var(Y, 0, 'all');
    theta = cfg.start;
    theta.signalVariance = interior_value_benchmark( ...
        0.8 * trainingVariance, cfg.bounds.signalVariance);
    theta.noiseVariance = interior_value_benchmark( ...
        0.2 * trainingVariance, cfg.bounds.noiseVariance);
    zFull = eq9_pack_parameters(theta, cfg);
    theta.epsilon = 1;
    zSeparable = pacific_pack_fixed_epsilon(theta, cfg);

    fullSeconds = nan(repetitions, 1);
    separableSeconds = nan(repetitions, 1);
    fullObjective = NaN;
    separableObjective = NaN;
    for j = 1:repetitions
        timerValue = tic;
        fullObjective = pacific_profiled_whittle_reml( ...
            zFull, prepared, cfg, []);
        fullSeconds(j) = toc(timerValue);
        timerValue = tic;
        separableObjective = pacific_profiled_whittle_reml( ...
            zSeparable, prepared, cfg, 1);
        separableSeconds(j) = toc(timerValue);
    end

    benchmark = table(string(cfg.version), string(cfg.preset), ...
        string(cfg.subset.name), height(subset.siteTable), ...
        numel(observed), numel(target), prepared.T, ...
        numel(prepared.usedIndices), repetitions, ...
        median(fullSeconds), median(separableSeconds), ...
        min(fullSeconds), min(separableSeconds), ...
        fullObjective, separableObjective, ...
        'VariableNames', {'version', 'preset', 'subsetName', ...
        'nSubsetSites', 'nTrainingSites', 'nTargetSites', ...
        'nTimePoints', 'nUsedFrequencies', 'repetitions', ...
        'medianFullSecondsPerEvaluation', ...
        'medianSeparableSecondsPerEvaluation', ...
        'minimumFullSecondsPerEvaluation', ...
        'minimumSeparableSecondsPerEvaluation', ...
        'fullObjective', 'separableObjective'});
    ensure_directory(cfg.paths.outputDir);
    writetable(benchmark, fullfile(cfg.paths.outputDir, ...
        'objective_benchmark.csv'));
    fprintf(['%s: %d training sites, %d frequencies; median full ', ...
        'objective %.3f s, separable %.3f s.\n'], cfg.preset, ...
        numel(observed), numel(prepared.usedIndices), ...
        median(fullSeconds), median(separableSeconds));
end

function value = interior_value_benchmark(value, bounds)
    value = min(bounds(2) / 1.01, max(bounds(1) * 1.01, value));
end
