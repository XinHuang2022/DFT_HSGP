function result = run_pacific_nugget_study(configOrPreset, useParallel)
%RUN_PACIFIC_NUGGET_STUDY Nested study of the Pacific nugget lower bound.

% The decision is based only on validation sites nested inside each outer
% training fold.  The outer CV targets provide an unbiased post-selection
% assessment and never choose sigma_n^2.

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
    ensure_directory(cfg.paths.nuggetDir);
    contract = load_pacific_wind_data(cfg);
    if any(~isfinite(contract.U), 'all')
        error('The nugget study requires a complete Pacific panel.');
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
            tasks{f} = run_pacific_nugget_fold(cfg, subset, ...
                folds(f).targetIndex, folds(f).foldID);
        end
    else
        for f = 1:nFold
            fprintf('Nugget study fold %d/%d.\n', f, nFold);
            tasks{f} = run_pacific_nugget_fold(cfg, subset, ...
                folds(f).targetIndex, folds(f).foldID);
        end
    end
    elapsedSeconds = toc(timerValue);

    profile = concatenate(tasks, 'profile');
    selections = concatenate(tasks, 'selection');
    outerMetrics = concatenate(tasks, 'outerMetrics');
    parameterManifest = concatenate(tasks, 'parameterManifest');
    innerManifest = concatenate(tasks, 'innerManifest');
    timing = concatenate(tasks, 'timing');
    [candidateSummary, decision] = summarize_decision(profile, cfg);
    runSummary = table(string(cfg.version), string(cfg.preset), ...
        height(subset.siteTable), nFold, ...
        cfg.scalable.parameterSubsetSize, ...
        numel(cfg.nugget.residualVarianceFractions) + 1, ...
        cfg.parallel.useParallel, elapsedSeconds, ...
        string(decision.recommendedMode), ...
        decision.recommendedResidualVarianceFraction, ...
        'VariableNames', {'version', 'preset', 'nSites', 'nOuterFolds', ...
        'parameterSubsetSize', 'nNuggetCandidates', 'parallelRequested', ...
        'elapsedSeconds', 'recommendedMode', ...
        'recommendedResidualVarianceFraction'});

    writetable(foldManifest, fullfile(cfg.paths.nuggetDir, ...
        'outer_fold_manifest.csv'));
    writetable(foldSummary, fullfile(cfg.paths.nuggetDir, ...
        'outer_fold_summary.csv'));
    writetable(profile, fullfile(cfg.paths.nuggetDir, ...
        'nugget_profile.csv'));
    writetable(candidateSummary, fullfile(cfg.paths.nuggetDir, ...
        'nugget_candidate_summary.csv'));
    writetable(decision, fullfile(cfg.paths.nuggetDir, ...
        'nugget_decision.csv'));
    writetable(selections, fullfile(cfg.paths.nuggetDir, ...
        'nugget_fold_selection.csv'));
    writetable(outerMetrics, fullfile(cfg.paths.nuggetDir, ...
        'outer_postselection_metrics.csv'));
    writetable(parameterManifest, fullfile(cfg.paths.nuggetDir, ...
        'parameter_subset_manifest.csv'));
    writetable(innerManifest, fullfile(cfg.paths.nuggetDir, ...
        'inner_split_manifest.csv'));
    writetable(timing, fullfile(cfg.paths.nuggetDir, 'task_timing.csv'));
    writetable(runSummary, fullfile(cfg.paths.nuggetDir, ...
        'run_summary.csv'));
    save(fullfile(cfg.paths.nuggetDir, 'nugget_stage_fits.mat'), ...
        'tasks', 'cfg', 'folds', '-v7.3');
    make_nugget_figure(candidateSummary, cfg);

    fprintf('\nPacific nugget study complete in %.1f minutes.\n', ...
        elapsedSeconds / 60);
    fprintf('  Recommendation: %s', char(decision.recommendedMode));
    if isfinite(decision.recommendedResidualVarianceFraction)
        fprintf(' at residual-variance fraction %.4g', ...
            decision.recommendedResidualVarianceFraction);
    end
    fprintf('.\n  Results: %s\n\n', cfg.paths.nuggetDir);

    result.cfg = cfg;
    result.profile = profile;
    result.candidateSummary = candidateSummary;
    result.decision = decision;
    result.selections = selections;
    result.outerMetrics = outerMetrics;
    result.runSummary = runSummary;
end

function value = concatenate(tasks, name)
    pieces = cellfun(@(x) x.(name), tasks, 'UniformOutput', false);
    value = vertcat(pieces{:});
end

function [summary, decision] = summarize_decision(profile, cfg)
    labels = unique(profile.candidate, 'stable');
    free = profile(profile.candidate == "free_current", :);
    rows = cell(numel(labels), 1);
    for j = 1:numel(labels)
        use = profile.candidate == labels(j);
        part = profile(use, :);
        [common, ia, ib] = intersect(part.foldID, free.foldID, 'stable');
        delta = part.validationNlpd(ia) - free.validationNlpd(ib);
        rows{j} = table(labels(j), part.residualVarianceFraction(1), ...
            numel(unique(part.foldID)), mean(part.validationNlpd), ...
            standard_error(part.validationNlpd), ...
            mean(part.validationRmse), mean(part.validationCoverage95), ...
            mean(delta), standard_error(delta), numel(common), ...
            'VariableNames', {'candidate', 'residualVarianceFraction', ...
            'nFolds', 'meanValidationNlpd', 'seValidationNlpd', ...
            'meanValidationRmse', 'meanValidationCoverage95', ...
            'meanDeltaNlpdVersusFree', 'seDeltaNlpdVersusFree', ...
            'nPairedFolds'});
    end
    summary = vertcat(rows{:});
    floors = summary(summary.candidate ~= "free_current", :);
    [~, best] = min(floors.meanValidationNlpd);
    bestFloor = floors(best, :);
    upper = bestFloor.meanDeltaNlpdVersusFree + ...
        cfg.nugget.approximateConfidenceMultiplier * ...
        bestFloor.seDeltaNlpdVersusFree;
    if isfinite(upper) && upper < 0
        mode = "larger_scale_adaptive_floor_supported";
        recommended = bestFloor.residualVarianceFraction;
    else
        mode = "retain_current_numerical_floor";
        recommended = NaN;
    end
    decision = table(mode, recommended, bestFloor.candidate, ...
        bestFloor.residualVarianceFraction, ...
        bestFloor.meanDeltaNlpdVersusFree, ...
        bestFloor.seDeltaNlpdVersusFree, upper, ...
        cfg.nugget.approximateConfidenceMultiplier, ...
        'VariableNames', {'recommendedMode', ...
        'recommendedResidualVarianceFraction', 'bestFloorCandidate', ...
        'bestFloorResidualVarianceFraction', ...
        'bestFloorMeanDeltaNlpdVersusFree', ...
        'bestFloorSeDeltaNlpdVersusFree', ...
        'bestFloorUpperUncertaintyBound', ...
        'uncertaintyMultiplier'});
end

function value = standard_error(x)
    x = x(isfinite(x));
    if numel(x) < 2
        value = NaN;
    else
        value = std(x, 0) / sqrt(numel(x));
    end
end

function make_nugget_figure(summary, cfg)
    floors = summary(summary.candidate ~= "free_current", :);
    figureHandle = figure('Visible', cfg.figures.visible, ...
        'Color', 'white', 'Position', [100, 100, 900, 650]);
    errorbar(floors.residualVarianceFraction, ...
        floors.meanValidationNlpd, 2 * floors.seValidationNlpd, ...
        '-o', 'LineWidth', 1.5, 'MarkerSize', 7);
    hold on;
    free = summary(summary.candidate == "free_current", :);
    yline(free.meanValidationNlpd, '--k', 'free nugget fit');
    set(gca, 'XScale', 'log');
    xlabel('Nugget lower bound / empirical residual variance');
    ylabel('Inner-validation NLPD (lower is better)');
    title('Training-only nugget sensitivity');
    grid on;
    figureDir = fullfile(cfg.paths.nuggetDir, 'figures');
    ensure_directory(figureDir);
    pacific_save_figure(figureHandle, figureDir, ...
        'nugget_validation_profile');
end

function available = parallel_available()
    available = false;
    try
        available = license('test', 'Distrib_Computing_Toolbox');
    catch
        available = false;
    end
end
