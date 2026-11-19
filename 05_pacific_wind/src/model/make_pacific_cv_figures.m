function make_pacific_cv_figures(subset, metrics, blockMetrics, ...
        predictions, fitDiagnostics, varianceSummary, cfg)
%MAKE_PACIFIC_CV_FIGURES Saved primitive prediction/calibration figures.

    figureDir = cfg.paths.figureDir;
    visibility = cfg.figures.visible;
    methods = unique(metrics.method, 'stable');
    labels = replace(methods, '_', ' ');
    meanRmse = nan(numel(methods), 1);
    seRmse = nan(numel(methods), 1);
    for j = 1:numel(methods)
        values = metrics.rmse(metrics.method == methods(j));
        meanRmse(j) = mean(values);
        if numel(values) > 1
            seRmse(j) = std(values) / sqrt(numel(values));
        end
    end
    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 850, 540]);
    errorbar(1:numel(methods), meanRmse, seRmse, 'o', ...
        'LineWidth', 1.6, 'MarkerSize', 8);
    xticks(1:numel(methods)); xticklabels(labels); xtickangle(20);
    ylabel('Held-out RMSE (m/s)'); grid on;
    title('Pacific local-subset reconstruction');
    pacific_save_figure(h, figureDir, 'rmse_by_method');

    probabilisticMethods = ["fitted_eq9"; "separable_eq9"];
    coverage = nan(2, 1);
    for j = 1:2
        coverage(j) = mean(metrics.coverage95( ...
            metrics.method == probabilisticMethods(j)));
    end
    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 760, 520]);
    plot(1:2, coverage, 'o', 'LineWidth', 1.6, 'MarkerSize', 9);
    hold on; yline(0.95, '--k', 'Nominal 95%');
    xticks(1:2); xticklabels(replace(probabilisticMethods, '_', ' '));
    xlim([0.5, 2.5]); ylim([max(0, min(coverage) - 0.08), 1]);
    ylabel('Empirical 95% coverage'); grid on;
    title('Pooled predictive interval calibration');
    pacific_save_figure(h, figureDir, 'coverage_by_method');

    exampleStation = metrics.targetStation(1);
    useTruth = predictions.targetStation == exampleStation & ...
        predictions.method == "fitted_eq9";
    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 1100, 560]);
    plot(predictions.relativeHours(useTruth) / 24, ...
        predictions.truthUWindMps(useTruth), 'k-', 'LineWidth', 1.2);
    hold on;
    for j = 1:numel(methods)
        use = predictions.targetStation == exampleStation & ...
            predictions.method == methods(j);
        plot(predictions.relativeHours(use) / 24, ...
            predictions.predictionUWindMps(use), 'LineWidth', 1.0);
    end
    yline(0, ':k'); grid on;
    xlabel('Relative day'); ylabel('u-wind (m/s)');
    title('Example held-out reconstruction: ' + exampleStation);
    legend(["truth"; labels], 'Location', 'best');
    pacific_save_figure(h, figureDir, ...
        'example_station_reconstruction');

    use = predictions.method == "fitted_eq9" & ...
        isfinite(predictions.standardizedResidual);
    z = predictions.standardizedResidual(use);
    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 1080, 470]);
    layout = tiledlayout(1, 2, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    nexttile;
    histogram(z, 45, 'Normalization', 'pdf', 'EdgeColor', 'none');
    hold on;
    xx = linspace(max(-5, min(z)), min(5, max(z)), 400);
    plot(xx, exp(-0.5 * xx .^ 2) / sqrt(2 * pi), ...
        'k-', 'LineWidth', 1.5);
    xlabel('Standardized residual'); ylabel('Density'); grid on;
    title('Histogram and N(0,1)');
    nexttile;
    ordered = sort(z);
    n = numel(ordered);
    probability = ((1:n)' - 0.5) / n;
    theoretical = sqrt(2) * erfinv(2 * probability - 1);
    plot(theoretical, ordered, '.', 'MarkerSize', 4); hold on;
    limits = [min([theoretical; ordered]), max([theoretical; ordered])];
    plot(limits, limits, 'k--'); axis equal; xlim(limits); ylim(limits);
    xlabel('N(0,1) quantile'); ylabel('Empirical quantile'); grid on;
    title('Normal QQ diagnostic');
    title(layout, 'Fitted Equation-9 standardized residuals');
    pacific_save_figure(h, figureDir, ...
        'standardized_residual_distribution');

    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 900, 540]);
    hold on;
    for j = 1:2
        blockMean = nan(cfg.cv.timeBlocks, 1);
        for b = 1:cfg.cv.timeBlocks
            use = blockMetrics.method == probabilisticMethods(j) & ...
                blockMetrics.timeBlock == b;
            blockMean(b) = mean(blockMetrics.coverage95(use));
        end
        plot(1:cfg.cv.timeBlocks, blockMean, 'o-', ...
            'LineWidth', 1.5, 'MarkerSize', 7);
    end
    yline(0.95, '--k'); xlim([1, cfg.cv.timeBlocks]);
    xticks(1:cfg.cv.timeBlocks); xlabel('Consecutive 30-day block');
    ylabel('Empirical 95% coverage'); grid on;
    title('Calibration stability through the four-month record');
    legend(replace(probabilisticMethods, '_', ' '), ...
        'Location', 'best');
    pacific_save_figure(h, figureDir, 'coverage_by_time_block');

    fullFits = fitDiagnostics(fitDiagnostics.method == "fitted_eq9", :);
    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 780, 530]);
    scatter(fullFits.alphaPer100Km, fullFits.epsilon, 70, ...
        fullFits.timeScaleHours, 'filled');
    colorbar; grid on; yline(1, ':k');
    xlabel('\alpha per 100 km'); ylabel('\epsilon');
    title('Fitted parameter stability; color = time scale (hours)');
    pacific_save_figure(h, figureDir, ...
        'fitted_alpha_epsilon_diagnostic');

    fullMetrics = metrics(metrics.method == "fitted_eq9", :);
    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 850, 560]);
    scatter(subset.siteTable.xKm, subset.siteTable.yKm, 20, ...
        [0.7, 0.7, 0.7], 'filled');
    hold on;
    scatter(fullMetrics.xKm, fullMetrics.yKm, 110, ...
        fullMetrics.rmse, 'filled');
    axis equal; grid on; colorbar;
    xlabel('East-west km'); ylabel('North-south km');
    title('Held-out fitted Equation-9 RMSE');
    pacific_save_figure(h, figureDir, 'heldout_rmse_map');

    manifest = subset.siteTable;
    [~, foldManifest] = make_pacific_cv_folds(subset, ...
        cfg.cv.holdoutPattern);
    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 780, 590]);
    scatter(manifest.xKm, manifest.yKm, 95, foldManifest.foldID, ...
        'filled');
    axis equal; grid on; colorbar;
    xlabel('East-west km'); ylabel('North-south km');
    title('Balanced whole-location CV folds');
    pacific_save_figure(h, figureDir, 'cv_fold_assignment');

    probabilisticVariance = varianceSummary(ismember( ...
        varianceSummary.method, probabilisticMethods), :);
    components = [probabilisticVariance.meanConditionalProcessVariance, ...
        probabilisticVariance.meanFixedEffectVariance, ...
        probabilisticVariance.meanNuggetVariance];
    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 860, 560]);
    bar(1:height(probabilisticVariance), components, 'stacked');
    hold on;
    plot(1:height(probabilisticVariance), ...
        probabilisticVariance.empiricalMse, 'ko', ...
        'LineWidth', 1.8, 'MarkerSize', 8);
    xticks(1:height(probabilisticVariance));
    xticklabels(replace(probabilisticVariance.method, '_', ' '));
    xtickangle(15); ylabel('Variance (m^2/s^2)'); grid on;
    legend({'conditional process', 'fixed-effect estimation', ...
        'nugget', 'empirical squared error'}, 'Location', 'best');
    title('Predictive variance decomposition');
    pacific_save_figure(h, figureDir, ...
        'predictive_variance_decomposition');

    fittedPredictions = predictions(predictions.method == ...
        "fitted_eq9", :);
    targets = unique(fittedPredictions.targetStation, 'stable');
    meanVariance = nan(numel(targets), 1);
    empiricalMse = nan(numel(targets), 1);
    for j = 1:numel(targets)
        use = fittedPredictions.targetStation == targets(j);
        meanVariance(j) = mean( ...
            fittedPredictions.observationVariance(use));
        empiricalMse(j) = mean(fittedPredictions.residual(use) .^ 2);
    end
    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 760, 570]);
    scatter(meanVariance, empiricalMse, 60, 'filled'); hold on;
    limits = [0, 1.05 * max([meanVariance; empiricalMse])];
    if limits(2) <= 0 || ~isfinite(limits(2)), limits = [0, 1]; end
    plot(limits, limits, 'k--'); xlim(limits); ylim(limits);
    axis square; grid on;
    xlabel('Mean predicted observation variance');
    ylabel('Empirical mean squared error');
    title('Variance calibration by held-out site');
    pacific_save_figure(h, figureDir, ...
        'variance_calibration_by_site');
end
