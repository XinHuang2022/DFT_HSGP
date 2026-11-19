function make_pacific_hsgp_figures(metrics, predictions, ...
        fitDiagnostics, methodSummary, cfg)
%MAKE_PACIFIC_HSGP_FIGURES Core real-data HSGP diagnostics.

    figureDir = fullfile(cfg.paths.hsgpCvDir, 'figures');
    ensure_directory(figureDir);
    visibility = cfg.figures.visible;

    methods = methodSummary.method;
    labels = replace(methods, '_', ' ');
    meanRmse = methodSummary.meanRmse;
    [ciLow, ciHigh] = bootstrap_mean_intervals(metrics, methods, cfg);
    h = figure('Visible', visibility, 'Color', 'white', ...
        'Position', [100, 100, 980, 560]);
    bars = bar(meanRmse, 'FaceColor', 'flat');
    bars.CData = lines(numel(methods)); hold on;
    errorbar(1:numel(methods), meanRmse, meanRmse - ciLow, ...
        ciHigh - meanRmse, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
    xticks(1:numel(methods)); xticklabels(labels); xtickangle(18);
    ylabel('Mean held-out-site RMSE (m/s)'); grid on;
    title('Pacific HSGP cross-validation');
    ylim([0, 1.18 * max(ciHigh)]);
    pacific_save_figure(h, figureDir, 'hsgp_cv_rmse');

    full = fitDiagnostics(fitDiagnostics.method == "hsgp_eq9", :);
    h = figure('Visible', visibility, 'Color', 'white', ...
        'Position', [100, 100, 1020, 520]);
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    plot(full.foldID, full.hsgpObjectivePerObservation, ...
        'o-', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
    if any(isfinite(full.exactObjectiveAtHsgpTheta))
        exactPerObservation = full.hsgpObjectivePerObservation - ...
            full.hsgpMinusExactObjectivePerObservation;
        plot(full.foldID, exactPerObservation, 's--', ...
            'LineWidth', 1.5, 'MarkerSize', 7);
        legend({'HSGP objective', 'exact objective at HSGP theta'}, ...
            'Location', 'best');
    end
    xlabel('CV fold'); ylabel('NLL per training observation'); grid on;
    title('Likelihood audit');
    nexttile;
    semilogy(full.foldID, full.noiseVariance, 'o-', ...
        'LineWidth', 1.5, 'MarkerSize', 7); hold on;
    semilogy(full.foldID, full.trainingOnlyNuggetFloor, 'k--', ...
        'LineWidth', 1.3);
    xlabel('CV fold'); ylabel('Noise variance'); grid on;
    legend({'estimate', 'training-only floor'}, 'Location', 'best');
    title('Nugget estimate and lower bound');
    pacific_save_figure(h, figureDir, 'hsgp_likelihood_and_nugget');

    station = choose_example_station(metrics, cfg);
    primary = predictions(predictions.targetStation == station & ...
        predictions.method == "hsgp_eq9", :);
    truth = primary.truthUWindMps;
    day = primary.relativeHours / 24;
    sd = sqrt(primary.observationVariance);
    lower = primary.predictionUWindMps - 1.96 * sd;
    upper = primary.predictionUWindMps + 1.96 * sd;
    h = figure('Visible', visibility, 'Color', 'white', ...
        'Position', [100, 100, 1250, 520]);
    fill([day; flipud(day)], [lower; flipud(upper)], ...
        [0.85, 0.92, 0.98], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.65); hold on;
    plot(day, truth, 'k-', 'LineWidth', 1.1);
    plot(day, primary.predictionUWindMps, '-', ...
        'Color', [0, 0.35, 0.75], 'LineWidth', 1.3);
    exact = predictions(predictions.targetStation == station & ...
        predictions.method == "hsgp_theta_exact_engine", :);
    legendText = ["HSGP 95% observation band"; "observed wind"; ...
        "HSGP reconstruction"];
    if ~isempty(exact)
        plot(exact.relativeHours / 24, exact.predictionUWindMps, ...
            '--', 'Color', [0.85, 0.3, 0.05], 'LineWidth', 1.2);
        legendText(end + 1) = "exact engine at HSGP parameters";
    end
    yline(0, ':k'); grid on;
    xlabel('Relative day'); ylabel('u-wind (m/s)');
    title('Held-out full-series reconstruction: ' + station);
    legend(legendText, 'Location', 'best');
    pacific_save_figure(h, figureDir, ...
        'hsgp_example_station_reconstruction');
end

function [low, high] = bootstrap_mean_intervals(metrics, methods, cfg)
    rng(cfg.hsgp.bootstrapSeed, 'twister');
    B = cfg.hsgp.bootstrapReplicates;
    low = nan(numel(methods), 1);
    high = nan(numel(methods), 1);
    for j = 1:numel(methods)
        values = metrics.rmse(metrics.method == methods(j));
        values = values(isfinite(values));
        n = numel(values);
        draws = mean(values(randi(n, n, B)), 1);
        interval = prctile(draws, [2.5, 97.5]);
        low(j) = interval(1);
        high(j) = interval(2);
    end
end

function station = choose_example_station(metrics, cfg)
    if strlength(string(cfg.hsgp.exampleStation)) > 0
        station = string(cfg.hsgp.exampleStation);
        if ~any(metrics.targetStation == station)
            error('Configured example station is not a CV target.');
        end
        return;
    end
    primary = metrics(metrics.method == "hsgp_eq9", :);
    [~, index] = min(primary.xKm .^ 2 + primary.yKm .^ 2);
    station = primary.targetStation(index);
end
