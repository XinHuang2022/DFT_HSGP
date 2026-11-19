function make_pacific_audit_figures(contract, siteSummary, timeSummary, ...
        blockSiteSummary, acfTable, cfg)
%MAKE_PACIFIC_AUDIT_FIGURES Primitive missingness/stationarity figures.

    figureDir = cfg.paths.figureDir;
    visibility = cfg.figures.visible;
    x = contract.siteTable.xKm;
    y = contract.siteTable.yKm;

    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 1100, 470]);
    layout = tiledlayout(1, 2, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    nexttile;
    scatter(x, y, 55, siteSummary.siteMean, 'filled');
    axis equal; grid on; colorbar;
    xlabel('Local east-west distance (km)');
    ylabel('Local north-south distance (km)');
    title('Long-run site mean (m/s)');
    nexttile;
    scatter(x, y, 55, siteSummary.siteSd, 'filled');
    axis equal; grid on; colorbar;
    xlabel('Local east-west distance (km)');
    ylabel('Local north-south distance (km)');
    title('Long-run site SD (m/s)');
    title(layout, 'Pacific u-wind spatial heterogeneity');
    pacific_save_figure(h, figureDir, 'site_mean_and_sd_maps');

    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 1100, 620]);
    layout = tiledlayout(2, 1, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    nexttile;
    plot(timeSummary.relativeHours / 24, ...
        timeSummary.crossSectionMean, 'LineWidth', 1.1);
    yline(0, '--k'); grid on;
    ylabel('Cross-sectional mean (m/s)');
    title('Spatially averaged u-wind');
    nexttile;
    plot(timeSummary.relativeHours / 24, ...
        timeSummary.crossSectionSd, 'LineWidth', 1.1);
    grid on; xlabel('Relative day');
    ylabel('Cross-sectional SD (m/s)');
    title('Spatial amplitude through time');
    title(layout, 'Primitive temporal nonstationarity diagnostic');
    pacific_save_figure(h, figureDir, 'temporal_mean_and_spatial_sd');

    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 920, 520]);
    plot(acfTable.lagHours, acfTable.averageSiteAcf, ...
        'LineWidth', 1.6);
    yline(0, '--k'); grid on;
    xlabel('Temporal lag (hours)');
    ylabel('Average site ACF');
    title('Average temporal autocorrelation before detrending');
    pacific_save_figure(h, figureDir, 'average_site_acf');

    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 1150, 820]);
    layout = tiledlayout(2, 2, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    colorLimit = [min(blockSiteSummary.blockSiteSd), ...
        max(blockSiteSummary.blockSiteSd)];
    for b = 1:4
        nexttile;
        use = blockSiteSummary.block == b;
        scatter(blockSiteSummary.xKm(use), blockSiteSummary.yKm(use), ...
            45, blockSiteSummary.blockSiteSd(use), 'filled');
        axis equal; grid on; colorbar; caxis(colorLimit);
        title(sprintf('Time block %d: site SD', b));
        xlabel('East-west km'); ylabel('North-south km');
    end
    title(layout, 'Spatial variance pattern across four time blocks');
    pacific_save_figure(h, figureDir, 'time_block_site_sd_maps');

    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 920, 520]);
    missingBySite = siteSummary.missingCount;
    scatter(x, y, 55, missingBySite, 'filled');
    axis equal; grid on; colorbar;
    xlabel('Local east-west distance (km)');
    ylabel('Local north-south distance (km)');
    title('Missing observations by site');
    pacific_save_figure(h, figureDir, 'missingness_map');

    exampleGrid = [9, 9; 8, 9; 10, 9; 9, 8; 9, 10];
    h = figure('Visible', visibility, 'Color', 'w', ...
        'Position', [100, 100, 1100, 540]);
    hold on;
    legendText = strings(size(exampleGrid, 1), 1);
    for j = 1:size(exampleGrid, 1)
        use = contract.siteTable.gridRow == exampleGrid(j, 1) & ...
            contract.siteTable.gridColumn == exampleGrid(j, 2);
        originalColumn = contract.siteTable.originalColumn(use);
        plot(contract.relativeHours / 24, ...
            contract.U(:, originalColumn), 'LineWidth', 0.9);
        legendText(j) = contract.siteTable.stationID(use);
    end
    yline(0, ':k'); grid on;
    xlabel('Relative day'); ylabel('u-wind (m/s)');
    title('Central grid-point time series');
    legend(legendText, 'Location', 'best');
    pacific_save_figure(h, figureDir, 'central_grid_time_series');
end
