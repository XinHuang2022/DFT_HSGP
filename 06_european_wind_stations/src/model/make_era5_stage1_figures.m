function make_era5_stage1_figures(summary,cfg)
%MAKE_ERA5_STAGE1_FIGURES Numerical convergence and same-theta prediction.
    f = figure('Visible',cfg.figures.visible,'Color','w','Position',[100 100 1100 430]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    nexttile; semilogy(summary.M,summary.p90LatentCovarianceError,'o-', ...
        summary.M,summary.maximumLatentCovarianceError,'s-','LineWidth',1.5);
    yline(cfg.hsgp.acceptance.maximumP90LatentCovarianceError,'--');
    xlabel('M=J^2'); ylabel('relative covariance error'); grid on;
    legend('90th percentile','maximum','90th-percentile tolerance','Location','best');
    title('Fixed-parameter HSGP covariance audit');
    nexttile; yyaxis left; plot(summary.M,summary.hsgpSameThetaRmse,'o-', ...
        'LineWidth',1.5); yline(summary.exactRmse(1),'--'); ylabel('RMSE (m/s)');
    yyaxis right; plot(summary.M,abs(summary.nllDifferencePerObservation),'s-', ...
        'LineWidth',1.5); ylabel('|NLL difference| per observation');
    xlabel('M=J^2'); grid on; title('Same-theta predictive and likelihood agreement');
    era5_save_figure(f,fullfile(cfg.paths.stage1Dir,'figures'),'stage1_J_audit');
end
