function make_cerra_dense_stage1_figure(summary,exactRmse,cfg)
%MAKE_CERRA_DENSE_STAGE1_FIGURE Fixed-theta screening diagnostics.

    figure('Visible',cfg.figures.visible,'Color','w','Position',[100 100 1500 560]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile;plot(summary.M,summary.pooledRmse,'-o','LineWidth',1.8);hold on;
    yline(exactRmse,'--','exact n=250');grid on;xlabel('M');ylabel('RMSE (m/s)');
    title(sprintf('Fixed-theta prediction, T=%d',cfg.data.expectedTimePoints));
    nexttile;semilogy(summary.M,summary.p90LatentCovarianceError,'-o','LineWidth',1.8);
    grid on;xlabel('M');ylabel('90th-percentile relative covariance error');
    title('Spatial covariance audit on 250 sites');
    cerra_save_figure(gcf,fullfile(cfg.paths.stage1Dir,'figures', ...
        'stage1_fixed_theta_dense_screen.png'));
    close(gcf);
end
