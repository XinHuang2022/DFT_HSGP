function make_era5_stage2_figures(summary,cfg)
%MAKE_ERA5_STAGE2_FIGURES Parameter and prediction stability after refitting.
    f=figure('Visible',cfg.figures.visible,'Color','w','Position',[100 100 1150 700]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    nexttile;plot(summary.M,summary.hsgpEngineRmse,'o-',summary.M, ...
        summary.hsgpThetaExactEngineRmse,'s-',summary.M, ...
        summary.exactFitExactEngineRmse,'--','LineWidth',1.4);grid on;
    xlabel('M=J^2');ylabel('RMSE (m/s)');legend('HSGP fit/HSGP engine', ...
        'HSGP theta/exact engine','exact fit/exact engine','Location','best');
    title('Prediction decomposition');
    nexttile;semilogy(summary.M,abs(summary.hsgpMinusExactNllPerObservation), ...
        'o-','LineWidth',1.4);grid on;xlabel('M=J^2');ylabel('|NLL difference|/observation');
    title('Likelihood engine displacement');
    nexttile;plot(summary.M,summary.alphaPer100Km,'o-',summary.M, ...
        summary.betaPerSixHours,'s-','LineWidth',1.4);grid on;xlabel('M=J^2');
    legend('\alpha per 100 km','\beta per 6 h','Location','best');title('Scale estimates');
    nexttile;plot(summary.M,summary.epsilon,'o-',summary.M, ...
        summary.noiseVariance,'s-','LineWidth',1.4);grid on;xlabel('M=J^2');
    legend('\epsilon','\sigma_n^2','Location','best');title('Interaction and nugget');
    era5_save_figure(f,fullfile(cfg.paths.stage2Dir,'figures'),'stage2_hsgp_refit_screen');
end
