function make_cerra_dense_stage3_figures(summary,matched,cfg)
%MAKE_CERRA_DENSE_STAGE3_FIGURES Accuracy-cost frontier and paired differences.

    isHsgp=startsWith(string(summary.method),"hsgp");
    figure('Visible',cfg.figures.visible,'Color','w','Position',[100 100 1450 560]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile;scatter(summary.totalStartFitSeconds/3600,summary.pooledRmse,90, ...
        double(isHsgp),'filled');grid on;xlabel('Total start CPU time (h)');
    ylabel('Held-site pooled RMSE (m/s)');title('Accuracy-cost frontier');
    for j=1:height(summary)
        text(summary.totalStartFitSeconds(j)/3600,summary.pooledRmse(j), ...
            "  "+summary.method(j),'Interpreter','none');
    end
    nexttile;
    if isempty(matched),axis off;text(.1,.5,'No matched dimensions completed.');
    else
        errorbar(matched.matchedDimension, ...
            matched.meanSiteRmseDifferenceHsgpMinusExact, ...
            matched.meanSiteRmseDifferenceHsgpMinusExact-matched.lower95, ...
            matched.upper95-matched.meanSiteRmseDifferenceHsgpMinusExact, ...
            'o-','LineWidth',1.6);yline(0,'--');grid on;xlabel('Matched M=n_c');
        ylabel('Mean site RMSE: HSGP - exact (m/s)');title('Spatial-block bootstrap');
    end
    cerra_save_figure(gcf,fullfile(cfg.paths.stage3Dir,'figures', ...
        'accuracy_cost_and_matched_comparison.png'));close(gcf);
end
