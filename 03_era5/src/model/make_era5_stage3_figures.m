function make_era5_stage3_figures(metrics,predictions,fitDiagnostics,summary,selectedJ,cfg)
%MAKE_ERA5_STAGE3_FIGURES Publication-oriented CV, timing, and trace plots.
    figureDir=fullfile(cfg.paths.stage3Dir,'figures');ensure_directory(figureDir);
    f=figure('Visible',cfg.figures.visible,'Color','w','Position',[100 100 1050 520]);
    x=1:height(summary);bar(x,summary.meanRmse,'FaceColor',[0.2 0.55 0.75]);hold on;
    errorbar(x,summary.meanRmse,summary.meanRmse-summary.bootstrapLower95, ...
        summary.bootstrapUpper95-summary.meanRmse,'k.','LineWidth',1.5);
    set(gca,'XTick',x,'XTickLabel',pretty(summary.method),'XTickLabelRotation',20);
    ylabel('Mean held-site RMSE (m/s)');title('North Sea ERA5 five-fold reconstruction');grid on;
    era5_save_figure(f,figureDir,'cv_rmse_with_bootstrap');

    f=figure('Visible',cfg.figures.visible,'Color','w','Position',[100 100 1150 470]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    nexttile;groups=unique(fitDiagnostics.method,'stable');
    for j=1:numel(groups),u=fitDiagnostics.method==groups(j);plot(fitDiagnostics.foldID(u), ...
            fitDiagnostics.objectivePerObservation(u),'o-','LineWidth',1.3);hold on;end
    xlabel('CV fold');ylabel('NLL per training observation');grid on;legend(pretty(groups),'Location','best');
    title('Likelihood diagnostics');
    nexttile;for j=1:numel(groups),u=fitDiagnostics.method==groups(j);plot(fitDiagnostics.foldID(u), ...
            fitDiagnostics.noiseVariance(u),'o-','LineWidth',1.3);hold on;end
    set(gca,'YScale','log');xlabel('CV fold');ylabel('estimated noise variance');grid on;legend(pretty(groups),'Location','best');
    title('Nugget estimates');era5_save_figure(f,figureDir,'likelihood_and_nugget');

    f=figure('Visible',cfg.figures.visible,'Color','w','Position',[100 100 900 500]);
    elapsed=fitDiagnostics.fitSeconds+fitDiagnostics.predictionSeconds;
    for j=1:numel(groups)
        u=fitDiagnostics.method==groups(j);
        scatter(fitDiagnostics.M(u),elapsed(u),48,'filled');hold on;
    end
    set(gca,'YScale','log');legend(pretty(groups),'Location','best');
    xlabel('M=J^2 (exact shown at M=0)');ylabel('task time (s, log scale)');grid on;
    title('Observed computational cost by fold');
    era5_save_figure(f,figureDir,'runtime_by_feature_count');

    exactName="exact";hsgpName=sprintf('hsgp_J%02d',selectedJ);
    sites=unique(predictions.targetStation,'stable');
    meta=metrics(metrics.method==exactName,:);
    [~,corner]=min(meta.latitude+meta.longitude);[~,center]=min((meta.latitude-mean(meta.latitude)).^2+ ...
        (meta.longitude-mean(meta.longitude)).^2);[~,edge]=max(meta.longitude);
    chosen=unique([meta.targetStation(corner);meta.targetStation(center);meta.targetStation(edge)],'stable');
    f=figure('Visible',cfg.figures.visible,'Color','w','Position',[100 100 1250 720]);
    tiledlayout(numel(chosen),1,'TileSpacing','compact','Padding','compact');
    for s=1:numel(chosen)
        e=predictions(predictions.method==exactName & predictions.targetStation==chosen(s),:);
        h=predictions(predictions.method==hsgpName & predictions.targetStation==chosen(s),:);
        nexttile;plot(e.relativeHours/24,e.truthUWindMps,'k-','LineWidth',0.8);hold on;
        plot(e.relativeHours/24,e.predictionUWindMps,'-','LineWidth',1.1);
        if ~isempty(h),plot(h.relativeHours/24,h.predictionUWindMps,'--','LineWidth',1.1);end
        ylabel('u_{10} (m/s)');grid on;title(strrep(chosen(s),'_','\_'));
        if s==1,legend('ERA5 held out','exact Equation 9',sprintf('HSGP J=%d',selectedJ), ...
                'Location','best','NumColumns',3);end
    end
    xlabel('relative day');sgtitle('Representative held-site reconstructions');
    era5_save_figure(f,figureDir,'representative_reconstructions');

    selected=metrics(metrics.method==hsgpName,:);exact=metrics(metrics.method==exactName,:);
    [same,ia,ib]=intersect(exact.targetStation,selected.targetStation,'stable');
    f=figure('Visible',cfg.figures.visible,'Color','w','Position',[100 100 1000 450]);
    scatter(exact.longitude(ia),exact.latitude(ia),45,selected.rmse(ib)-exact.rmse(ia),'filled');
    colorbar;xlabel('longitude');ylabel('latitude');axis equal tight;grid on;
    title(sprintf('Held-site RMSE difference: HSGP J=%d minus exact',selectedJ));
    era5_save_figure(f,figureDir,'spatial_rmse_difference');
end

function out=pretty(values)
    out=strrep(string(values),'_',' ');
end
