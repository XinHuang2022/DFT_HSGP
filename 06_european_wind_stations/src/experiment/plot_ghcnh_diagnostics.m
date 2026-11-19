function files = plot_ghcnh_diagnostics()
%PLOT_GHCNH_DIAGNOSTICS Render data and Algorithm-2 pilot diagnostics.

    C=ghcnh_trace_config();
    if ~isfolder(fullfile(C.resultsDir,'figures'))
        mkdir(fullfile(C.resultsDir,'figures'));
    end
    S=load(fullfile(C.dataProcessedDir,'wind_panel_1200.mat'), ...
        'coverage','Y');
    P=load(fullfile(C.dataProcessedDir,'ghcnh_complete_panel.mat'),'subset');
    folds=make_ghcnh_cv_folds(P.subset,C.nFolds);
    foldID=zeros(height(P.subset.siteTable),1);
    for f=1:numel(folds), foldID(folds(f).targetIndex)=f; end

    fig=figure('Visible','off','Color','w','Position',[100 100 1100 760]);
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
    nexttile; hold on;
    scatter(S.coverage.longitude(~S.coverage.nearComplete), ...
        S.coverage.latitude(~S.coverage.nearComplete),8,[0.75 0.75 0.75],'filled');
    scatter(S.coverage.longitude(S.coverage.nearComplete & ~S.coverage.complete), ...
        S.coverage.latitude(S.coverage.nearComplete & ~S.coverage.complete), ...
        10,[0.90 0.55 0.10],'filled');
    scatter(S.coverage.longitude(S.coverage.complete), ...
        S.coverage.latitude(S.coverage.complete),12,[0.05 0.35 0.65],'filled');
    xlabel('Longitude'); ylabel('Latitude'); title('(a) Completeness screening');
    legend({'Other','At most three gaps','Complete'},'Location','best'); grid on; box on;

    nexttile;
    scatter(P.subset.siteTable.longitude,P.subset.siteTable.latitude,20,foldID,'filled');
    xlabel('Longitude'); ylabel('Latitude'); title('(b) Complete-panel folds');
    cb=colorbar; cb.Label.String='Fold'; colormap(gca,lines(C.nFolds)); grid on; box on;

    nexttile;
    histogram(S.coverage.nMissing,'BinEdges',[-0.5:1:10.5,20.5:10:100.5,360.5], ...
        'FaceColor',[0.25 0.5 0.75]);
    xlabel('Missing target times'); ylabel('Stations'); title('(c) Record completeness');
    grid on; box on;

    nexttile;
    values=P.subset.U(:);
    histogram(values,40,'Normalization','pdf','FaceColor',[0.2 0.55 0.35]);
    xlabel('Wind speed (m/s)'); ylabel('Density'); title('(d) Complete-panel wind speed');
    grid on; box on;
    dataPdf=fullfile(C.resultsDir,'figures','ghcnh_data_diagnostics.pdf');
    dataPng=fullfile(C.resultsDir,'figures','ghcnh_data_diagnostics.png');
    dataFig=fullfile(C.resultsDir,'figures','ghcnh_data_diagnostics.fig');
    exportgraphics(fig,dataPdf,'ContentType','vector');
    exportgraphics(fig,dataPng,'Resolution',220); savefig(fig,dataFig); close(fig);

    continuation=load(fullfile(C.resultsDir,'fold_01','algorithm2_continuation.mat'), ...
        'continuation');
    rank=continuation.continuation.selectedJ;
    checks=readtable(fullfile(C.resultsDir,'fold_01', ...
        'continuation_1_exact_spot_checks.csv'));
    fig=figure('Visible','off','Color','w','Position',[100 100 1050 390]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile;
    stairs(0:numel(rank)-1,rank,'LineWidth',1.6,'Color',[0.05 0.35 0.65]);
    yline(C.referenceJ,'--','Reference rank','LineWidth',1.2);
    xlabel('Temporal frequency index k'); ylabel('Selected J_k');
    title('(a) Trace-certified rank profile'); grid on; box on;
    nexttile;
    semilogy(checks.frequencyIndex,max(checks.exactKlPerSite,1e-12), ...
        'o-','LineWidth',1.5,'MarkerSize',6); hold on;
    semilogy(checks.frequencyIndex,max(checks.traceBoundPerSite,1e-12), ...
        's--','LineWidth',1.5,'MarkerSize',6);
    yline(C.epsilonKlPerSite,':','Tolerance','LineWidth',1.2);
    xlabel('Temporal frequency index k'); ylabel('Divergence per site');
    title('(b) Exact KL and trace certificate');
    legend({'Exact KL','Trace bound','Tolerance'},'Location','best'); grid on; box on;
    rankPdf=fullfile(C.resultsDir,'figures','ghcnh_algorithm2_diagnostics.pdf');
    rankPng=fullfile(C.resultsDir,'figures','ghcnh_algorithm2_diagnostics.png');
    rankFig=fullfile(C.resultsDir,'figures','ghcnh_algorithm2_diagnostics.fig');
    exportgraphics(fig,rankPdf,'ContentType','vector');
    exportgraphics(fig,rankPng,'Resolution',220); savefig(fig,rankFig); close(fig);
    files=string({dataPdf;dataPng;dataFig;rankPdf;rankPng;rankFig});
    ghcnh_log(C.progressLog,'figures','saved data and Algorithm-2 diagnostics');
end
