function make_cerra_dense_stage2_figures(summary,predictions,R,cfg)
%MAKE_CERRA_DENSE_STAGE2_FIGURES Accuracy/cost and three reconstructions.

    figure('Visible',cfg.figures.visible,'Color','w','Position',[100 100 1400 560]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile;plot(summary.M,summary.pooledRmse,'-o','LineWidth',1.8);hold on;
    yline(R.exactMetric.rmse,'--','exact n=250 reference');grid on;
    xlabel('M');ylabel('RMSE (m/s)');title('Refitted full-data HSGP');
    nexttile;plot(summary.M,summary.totalStartFitSeconds/3600,'-o','LineWidth',1.8);
    grid on;xlabel('M');ylabel('Total start CPU time (h)');title('Observed fitting cost');
    cerra_save_figure(gcf,fullfile(cfg.paths.stage2Dir,'figures', ...
        'hsgp_accuracy_and_runtime.png'));close(gcf);

    [~,best]=min(summary.objectivePerObservation);
    prediction=predictions{best};x=R.targetSiteTable.xKm;
    center=mean([min(x),max(x)]);[~,middle]=min(abs(x-center));
    [~,west]=min(x);[~,east]=max(x);chosen=unique([west,middle,east],'stable');
    figure('Visible',cfg.figures.visible,'Color','w','Position',[100 100 1500 850]);
    tiledlayout(numel(chosen),1,'Padding','compact','TileSpacing','compact');
    for q=1:numel(chosen)
        j=chosen(q);nexttile;
        plot(R.contract.relativeHours/24,R.Ytarget(:,j),'k-','LineWidth',1);hold on;
        plot(R.contract.relativeHours/24,prediction.mean(:,j),'-','LineWidth',1.2);
        grid on;ylabel('u_{10} (m/s)');title(string(R.targetSiteTable.stationID(j)));
        if q==1,legend('CERRA held out',char(summary.method(best)),'Location','best');end
    end
    xlabel('Relative day');
    cerra_save_figure(gcf,fullfile(cfg.paths.stage2Dir,'figures', ...
        'representative_reconstructions.png'));close(gcf);
end
