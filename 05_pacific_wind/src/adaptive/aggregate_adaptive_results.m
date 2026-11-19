function result=aggregate_adaptive_results(C,cfg,subset,folds,tasks,elapsedSeconds)
%AGGREGATE_ADAPTIVE_RESULTS Write manuscript-ready results and timing data.
    metrics=combine(tasks,'metrics'); predictions=combine(tasks,'predictions');
    rankProfile=combine(tasks,'rankProfile'); iterations=combine(tasks,'iterations');
    audits=combine(tasks,'audits'); references=combine(tasks,'references');
    timing=combine(tasks,'timing'); parameters=combine(tasks,'parameters');
    methodSummary=summarize_metrics(metrics,predictions,C);
    rankSummary=summarize_ranks(rankProfile);
    timingSummary=summarize_timing(timing);
    runSummary=table(C.version,C.dataset,height(subset.siteTable),numel(folds), ...
        elapsedSeconds,mean(timing.actualRunSeconds),sum(timing.actualRunSeconds), ...
        mean(timing.equivalentFirstUseSeconds), ...
        all(timing.profileStabilized),sum(timing.runStatus=="complete"), ...
        sum(contains(timing.runStatus,"not_stabilized")), ...
        'VariableNames',{'version','dataset', ...
        'nSites','nFolds','elapsedSeconds','meanActualFoldSeconds', ...
        'sumActualFoldSeconds','meanEquivalentFirstUseSeconds', ...
        'allProfilesStabilized','nCompletedFolds','nUnstableReferenceFolds'});
    writetable(metrics,fullfile(C.outputRoot,'prediction_metrics_by_site.csv'));
    writetable(predictions,fullfile(C.outputRoot,'predictions.csv'));
    writetable(rankProfile,fullfile(C.outputRoot,'frequency_rank_profiles.csv'));
    writetable(rankSummary,fullfile(C.outputRoot,'rank_summary_by_fold.csv'));
    writetable(iterations,fullfile(C.outputRoot,'adaptive_iterations.csv'));
    writetable(audits,fullfile(C.outputRoot,'rank_audits.csv'));
    writetable(references,fullfile(C.outputRoot,'reference_stability.csv'));
    writetable(timing,fullfile(C.outputRoot,'timing_by_fold.csv'));
    writetable(timingSummary,fullfile(C.outputRoot,'timing_summary.csv'));
    writetable(parameters,fullfile(C.outputRoot,'parameter_estimates.csv'));
    writetable(methodSummary,fullfile(C.outputRoot,'method_summary.csv'));
    writetable(runSummary,fullfile(C.outputRoot,'run_summary.csv'));
    result=struct('C',C,'metrics',metrics,'predictions',predictions, ...
        'rankProfile',rankProfile,'rankSummary',rankSummary, ...
        'iterations',iterations,'audits',audits,'references',references, ...
        'timing',timing,'timingSummary',timingSummary, ...
        'parameters',parameters,'methodSummary',methodSummary, ...
        'runSummary',runSummary);
    save(fullfile(C.outputRoot,'complete_result.mat'),'result','cfg','-v7.3');
end

function value=combine(tasks,name)
    cells=cellfun(@(x)x.(name),tasks,'UniformOutput',false); value=vertcat(cells{:});
end

function summary=summarize_metrics(metrics,predictions,C)
    if isempty(metrics)
        summary=table("adaptive_hsgp",0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
            'VariableNames',{'method','nTargetSites','meanSiteRmse', ...
            'bootstrapLower95','bootstrapUpper95','pooledRmse','meanMae', ...
            'meanCorrelation','meanNlpd','meanCoverage95'});
        return;
    end
    siteRmse=metrics.rmse; foldId=metrics.foldID;
    rng(C.bootstrapSeed,'twister'); B=C.bootstrapReplicates;
    levels=unique(foldId); boot=zeros(B,1);
    for b=1:B
        sampled=levels(randi(numel(levels),numel(levels),1)); values=[];
        for q=1:numel(sampled),values=[values;siteRmse(foldId==sampled(q))];end %#ok<AGROW>
        boot(b)=mean(values);
    end
    ci=prctile(boot,[2.5 97.5]);
    pooled=sqrt(mean((predictions.truthUWindMps-predictions.predictionUWindMps).^2));
    summary=table("adaptive_hsgp",height(metrics),mean(siteRmse),ci(1),ci(2), ...
        pooled,mean(metrics.mae),mean(metrics.correlation,'omitnan'), ...
        mean(metrics.nlpd,'omitnan'),mean(metrics.coverage95,'omitnan'), ...
        'VariableNames',{'method','nTargetSites','meanSiteRmse', ...
        'bootstrapLower95','bootstrapUpper95','pooledRmse','meanMae', ...
        'meanCorrelation','meanNlpd','meanCoverage95'});
end

function summary=summarize_ranks(x)
    fold=unique(x.foldID); rows=cell(numel(fold),1);
    for j=1:numel(fold)
        v=x.selectedM(x.foldID==fold(j));
        rows{j}=table(fold(j),min(v),median(v),mean(v),max(v), ...
            numel(unique(v)),'VariableNames',{'foldID','minimumM', ...
            'medianM','meanM','maximumM','nDistinctRanks'});
    end
    summary=vertcat(rows{:});
end

function summary=summarize_timing(x)
    names={'basisSeconds','pilotFitSeconds','initialAuditSeconds', ...
        'adaptiveFitSeconds','verificationSeconds','finalCorrectionSeconds', ...
        'predictionSeconds','actualRunSeconds','equivalentFirstUseSeconds'};
    stage=string(names(:)); meanSeconds=zeros(numel(names),1); ...
        sumSeconds=zeros(numel(names),1); maximumSeconds=zeros(numel(names),1);
    for j=1:numel(names)
        v=x.(names{j}); meanSeconds(j)=mean(v); sumSeconds(j)=sum(v); ...
            maximumSeconds(j)=max(v);
    end
    summary=table(stage,meanSeconds,sumSeconds,maximumSeconds);
end
