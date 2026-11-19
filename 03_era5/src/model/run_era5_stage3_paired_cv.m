function result=run_era5_stage3_paired_cv(cfg)
%RUN_ERA5_STAGE3_PAIRED_CV Exact and selected-J HSGP five-fold CV.
    if nargin<1,cfg=era5_northsea_config();end
    auditFile=fullfile(cfg.paths.stage1Dir,'hsgp_J_audit_summary.csv');
    if ~isfile(auditFile),error('Complete Stage 1 before final CV.');end
    audit=readtable(auditFile);candidate=cfg.hsgp.candidateJ(:);
    selected=cfg.final.selectedJ;
    if isempty(selected)
        passMask=era5_logical_column(audit.passesAudit);
        pass=audit.J(passMask);
        if isempty(pass),selected=audit.J(end);else,selected=min(pass);end
    end
    idx=find(candidate==selected,1);if isempty(idx),error('selectedJ must be on candidateJ ladder.');end
    Jvalues=selected;
    if cfg.final.includeAdjacentJ
        Jvalues=unique(candidate(max(1,idx-1):min(numel(candidate),idx+1)),'stable');
    end
    ensure_directory(cfg.paths.stage3Dir);ensure_directory(fullfile(cfg.paths.stage3Dir,'tasks'));
    ensure_directory(fullfile(cfg.paths.stage3Dir,'figures'));
    contract=load_era5_wind_data(cfg);subset=make_era5_subset(contract,'full17');
    [folds,foldManifest,foldSummary]=make_era5_cv_folds(subset,cfg.cv.holdoutPattern);
    writetable(subset.siteTable,fullfile(cfg.paths.stage3Dir,'selected_subset_manifest.csv'));
    writetable(foldManifest,fullfile(cfg.paths.stage3Dir,'cv_fold_manifest.csv'));
    writetable(foldSummary,fullfile(cfg.paths.stage3Dir,'cv_fold_summary.csv'));
    nTrain=height(subset.siteTable)-height(subset.siteTable)/cfg.cv.nFolds;
    methods=["exact";repmat("hsgp",numel(Jvalues),1)];Js=[0;Jvalues(:)];
    rows=cell(numel(methods)*numel(folds),1);c=0;
    for m=1:numel(methods)
        for f=1:numel(folds)
            c=c+1;M=Js(m)^2;
            if methods(m)=="exact",cost=nTrain^3;elseif M<nTrain,cost=M^3+nTrain*M^2;
            else,cost=nTrain^3+nTrain^2*M;end
            rows{c}=table(methods(m),Js(m),M,folds(f).foldID,cost, ...
                'VariableNames',{'engine','J','M','foldID','estimatedCostScore'});
        end
    end
    plan=vertcat(rows{:});plan=sortrows(plan,'estimatedCostScore','descend');
    plan.queuePosition=(1:height(plan))';
    writetable(plan,fullfile(cfg.paths.stage3Dir,'task_plan.csv'));
    fullDomainXy=(subset.xyKm-mean(subset.xyKm,1))./cfg.spatial.coordinateScaleKm;
    tasks=cell(height(plan),1);[parallelUsed,poolWorkers]=era5_open_pool(cfg,height(plan));
    if parallelUsed
        parfor q=1:height(plan)
            fold=folds([folds.foldID]==plan.foldID(q));
            tasks{q}=era5_cv_task(cfg,subset,fold,plan.engine(q),plan.J(q),fullDomainXy);
        end
    else
        for q=1:height(plan)
            fold=folds([folds.foldID]==plan.foldID(q));
            tasks{q}=era5_cv_task(cfg,subset,fold,plan.engine(q),plan.J(q),fullDomainXy);
        end
    end
    metrics=combine(tasks,'metrics');predictions=combine(tasks,'predictions');
    fitDiagnostics=combine(tasks,'fitDiagnostics');checks=combine(tasks,'checks');
    timing=combine(tasks,'timing');methodSummary=summarize_methods(metrics,cfg);
    parameterSummary=summarize_parameters(fitDiagnostics);
    likelihoodSummary=summarize_likelihood(fitDiagnostics);
    writetable(metrics,fullfile(cfg.paths.stage3Dir,'prediction_metrics_by_site.csv'));
    writetable(predictions,fullfile(cfg.paths.stage3Dir,'predictions.csv'));
    writetable(fitDiagnostics,fullfile(cfg.paths.stage3Dir,'fit_diagnostics.csv'));
    writetable(checks,fullfile(cfg.paths.stage3Dir,'implementation_checks.csv'));
    writetable(timing,fullfile(cfg.paths.stage3Dir,'task_timing.csv'));
    writetable(methodSummary,fullfile(cfg.paths.stage3Dir,'method_summary.csv'));
    writetable(parameterSummary,fullfile(cfg.paths.stage3Dir,'parameter_summary.csv'));
    writetable(likelihoodSummary,fullfile(cfg.paths.stage3Dir,'likelihood_summary.csv'));
    runSummary=table(string(cfg.version),height(subset.siteTable),numel(folds), ...
        selected,strjoin(string(Jvalues),';'),height(plan),parallelUsed,poolWorkers, ...
        all(checks.passed),'VariableNames',{'version','nSites','nFolds', ...
        'selectedJ','evaluatedJ','nTasks','parallelUsed','poolWorkers', ...
        'allImplementationChecksPassed'});
    writetable(runSummary,fullfile(cfg.paths.stage3Dir,'run_summary.csv'));
    make_era5_stage3_figures(metrics,predictions,fitDiagnostics,methodSummary,selected,cfg);
    save(fullfile(cfg.paths.stage3Dir,'complete_result.mat'),'cfg','subset','folds', ...
        'plan','metrics','predictions','fitDiagnostics','checks','timing', ...
        'methodSummary','parameterSummary','likelihoodSummary','runSummary','-v7.3');
    result=struct('metrics',metrics,'methodSummary',methodSummary, ...
        'fitDiagnostics',fitDiagnostics,'runSummary',runSummary);
end

function value=combine(tasks,name)
    parts=cellfun(@(x)x.(name),tasks,'UniformOutput',false);value=vertcat(parts{:});
end
function out=summarize_methods(x,cfg)
    methods=unique(x.method,'stable');rows=cell(numel(methods),1);rng(90477,'twister');
    for j=1:numel(methods)
        use=x.method==methods(j);v=x.rmse(use);foldId=x.foldID(use);
        foldLevels=unique(foldId);B=2000;boot=zeros(B,1);nFold=numel(foldLevels);
        for b=1:B
            selectedFolds=foldLevels(randi(nFold,nFold,1));sample=[];
            for q=1:nFold,sample=[sample;v(foldId==selectedFolds(q))];end %#ok<AGROW>
            boot(b)=mean(sample);
        end
        boot=sort(boot);ci=boot([max(1,ceil(0.025*B)),min(B,ceil(0.975*B))]);
        rows{j}=table(methods(j),sum(use),mean(v),ci(1),ci(2), ...
            mean(x.mae(use)),mean(x.correlation(use),'omitnan'), ...
            mean(x.nlpd(use),'omitnan'),mean(x.coverage95(use),'omitnan'), ...
            'VariableNames',{'method','nSites','meanRmse','bootstrapLower95', ...
            'bootstrapUpper95','meanMae','meanCorrelation','meanNlpd','meanCoverage95'});
    end
    out=vertcat(rows{:});
end
function out=summarize_parameters(x)
    methods=unique(x.method,'stable');rows=cell(numel(methods),1);
    for j=1:numel(methods),u=x.method==methods(j);rows{j}=table(methods(j),sum(u), ...
        mean(x.signalVariance(u)),mean(x.alphaPer100Km(u)),mean(x.betaPerSixHours(u)), ...
        mean(x.epsilon(u)),mean(x.noiseVariance(u)), ...
        'VariableNames',{'method','nFolds','meanSignalVariance','meanAlphaPer100Km', ...
        'meanBetaPerSixHours','meanEpsilon','meanNoiseVariance'});end
    out=vertcat(rows{:});
end
function out=summarize_likelihood(x)
    methods=unique(x.method,'stable');rows=cell(numel(methods),1);
    for j=1:numel(methods),u=x.method==methods(j);rows{j}=table(methods(j),sum(u), ...
        mean(x.objectivePerObservation(u)),std(x.objectivePerObservation(u)), ...
        mean(x.convergenceCertified(u)),mean(x.boundaryHit(u)), ...
        'VariableNames',{'method','nFolds','meanNllPerObservation', ...
        'sdNllPerObservation','convergenceFraction','boundaryHitFraction'});end
    out=vertcat(rows{:});
end
