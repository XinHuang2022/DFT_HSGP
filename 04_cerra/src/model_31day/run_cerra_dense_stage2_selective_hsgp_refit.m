function result=run_cerra_dense_stage2_selective_hsgp_refit(cfg)
%RUN_CERRA_DENSE_STAGE2_SELECTIVE_HSGP_REFIT One start for all M, then expand two.

    if nargin<1,cfg=cerra_northsea_config();end
    referenceFile=fullfile(cfg.paths.stage1Dir,'reference.mat');
    if ~isfile(referenceFile),error('Complete dense Stage 1 before Stage 2.');end
    ensure_directory(cfg.paths.stage2Dir);
    ensure_directory(fullfile(cfg.paths.stage2Dir,'tasks'));
    ensure_directory(fullfile(cfg.paths.stage2Dir,'best_models'));
    ensure_directory(fullfile(cfg.paths.stage2Dir,'figures'));
    Mvalues=cfg.hsgp.candidateM(:);
    phase1=cell(numel(Mvalues),1);
    [parallelUsed,poolWorkers]=cerra_open_pool(cfg,numel(Mvalues));
    if parallelUsed
        parfor q=1:numel(Mvalues)
            phase1{q}=cerra_dense_hsgp_start_task(cfg,referenceFile,Mvalues(q),1);
        end
    else
        for q=1:numel(Mvalues)
            phase1{q}=cerra_dense_hsgp_start_task(cfg,referenceFile,Mvalues(q),1);
        end
    end
    phase1Summary=vertcat(phase1{:});
    ranked=sortrows(phase1Summary,{'objectivePerObservation','M'}, ...
        {'ascend','ascend'});
    selectedM=ranked.M(1:min(cfg.stage2.nMultistartCandidates,height(ranked)));
    additionalStarts=(2:cfg.stage2.nIndependentStarts)';
    [MM,SS]=ndgrid(selectedM,additionalStarts);
    plan=table(MM(:),SS(:),'VariableNames',{'M','startID'});
    phase2=cell(height(plan),1);
    if parallelUsed
        parfor q=1:height(plan)
            phase2{q}=cerra_dense_hsgp_start_task(cfg,referenceFile, ...
                plan.M(q),plan.startID(q));
        end
    else
        for q=1:height(plan)
            phase2{q}=cerra_dense_hsgp_start_task(cfg,referenceFile, ...
                plan.M(q),plan.startID(q));
        end
    end
    if isempty(phase2),startSummary=phase1Summary;
    else,startSummary=[phase1Summary;vertcat(phase2{:})];end
    writetable(plan,fullfile(cfg.paths.stage2Dir,'additional_start_plan.csv'));
    writetable(startSummary,fullfile(cfg.paths.stage2Dir,'hsgp_start_summary.csv'));
    S=load(referenceFile,'reference');R=S.reference;
    modelRows=cell(numel(Mvalues),1);siteRows=cell(numel(Mvalues),1);
    bestFits=cell(numel(Mvalues),1);bestPredictions=cell(numel(Mvalues),1);
    for j=1:numel(Mvalues)
        M=Mvalues(j);use=startSummary.M==M;
        candidates=startSummary(use,:);
        [~,which]=min(candidates.objective);
        bestStart=candidates.startID(which);
        taskFile=fullfile(cfg.paths.stage2Dir,'tasks', ...
            sprintf('M_%04d_start_%02d.mat',M,bestStart));
        C=load(taskFile,'checkpoint');fit=C.checkpoint.fit;basis=C.checkpoint.basis;
        predictionTimer=tic;
        prediction=cerra_predict_fit("hsgp",R.Yfull,R.fullObserved, ...
            R.fullTarget,fit,R.fullCfg,basis);
        predictionSeconds=toc(predictionTimer);
        [row,siteRows{j}]=cerra_dense_prediction_summary( ...
            "hsgp_full_M"+string(M),M,numel(R.fullIndex),R.Ytarget, ...
            prediction,R.targetSiteTable,fit,R.fullCfg);
        row.bestStartID=bestStart;
        row.nCompletedStarts=sum(use);
        row.totalStartFitSeconds=sum(candidates.fitSeconds);
        row.predictionSeconds=predictionSeconds;
        row.maximumRadialFrequency=basis.maximumRadialFrequency;
        modelRows{j}=row;bestFits{j}=fit;bestPredictions{j}=prediction;
        save(fullfile(cfg.paths.stage2Dir,'best_models',sprintf('M_%04d.mat',M)), ...
            'fit','basis','prediction','bestStart','-v7.3');
    end
    modelSummary=vertcat(modelRows{:});bySite=vertcat(siteRows{:});
    writetable(modelSummary,fullfile(cfg.paths.stage2Dir,'hsgp_model_summary.csv'));
    writetable(bySite,fullfile(cfg.paths.stage2Dir,'hsgp_metrics_by_site.csv'));
    make_cerra_dense_stage2_figures(modelSummary,bestPredictions,R,cfg);
    runSummary=table(string(cfg.version),cfg.data.expectedTimePoints, ...
        numel(R.fullIndex),numel(R.targetIndex),strjoin(string(Mvalues),';'), ...
        strjoin(string(selectedM),';'),height(startSummary),parallelUsed,poolWorkers, ...
        'VariableNames',{'version','nTimePoints','nFullTrainingSites','nTestSites', ...
        'allM','multistartM','nCompletedStartTasks','parallelUsed','poolWorkers'});
    writetable(runSummary,fullfile(cfg.paths.stage2Dir,'run_summary.csv'));
    save(fullfile(cfg.paths.stage2Dir,'complete_result.mat'),'cfg','startSummary', ...
        'modelSummary','bySite','runSummary','selectedM','-v7.3');
    result=struct('modelSummary',modelSummary,'bySite',bySite, ...
        'startSummary',startSummary,'runSummary',runSummary);
end
