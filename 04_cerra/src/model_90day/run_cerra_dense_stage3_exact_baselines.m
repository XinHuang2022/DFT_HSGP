function result=run_cerra_dense_stage3_exact_baselines(cfg)
%RUN_CERRA_DENSE_STAGE3_EXACT_BASELINES Coarse exact fits and final frontier.

    if nargin<1,cfg=cerra_northsea_config();end
    referenceFile=fullfile(cfg.paths.stage1Dir,'reference.mat');
    hsgpSummaryFile=fullfile(cfg.paths.stage2Dir,'hsgp_model_summary.csv');
    hsgpSiteFile=fullfile(cfg.paths.stage2Dir,'hsgp_metrics_by_site.csv');
    if ~isfile(referenceFile)||~isfile(hsgpSummaryFile)||~isfile(hsgpSiteFile)
        error('Complete dense Stages 1 and 2 before Stage 3.');
    end
    ensure_directory(cfg.paths.stage3Dir);
    ensure_directory(fullfile(cfg.paths.stage3Dir,'tasks'));
    ensure_directory(fullfile(cfg.paths.stage3Dir,'best_models'));
    ensure_directory(fullfile(cfg.paths.stage3Dir,'figures'));
    sizes=cfg.stage3.exactSizes(:);starts=(1:cfg.stage3.nIndependentStarts)';
    [NN,SS]=ndgrid(sizes,starts);
    plan=table(NN(:),SS(:),'VariableNames',{'nCoarseSites','startID'});
    writetable(plan,fullfile(cfg.paths.stage3Dir,'task_plan.csv'));
    tasks=cell(height(plan),1);
    [parallelUsed,poolWorkers]=cerra_open_pool(cfg,height(plan));
    if parallelUsed
        parfor q=1:height(plan)
            tasks{q}=cerra_dense_exact_start_task(cfg,referenceFile, ...
                plan.nCoarseSites(q),plan.startID(q));
        end
    else
        for q=1:height(plan)
            tasks{q}=cerra_dense_exact_start_task(cfg,referenceFile, ...
                plan.nCoarseSites(q),plan.startID(q));
        end
    end
    startSummary=vertcat(tasks{:});
    writetable(startSummary,fullfile(cfg.paths.stage3Dir,'exact_start_summary.csv'));
    S=load(referenceFile,'reference');R=S.reference;
    modelRows=cell(numel(sizes),1);siteRows=cell(numel(sizes),1);
    for j=1:numel(sizes)
        nCoarse=sizes(j);use=startSummary.nCoarseSites==nCoarse;
        candidates=startSummary(use,:);[~,which]=min(candidates.objective);
        bestStart=candidates.startID(which);
        C=load(fullfile(cfg.paths.stage3Dir,'tasks', ...
            sprintf('n_%04d_start_%02d.mat',nCoarse,bestStart)),'checkpoint');
        observedIndex=C.checkpoint.observedIndex;fitCfg=C.checkpoint.fitCfg;
        fit=C.checkpoint.fit;
        [observed,target]=cerra_fold_geometry(R.contract,observedIndex, ...
            R.targetIndex,cfg);
        Y=R.contract.U(:,observedIndex);
        predictionTimer=tic;
        prediction=cerra_predict_fit("exact",Y,observed,target,fit,fitCfg,[]);
        predictionSeconds=toc(predictionTimer);
        [row,siteRows{j}]=cerra_dense_prediction_summary( ...
            "exact_n"+string(nCoarse),0,nCoarse,R.Ytarget,prediction, ...
            R.targetSiteTable,fit,fitCfg);
        row.bestStartID=bestStart;row.nCompletedStarts=sum(use);
        row.totalStartFitSeconds=sum(candidates.fitSeconds);
        row.predictionSeconds=predictionSeconds;modelRows{j}=row;
        save(fullfile(cfg.paths.stage3Dir,'best_models', ...
            sprintf('n_%04d.mat',nCoarse)),'fit','prediction','observedIndex', ...
            'bestStart','-v7.3');
    end
    exactSummary=vertcat(modelRows{:});exactBySite=vertcat(siteRows{:});
    writetable(exactSummary,fullfile(cfg.paths.stage3Dir,'exact_model_summary.csv'));
    writetable(exactBySite,fullfile(cfg.paths.stage3Dir,'exact_metrics_by_site.csv'));
    hsgpSummary=readtable(hsgpSummaryFile,'TextType','string');
    hsgpBySite=readtable(hsgpSiteFile,'TextType','string');
    combinedSummary=[hsgpSummary;align_exact(exactSummary,hsgpSummary)];
    writetable(combinedSummary,fullfile(cfg.paths.stage3Dir,'accuracy_cost_frontier.csv'));
    matched=matched_comparisons(hsgpBySite,exactBySite,sizes,cfg);
    writetable(matched,fullfile(cfg.paths.stage3Dir,'matched_comparisons.csv'));
    make_cerra_dense_stage3_figures(combinedSummary,matched,cfg);
    runSummary=table(string(cfg.version),cfg.data.expectedTimePoints, ...
        numel(R.fullIndex),numel(R.targetIndex),strjoin(string(sizes),';'), ...
        height(plan),parallelUsed,poolWorkers, ...
        'VariableNames',{'version','nTimePoints','nFullTrainingSites', ...
        'nTestSites','exactSizes','nStartTasks','parallelUsed','poolWorkers'});
    writetable(runSummary,fullfile(cfg.paths.stage3Dir,'run_summary.csv'));
    save(fullfile(cfg.paths.stage3Dir,'complete_result.mat'),'cfg','startSummary', ...
        'exactSummary','exactBySite','combinedSummary','matched','runSummary','-v7.3');
    result=struct('exactSummary',exactSummary,'combinedSummary',combinedSummary, ...
        'matched',matched,'runSummary',runSummary);
end

function exactAligned=align_exact(exactSummary,hsgpSummary)
    missing=setdiff(hsgpSummary.Properties.VariableNames, ...
        exactSummary.Properties.VariableNames,'stable');
    for j=1:numel(missing)
        name=missing{j};template=hsgpSummary.(name);
        if isstring(template),exactSummary.(name)=repmat("",height(exactSummary),1);
        elseif islogical(template),exactSummary.(name)=false(height(exactSummary),1);
        else,exactSummary.(name)=nan(height(exactSummary),1);end
    end
    exactAligned=exactSummary(:,hsgpSummary.Properties.VariableNames);
end

function out=matched_comparisons(hsgp,exact,sizes,cfg)
    available=intersect(unique(hsgp.M),sizes);rows=cell(numel(available),1);
    for j=1:numel(available)
        M=available(j);a=hsgp(hsgp.M==M,:);
        b=exact(exact.nTrainingSites==M,:);
        boot=cerra_spatial_block_bootstrap_difference(a,b,5,2000, ...
            cfg.seed.split+M);
        rows{j}=table(M,boot.nSites,boot.nSpatialBlocks,boot.meanDifference, ...
            boot.lower95,boot.upper95, ...
            'VariableNames',{'matchedDimension','nTestSites','nSpatialBlocks', ...
            'meanSiteRmseDifferenceHsgpMinusExact','lower95','upper95'});
    end
    if isempty(rows)
        out=table([],[],[],[],[],[],'VariableNames', ...
            {'matchedDimension','nTestSites','nSpatialBlocks', ...
            'meanSiteRmseDifferenceHsgpMinusExact','lower95','upper95'});
    else,out=vertcat(rows{:});end
end
