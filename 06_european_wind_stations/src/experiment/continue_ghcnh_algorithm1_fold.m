function task=continue_ghcnh_algorithm1_fold(foldID,maxAdditionalCycles)
%CONTINUE_GHCNH_ALGORITHM1_FOLD Resume an unfinished exact-KL rank fit.

    if nargin<1, foldID=1; end
    if nargin<2, maxAdditionalCycles=3; end
    C=ghcnh_trace_config();
    [cfg,subset,folds,A]=ghcnh_algorithm2_setup();
    foldDir=fullfile(C.resultsDir,sprintf('fold_%02d',foldID),'algorithm1');
    source=fullfile(foldDir,'complete_task.mat');
    if ~isfile(source)
        error('Algorithm1:MissingCheckpoint','Missing checkpoint: %s',source);
    end
    loaded=load(source,'task'); previous=loaded.task;
    if previous.rankStable
        task=previous;
        return
    end

    targetIndex=folds(foldID).targetIndex(:);
    observedIndex=setdiff((1:size(subset.U,2))',targetIndex,'stable');
    Yo=subset.U(:,observedIndex); Yt=subset.U(:,targetIndex);
    [observed,target,domainXy]=fold_geometry(subset,observedIndex, ...
        targetIndex,cfg,A);
    master=adaptive_basis_bank(domainXy,cfg.hsgp.boundaryFactor, ...
        C.referenceJ,A,cfg);
    selectedJ=previous.rankProfile.selectedJ;
    fit=previous.fit;
    logFile=fullfile(foldDir,'progress.log'); totalTimer=tic;
    priorCycles=height(previous.history); rows=cell(maxAdditionalCycles,1);
    fitSecondsTotal=0; auditSecondsTotal=0; stable=false;

    for step=1:maxAdditionalCycles
        cycle=priorCycles+step; timer=tic;
        ghcnh_log(logFile,sprintf('algorithm1_fold%d_cycle%d',foldID,cycle), ...
            sprintf('continuation fit starting; J range %d--%d', ...
            min(selectedJ),max(selectedJ)));
        fit=adaptive_hsgp_fit(Yo,observed.xyCovariance,observed.design, ...
            cfg,master,selectedJ,fit.theta,cfg.hsgp.optimizer, ...
            cfg.seed.optimizer+20000*foldID+100*cycle,A);
        fitSeconds=toc(timer); fitSecondsTotal=fitSecondsTotal+fitSeconds;

        timer=tic;
        [requiredJ,audit,auditSummary]=exact_rank_audit(fit.theta, ...
            observed.xyCovariance,cfg,master,C,A,logFile, ...
            sprintf('algorithm1_fold%d_cycle%d_audit',foldID,cycle));
        auditSeconds=toc(timer); auditSecondsTotal=auditSecondsTotal+auditSeconds;
        enlarged=max(selectedJ,requiredJ); nEnlarged=sum(enlarged>selectedJ);
        stable=nEnlarged==0;
        write_audit(foldDir,sprintf('cycle_%d',cycle),requiredJ,audit,auditSummary);
        rows{step}=table(cycle,fitSeconds,auditSeconds,fit.objective, ...
            fit.theta.noiseVariance,fit.exitFlag,fit.convergenceCertified, ...
            fit.boundaryHit,min(selectedJ),median(selectedJ),max(selectedJ), ...
            nEnlarged,stable,'VariableNames',{'cycle','fitSeconds', ...
            'auditSeconds','objective','noiseVariance','exitFlag', ...
            'convergenceCertified','boundaryHit','minimumJ','medianJ', ...
            'maximumJ','nFrequenciesEnlarged','rankStable'});
        ghcnh_log(logFile,sprintf('algorithm1_fold%d_cycle%d',foldID,cycle), ...
            sprintf('continuation fit %.2f min; audit %.2fs; enlarged=%d; objective %.6g', ...
            fitSeconds/60,auditSeconds,nEnlarged,fit.objective));
        selectedJ=enlarged;
        if stable, break; end
    end

    history=[previous.history; vertcat(rows{1:step})];
    profile=table((0:numel(selectedJ)-1)',selectedJ,selectedJ.^2, ...
        'VariableNames',{'frequencyIndex','selectedJ','selectedM'});
    writetable(history,fullfile(foldDir,'iterations_continued.csv'));
    writetable(profile,fullfile(foldDir,'final_rank_profile_continued.csv'));

    if stable
        fit.selectedJ=selectedJ; timer=tic;
        prediction=predict_adaptive_hsgp(Yo,observed.xyCovariance, ...
            target.xyCovariance,observed.design,target.design,fit,cfg,master,A);
        predictionSeconds=toc(timer);
        [metrics,siteMetrics,predictions]=score_ghcnh_predictions( ...
            Yt,prediction,subset,targetIndex,foldID,cfg,'algorithm1_exact_kl');
        writetable(metrics,fullfile(foldDir,'metrics_continued.csv'));
        writetable(siteMetrics,fullfile(foldDir,'site_metrics_continued.csv'));
        writetable(predictions,fullfile(foldDir,'predictions_continued.csv'));
    else
        prediction=[]; predictionSeconds=0;
        metrics=table(); siteMetrics=table(); predictions=table();
    end

    continuationSeconds=toc(totalTimer);
    timing=previous.timing;
    timing.adaptiveFitSeconds=timing.adaptiveFitSeconds+fitSecondsTotal;
    timing.totalAuditSeconds=timing.totalAuditSeconds+auditSecondsTotal;
    timing.predictionSeconds=timing.predictionSeconds+predictionSeconds;
    timing.actualSeconds=timing.actualSeconds+continuationSeconds;
    timing.equivalentFirstUseSeconds=timing.equivalentFirstUseSeconds+continuationSeconds;
    timing.continuationSeconds=repmat(continuationSeconds,height(timing),1);
    parameters=struct2table(fit.theta,'AsArray',true);
    parameters.foldID=foldID; parameters.objective=fit.objective;
    parameters.exitFlag=fit.exitFlag;
    parameters.convergenceCertified=fit.convergenceCertified;
    parameters.boundaryHit=fit.boundaryHit;
    writetable(timing,fullfile(foldDir,'timing_continued.csv'));
    writetable(parameters,fullfile(foldDir,'parameters_continued.csv'));

    task=struct('foldID',foldID,'metrics',metrics,'siteMetrics',siteMetrics, ...
        'timing',timing,'parameters',parameters,'rankProfile',profile, ...
        'history',history,'fit',fit,'prediction',prediction, ...
        'predictions',predictions,'rankStable',stable, ...
        'observedIndex',observedIndex,'targetIndex',targetIndex, ...
        'version',C.version);
    save(fullfile(foldDir,'continuation_task.mat'),'task','C','cfg','-v7.3');
    if stable
        ghcnh_log(logFile,sprintf('algorithm1_fold%d_complete',foldID), ...
            sprintf('continued run stable: RMSE %.4f, objective %.6g, first-use %.2f min', ...
            metrics.rmse,fit.objective,timing.equivalentFirstUseSeconds/60));
    else
        ghcnh_log(logFile,sprintf('algorithm1_fold%d_incomplete',foldID), ...
            'continued rank profile did not stabilize; predictions withheld');
    end
end

function [observed,target,domainXy]=fold_geometry(subset,observedIndex, ...
        targetIndex,cfg,A)
    centerAll=mean(subset.xyKm,1);
    domainXy=(subset.xyKm-centerAll)/cfg.spatial.coordinateScaleKm;
    observed.xyCovariance=domainXy(observedIndex,:);
    target.xyCovariance=domainXy(targetIndex,:);
    center=mean(subset.xyKm(observedIndex,:),1);
    scale=std(subset.xyKm(observedIndex,:),0,1); scale(scale==0)=1;
    observed.design=A.makeDesign(size(subset.U,1), ...
        (subset.xyKm(observedIndex,:)-center)./scale);
    target.design=A.makeDesign(size(subset.U,1), ...
        (subset.xyKm(targetIndex,:)-center)./scale);
end

function write_audit(folder,prefix,selectedJ,audit,summary)
    writetable(audit,fullfile(folder,[prefix '_audit.csv']));
    writetable(summary,fullfile(folder,[prefix '_audit_timing.csv']));
    profile=table((0:numel(selectedJ)-1)',selectedJ,selectedJ.^2, ...
        'VariableNames',{'frequencyIndex','selectedJ','selectedM'});
    writetable(profile,fullfile(folder,[prefix '_rank_profile.csv']));
end
