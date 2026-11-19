function task=run_ghcnh_algorithm1_fold(foldID,maxCycles)
%RUN_GHCNH_ALGORITHM1_FOLD Exact-KL adaptive-rank fit and reconstruction.

    if nargin<1, foldID=1; end
    if nargin<2, maxCycles=3; end
    C=ghcnh_trace_config();
    [cfg,subset,folds,A]=ghcnh_algorithm2_setup();
    fold=folds(foldID); targetIndex=fold.targetIndex(:);
    observedIndex=setdiff((1:size(subset.U,2))',targetIndex,'stable');
    Yo=subset.U(:,observedIndex); Yt=subset.U(:,targetIndex);
    [observed,target,domainXy]=fold_geometry(subset,observedIndex, ...
        targetIndex,cfg,A);
    foldDir=fullfile(C.resultsDir,sprintf('fold_%02d',foldID), ...
        'algorithm1');
    if ~isfolder(foldDir), mkdir(foldDir); end
    logFile=fullfile(foldDir,'progress.log'); totalTimer=tic;
    ghcnh_log(logFile,sprintf('algorithm1_fold%d',foldID),sprintf( ...
        'started with %d training and %d test stations', ...
        numel(observedIndex),numel(targetIndex)));

    timer=tic;
    master=adaptive_basis_bank(domainXy,cfg.hsgp.boundaryFactor, ...
        C.referenceJ,A,cfg);
    basisSeconds=toc(timer);
    pilotSource=fullfile(C.resultsDir,sprintf('fold_%02d',foldID), ...
        'algorithm2_pilot_fit.mat');
    if foldID==1 && isfile(pilotSource)
        loaded=load(pilotSource,'pilotFit','pilotSeconds');
        pilotFit=loaded.pilotFit; pilotSeconds=loaded.pilotSeconds;
        pilotReused=true;
        ghcnh_log(logFile,sprintf('algorithm1_fold%d',foldID), ...
            'loaded the identical Algorithm-2 pilot checkpoint');
    else
        pilotBasis=subset_adaptive_basis(master,C.initialJ);
        timer=tic;
        pilotFit=A.fitFixed(Yo,observed.xyCovariance,observed.design, ...
            cfg,pilotBasis,cfg.start,[],[],cfg.hsgp.optimizer, ...
            cfg.seed.optimizer+1000*foldID);
        pilotSeconds=toc(timer); pilotReused=false;
    end

    timer=tic;
    [selectedJ,initialAudit,initialSummary]=exact_rank_audit( ...
        pilotFit.theta,observed.xyCovariance,cfg,master,C,A,logFile, ...
        sprintf('algorithm1_fold%d_initial',foldID));
    initialAuditSeconds=toc(timer);
    write_audit(foldDir,'initial',selectedJ,initialAudit,initialSummary);

    fit=pilotFit; rows=cell(maxCycles,1); fitSecondsTotal=0;
    auditSecondsTotal=initialAuditSeconds; stable=false;
    for cycle=1:maxCycles
        timer=tic;
        ghcnh_log(logFile,sprintf('algorithm1_fold%d_cycle%d',foldID,cycle), ...
            sprintf('fit starting; J range %d--%d',min(selectedJ),max(selectedJ)));
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
        rows{cycle}=table(cycle,fitSeconds,auditSeconds,fit.objective, ...
            fit.theta.noiseVariance,fit.exitFlag,fit.convergenceCertified, ...
            fit.boundaryHit,min(selectedJ),median(selectedJ),max(selectedJ), ...
            nEnlarged,stable,'VariableNames',{'cycle','fitSeconds', ...
            'auditSeconds','objective','noiseVariance','exitFlag', ...
            'convergenceCertified','boundaryHit','minimumJ','medianJ', ...
            'maximumJ','nFrequenciesEnlarged','rankStable'});
        ghcnh_log(logFile,sprintf('algorithm1_fold%d_cycle%d',foldID,cycle), ...
            sprintf('fit %.2f min; audit %.2fs; enlarged=%d; objective %.6g', ...
            fitSeconds/60,auditSeconds,nEnlarged,fit.objective));
        selectedJ=enlarged;
        if stable, break; end
    end
    history=vertcat(rows{1:cycle});
    profile=table((0:numel(selectedJ)-1)',selectedJ,selectedJ.^2, ...
        'VariableNames',{'frequencyIndex','selectedJ','selectedM'});
    writetable(history,fullfile(foldDir,'iterations.csv'));
    writetable(profile,fullfile(foldDir,'final_rank_profile.csv'));

    if stable
        fit.selectedJ=selectedJ;
        timer=tic;
        prediction=predict_adaptive_hsgp(Yo,observed.xyCovariance, ...
            target.xyCovariance,observed.design,target.design,fit,cfg,master,A);
        predictionSeconds=toc(timer);
        [metrics,siteMetrics,predictions]=score_ghcnh_predictions( ...
            Yt,prediction,subset,targetIndex,foldID,cfg,'algorithm1_exact_kl');
        writetable(metrics,fullfile(foldDir,'metrics.csv'));
        writetable(siteMetrics,fullfile(foldDir,'site_metrics.csv'));
        writetable(predictions,fullfile(foldDir,'predictions.csv'));
    else
        prediction=[]; predictionSeconds=0;
        metrics=table(); siteMetrics=table(); predictions=table();
    end
    actualSeconds=toc(totalTimer);
    equivalentFirstUseSeconds=actualSeconds+pilotReused*pilotSeconds;
    timing=table(foldID,basisSeconds,pilotSeconds,initialAuditSeconds, ...
        fitSecondsTotal,auditSecondsTotal,predictionSeconds,actualSeconds, ...
        equivalentFirstUseSeconds,pilotReused,'VariableNames',{'foldID', ...
        'basisSeconds','pilotFitSeconds','initialAuditSeconds', ...
        'adaptiveFitSeconds','totalAuditSeconds','predictionSeconds', ...
        'actualSeconds','equivalentFirstUseSeconds','pilotReused'});
    parameters=struct2table(fit.theta,'AsArray',true);
    parameters.foldID=foldID; parameters.objective=fit.objective;
    parameters.exitFlag=fit.exitFlag;
    parameters.convergenceCertified=fit.convergenceCertified;
    parameters.boundaryHit=fit.boundaryHit;
    writetable(timing,fullfile(foldDir,'timing.csv'));
    writetable(parameters,fullfile(foldDir,'parameters.csv'));
    task=struct('foldID',foldID,'metrics',metrics,'siteMetrics',siteMetrics, ...
        'timing',timing,'parameters',parameters,'rankProfile',profile, ...
        'history',history,'fit',fit,'prediction',prediction, ...
        'predictions',predictions,'rankStable',stable, ...
        'observedIndex',observedIndex,'targetIndex',targetIndex, ...
        'version',C.version);
    save(fullfile(foldDir,'complete_task.mat'),'task','C','cfg','-v7.3');
    if stable
        ghcnh_log(logFile,sprintf('algorithm1_fold%d_complete',foldID), ...
            sprintf('RMSE %.4f, objective %.6g, first-use %.2f min', ...
            metrics.rmse,fit.objective,equivalentFirstUseSeconds/60));
    else
        ghcnh_log(logFile,sprintf('algorithm1_fold%d_incomplete',foldID), ...
            'rank profile did not stabilize; predictions withheld');
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
