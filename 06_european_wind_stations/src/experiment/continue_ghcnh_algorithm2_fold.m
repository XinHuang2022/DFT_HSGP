function continuation = continue_ghcnh_algorithm2_fold(foldID,maxCycles)
%CONTINUE_GHCNH_ALGORITHM2_FOLD Resume correction until ranks stabilize.

    if nargin<1, foldID=1; end
    if nargin<2, maxCycles=3; end
    C=ghcnh_trace_config();
    [cfg,subset,folds,A]=ghcnh_algorithm2_setup();
    foldDir=fullfile(C.resultsDir,sprintf('fold_%02d',foldID));
    source=fullfile(foldDir,'algorithm2_complete_task.mat');
    if ~isfile(source), error('Initial Algorithm-2 task is absent: %s',source); end
    loaded=load(source,'task'); previous=loaded.task;
    priorContinuationFile=fullfile(foldDir,'algorithm2_continuation.mat');
    if isfile(priorContinuationFile)
        prior=load(priorContinuationFile,'continuation');
        if isfield(prior.continuation,'fit') && ...
                isfield(prior.continuation,'selectedJ')
            previous.fit=prior.continuation.fit;
            previous.rankProfile=table((0:numel(prior.continuation.selectedJ)-1)', ...
                prior.continuation.selectedJ(:),prior.continuation.selectedJ(:).^2, ...
                'VariableNames',{'frequencyIndex','selectedJ','selectedM'});
        end
    end
    fold=folds(foldID); targetIndex=fold.targetIndex(:);
    observedIndex=setdiff((1:size(subset.U,2))',targetIndex,'stable');
    Yo=subset.U(:,observedIndex); Yt=subset.U(:,targetIndex);
    allCenter=mean(subset.xyKm,1);
    domainXy=(subset.xyKm-allCenter)/cfg.spatial.coordinateScaleKm;
    observed.xyCovariance=domainXy(observedIndex,:);
    target.xyCovariance=domainXy(targetIndex,:);
    meanCenter=mean(subset.xyKm(observedIndex,:),1);
    meanScale=std(subset.xyKm(observedIndex,:),0,1); meanScale(meanScale==0)=1;
    observed.design=A.makeDesign(size(Yo,1), ...
        (subset.xyKm(observedIndex,:)-meanCenter)./meanScale);
    target.design=A.makeDesign(size(Yo,1), ...
        (subset.xyKm(targetIndex,:)-meanCenter)./meanScale);
    master=adaptive_basis_bank(domainXy,cfg.hsgp.boundaryFactor, ...
        C.referenceJ,A,cfg);
    logFile=fullfile(foldDir,'progress.log');

    requiredFile=fullfile(foldDir,'final_parameter_required_rank_profile.csv');
    required=readtable(requiredFile);
    currentJ=max(previous.rankProfile.selectedJ,required.selectedJ);
    fit=previous.fit; rows=cell(maxCycles,1); stabilized=false;
    timerTotal=tic;
    for cycle=1:maxCycles
        timer=tic;
        ghcnh_log(logFile,sprintf('fold%d_continuation%d',foldID,cycle), ...
            sprintf('fit starting; J range %d--%d; reference frequencies %d', ...
            min(currentJ),max(currentJ),sum(currentJ==C.referenceJ)));
        fit=adaptive_hsgp_fit(Yo,observed.xyCovariance,observed.design, ...
            cfg,master,currentJ,fit.theta,cfg.hsgp.optimizer, ...
            cfg.seed.optimizer+10000*foldID+cycle,A);
        fitSeconds=toc(timer);
        timer=tic;
        [requiredJ,audit,checks]=trace_rank_audit(fit.theta, ...
            observed.xyCovariance,cfg,master,C,A,logFile, ...
            sprintf('fold%d_continuation%d_trace',foldID,cycle));
        auditSeconds=toc(timer);
        enlarged=max(currentJ,requiredJ); nEnlarged=sum(enlarged>currentJ);
        stabilized=nEnlarged==0;
        prefix=sprintf('continuation_%d',cycle);
        writetable(audit,fullfile(foldDir,[prefix '_trace_audit.csv']));
        writetable(checks,fullfile(foldDir,[prefix '_exact_spot_checks.csv']));
        rows{cycle}=table(cycle,fitSeconds,auditSeconds,fit.objective, ...
            fit.theta.noiseVariance,fit.convergenceCertified,min(currentJ), ...
            median(currentJ),max(currentJ),sum(currentJ==C.referenceJ), ...
            nEnlarged,stabilized,'VariableNames',{'cycle','fitSeconds', ...
            'auditSeconds','objective','noiseVariance','convergenceCertified', ...
            'minimumJ','medianJ','maximumJ','referenceFrequencies', ...
            'nFrequenciesEnlarged','rankStable'});
        ghcnh_log(logFile,sprintf('fold%d_continuation%d',foldID,cycle), ...
            sprintf('fit %.2f min; audit %.2fs; enlarged=%d; certified=%d', ...
            fitSeconds/60,auditSeconds,nEnlarged,fit.convergenceCertified));
        currentJ=enlarged;
        if stabilized, break; end
    end
    history=vertcat(rows{1:cycle});
    writetable(history,fullfile(foldDir,'continuation_history.csv'));
    profile=table((0:numel(currentJ)-1)',currentJ,currentJ.^2, ...
        'VariableNames',{'frequencyIndex','selectedJ','selectedM'});
    writetable(profile,fullfile(foldDir,'continuation_rank_profile.csv'));
    continuation=struct('foldID',foldID,'fit',fit,'selectedJ',currentJ, ...
        'history',history,'rankStable',stabilized,'elapsedSeconds',toc(timerTotal));
    if stabilized
        fit.selectedJ=currentJ;
        timer=tic;
        prediction=predict_adaptive_hsgp(Yo,observed.xyCovariance, ...
            target.xyCovariance,observed.design,target.design,fit,cfg,master,A);
        predictionSeconds=toc(timer);
        metrics=pooled_metrics(Yt,prediction,cfg,foldID);
        continuation.fit=fit; continuation.prediction=prediction;
        continuation.metrics=metrics;
        continuation.predictionSeconds=predictionSeconds;
        writetable(metrics,fullfile(foldDir,'continuation_metrics.csv'));
        ghcnh_log(logFile,sprintf('fold%d_continuation_complete',foldID), ...
            sprintf('rank stable; RMSE %.4f, MAE %.4f, NLPD %.4f', ...
            metrics.rmse,metrics.mae,metrics.nlpd));
    else
        ghcnh_log(logFile,sprintf('fold%d_continuation_incomplete',foldID), ...
            'rank profile did not stabilize; no final prediction reported');
    end
    save(fullfile(foldDir,'algorithm2_continuation.mat'),'continuation','-v7.3');
end

function metrics=pooled_metrics(truth,prediction,cfg,foldID)
    error=prediction.mean-truth;
    variance=max(prediction.observationVariance,1e-10);
    coverage=abs(error)<=norminv(0.5+cfg.cv.nominalCoverage/2).*sqrt(variance);
    correlation=zeros(size(truth,2),1);
    for j=1:size(truth,2)
        value=corrcoef(truth(:,j),prediction.mean(:,j));
        if all(size(value)==[2 2]), correlation(j)=value(1,2); else, correlation(j)=NaN; end
    end
    metrics=table(foldID,sqrt(mean(error.^2,'all')),mean(abs(error),'all'), ...
        mean(0.5*(log(2*pi*variance)+error.^2./variance),'all'), ...
        mean(coverage,'all'),mean(correlation,'omitnan'),'VariableNames', ...
        {'foldID','rmse','mae','nlpd','coverage95','meanSiteCorrelation'});
end
