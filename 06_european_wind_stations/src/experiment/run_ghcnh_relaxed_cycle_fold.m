function task=run_ghcnh_relaxed_cycle_fold(method,foldID,maxCycles,maxEnlargedFraction)
%RUN_GHCNH_RELAXED_CYCLE_FOLD Fold-level early-stopping sensitivity.
%
% This diagnostic keeps every data, model, optimizer, rank-ladder, and seed
% setting from the expanded European-wind experiment.  It stops when the
% post-fit audit asks to enlarge no more than maxEnlargedFraction of the
% retained frequencies, subject to maxCycles adaptive fits.  Prediction uses
% the rank profile under which the final covariance parameters were fitted.
% The requested enlargements and an exact post-hoc audit are saved separately.

    if nargin<1, method="exact_kl"; end
    if nargin<2, foldID=1; end
    if nargin<3, maxCycles=3; end
    if nargin<4, maxEnlargedFraction=0.05; end
    method=string(method);
    if ~ismember(method,["exact_kl","solve_trace"])
        error('method must be exact_kl or solve_trace.');
    end
    if maxCycles<1 || maxCycles~=floor(maxCycles)
        error('maxCycles must be a positive integer.');
    end
    if maxEnlargedFraction<0 || maxEnlargedFraction>1
        error('maxEnlargedFraction must lie in [0,1].');
    end

    C=ghcnh_trace_config(); [cfg,subset,folds,A,E]=ghcnh_expanded_setup();
    foldRoot=fullfile(E.resultsDir,sprintf('fold_%02d',foldID));
    configurationTag=sprintf('maxcycles_%d_fraction_%04d',maxCycles, ...
        round(10000*maxEnlargedFraction));
    sensitivityRoot=fullfile(E.resultsDir,'relaxed_cycle_sensitivity', ...
        sprintf('fold_%02d',foldID),configurationTag);
    methodDir=fullfile(sensitivityRoot,char(method));
    if ~isfolder(methodDir), mkdir(methodDir); end
    logFile=fullfile(methodDir,'progress.log');
    totalTimer=tic;

    targetIndex=folds(foldID).targetIndex(:);
    observedIndex=setdiff((1:size(subset.U,2))',targetIndex,'stable');
    Yo=subset.U(:,observedIndex); Yt=subset.U(:,targetIndex);
    [observed,target,domainXy]=fold_geometry(subset,observedIndex, ...
        targetIndex,cfg,A);
    K=floor(size(Yo,1)/2)+1;
    allowedEnlargements=floor(maxEnlargedFraction*K+1e-12);
    ghcnh_log(logFile,char(method),sprintf([ ...
        'relaxed fold %d started: nTrain=%d nTest=%d T=%d; ' ...
        'maxCycles=%d allowedEnlargements=%d/%d'],foldID, ...
        numel(observedIndex),numel(targetIndex),size(Yo,1),maxCycles, ...
        allowedEnlargements,K));

    timer=tic;
    master=adaptive_basis_bank(domainXy,cfg.hsgp.boundaryFactor, ...
        C.referenceJ,A,cfg);
    basisSeconds=toc(timer);

    pilotFile=fullfile(foldRoot,'shared_pilot_fit.mat');
    if ~isfile(pilotFile)
        error('Shared pilot checkpoint is absent: %s',pilotFile);
    end
    S=load(pilotFile,'pilotFit','pilotSeconds');
    pilotFit=S.pilotFit; pilotSeconds=S.pilotSeconds;
    ghcnh_log(logFile,char(method),sprintf( ...
        'loaded shared pilot checkpoint (%.2f min)',pilotSeconds/60));

    timer=tic;
    [selectedJ,initialAudit,initialAux]=run_audit(method,pilotFit.theta, ...
        observed.xyCovariance,cfg,master,C,A,logFile, ...
        sprintf('relaxed_%s_fold%d_initial',method,foldID));
    initialAuditSeconds=toc(timer);
    write_audit(methodDir,'initial',selectedJ,initialAudit,initialAux);

    fit=pilotFit; rows=cell(maxCycles,1); fitTotal=0;
    auditTotal=initialAuditSeconds; relaxedStop=false; strictStable=false;
    finalRequestedJ=selectedJ; finalAudit=initialAudit;
    for cycle=1:maxCycles
        currentJ=selectedJ;
        ghcnh_log(logFile,char(method),sprintf( ...
            'cycle %d fit starting; J=%d--%d',cycle,min(currentJ),max(currentJ)));
        timer=tic;
        fit=adaptive_hsgp_fit(Yo,observed.xyCovariance,observed.design, ...
            cfg,master,currentJ,fit.theta,cfg.hsgp.optimizer, ...
            cfg.seed.optimizer+1000*foldID+100*cycle,A);
        fitSeconds=toc(timer); fitTotal=fitTotal+fitSeconds;

        timer=tic;
        [requiredJ,audit,aux]=run_audit(method,fit.theta, ...
            observed.xyCovariance,cfg,master,C,A,logFile, ...
            sprintf('relaxed_%s_fold%d_cycle%d',method,foldID,cycle));
        auditSeconds=toc(timer); auditTotal=auditTotal+auditSeconds;
        enlarged=max(currentJ,requiredJ);
        nEnlarged=sum(enlarged>currentJ);
        strictStable=nEnlarged==0;
        relaxedStop=nEnlarged<=allowedEnlargements;
        rows{cycle}=table(cycle,fitSeconds,auditSeconds,fit.objective, ...
            fit.totalFunctionEvaluations,fit.theta.noiseVariance,fit.exitFlag, ...
            fit.convergenceCertified,fit.boundaryHit,min(currentJ), ...
            median(currentJ),max(currentJ),nEnlarged,nEnlarged/K, ...
            strictStable,relaxedStop,'VariableNames',{'cycle','fitSeconds', ...
            'auditSeconds','objective','functionEvaluations','noiseVariance', ...
            'exitFlag','convergenceCertified','boundaryHit','minimumJ', ...
            'medianJ','maximumJ','nFrequenciesEnlarged', ...
            'fractionFrequenciesEnlarged','strictStable','relaxedStop'});
        write_audit(methodDir,sprintf('cycle_%d',cycle),requiredJ,audit,aux);
        ghcnh_log(logFile,char(method),sprintf([ ...
            'cycle %d fit %.2f min audit %.2f s enlarged=%d/%d ' ...
            'strict=%d relaxed=%d objective=%.7g'],cycle,fitSeconds/60, ...
            auditSeconds,nEnlarged,K,strictStable,relaxedStop,fit.objective));
        selectedJ=currentJ;
        finalRequestedJ=requiredJ; finalAudit=audit;
        if relaxedStop, break; end
        selectedJ=enlarged;
    end
    history=vertcat(rows{1:cycle});
    stoppedAtCap=~relaxedStop && cycle==maxCycles;
    fit.selectedJ=selectedJ;

    timer=tic;
    prediction=predict_adaptive_hsgp(Yo,observed.xyCovariance, ...
        target.xyCovariance,observed.design,target.design,fit,cfg,master,A);
    predictionSeconds=toc(timer);
    [metrics,siteMetrics,predictions]=score_ghcnh_predictions( ...
        Yt,prediction,subset,targetIndex,foldID,cfg, ...
        "relaxed_"+method);
    procedureSeconds=toc(totalTimer);
    equivalentFirstUseSeconds=procedureSeconds+pilotSeconds;

    % Validation-only exact audit.  It is excluded from the procedure time.
    if method=="exact_kl"
        exactRequiredJ=finalRequestedJ;
        exactAudit=finalAudit;
        validationAuditSeconds=0;
    else
        timer=tic;
        [exactRequiredJ,exactAudit,exactAux]=exact_rank_audit(fit.theta, ...
            observed.xyCovariance,cfg,master,C,A,logFile, ...
            sprintf('relaxed_%s_fold%d_posthoc_exact',method,foldID));
        validationAuditSeconds=toc(timer);
        writetable(exactAudit,fullfile(methodDir,'posthoc_exact_audit.csv'));
        writetable(exactAux,fullfile(methodDir,'posthoc_exact_audit_aux.csv'));
    end
    exactViolationCount=sum(exactRequiredJ>selectedJ);
    maxExactKlPerSite=maximum_selected_exact_kl(exactAudit,selectedJ,C.referenceJ);

    profile=table((0:K-1)',selectedJ,selectedJ.^2,finalRequestedJ, ...
        finalRequestedJ.^2,exactRequiredJ,exactRequiredJ.^2, ...
        'VariableNames',{'frequencyIndex','fittedJ','fittedM', ...
        'ruleRequestedJ','ruleRequestedM','exactRequiredJ','exactRequiredM'});
    timing=table(foldID,basisSeconds,pilotSeconds,initialAuditSeconds, ...
        fitTotal,auditTotal,predictionSeconds,procedureSeconds, ...
        equivalentFirstUseSeconds,validationAuditSeconds, ...
        'VariableNames',{'foldID','basisSeconds','pilotFitSeconds', ...
        'initialAuditSeconds','adaptiveFitSeconds','totalAuditSeconds', ...
        'predictionSeconds','procedureSeconds','equivalentFirstUseSeconds', ...
        'validationAuditSeconds'});
    stopping=table(foldID,string(method),maxCycles,maxEnlargedFraction, ...
        allowedEnlargements,cycle,strictStable,relaxedStop,stoppedAtCap, ...
        history.nFrequenciesEnlarged(end),exactViolationCount, ...
        maxExactKlPerSite,'VariableNames',{'foldID','method','maxCycles', ...
        'maxEnlargedFraction','allowedEnlargements','cyclesUsed', ...
        'strictStable','relaxedStop','stoppedAtCap','ruleViolationCount', ...
        'exactKlViolationCount','maximumExactKlPerSiteAtFittedRanks'});
    parameters=struct2table(fit.theta,'AsArray',true);
    parameters.foldID=foldID; parameters.objective=fit.objective;
    parameters.exitFlag=fit.exitFlag;
    parameters.convergenceCertified=fit.convergenceCertified;
    parameters.boundaryHit=fit.boundaryHit;

    writetable(history,fullfile(methodDir,'iterations.csv'));
    writetable(profile,fullfile(methodDir,'rank_profile.csv'));
    writetable(timing,fullfile(methodDir,'timing.csv'));
    writetable(stopping,fullfile(methodDir,'stopping_summary.csv'));
    writetable(parameters,fullfile(methodDir,'parameters.csv'));
    writetable(metrics,fullfile(methodDir,'metrics.csv'));
    writetable(siteMetrics,fullfile(methodDir,'site_metrics.csv'));
    writetable(predictions,fullfile(methodDir,'predictions.csv'));
    task=struct('method',method,'foldID',foldID,'metrics',metrics, ...
        'siteMetrics',siteMetrics,'predictions',predictions,'timing',timing, ...
        'stopping',stopping,'parameters',parameters,'rankProfile',profile, ...
        'history',history,'fit',fit,'prediction',prediction, ...
        'observedIndex',observedIndex,'targetIndex',targetIndex);
    checkpointFile=fullfile(methodDir,'complete_task.mat');
    try
        save(checkpointFile,'task','C','cfg','-v7.3');
    catch ME
        ghcnh_log(logFile,char(method),sprintf([ ...
            'nonessential complete-task checkpoint failed after CSV outputs ' ...
            'were written: %s'],ME.message));
        warning('ghcnh:checkpointSaveFailed', ...
            'CSV outputs are complete, but checkpoint save failed: %s', ...
            ME.message);
    end
    ghcnh_log(logFile,char(method),sprintf([ ...
        'relaxed fold %d complete: cycles=%d RMSE=%.6f first-use=%.2f min ' ...
        'ruleViolations=%d exactViolations=%d maxExactKL=%.6g'],foldID, ...
        cycle,metrics.rmse,equivalentFirstUseSeconds/60, ...
        history.nFrequenciesEnlarged(end),exactViolationCount,maxExactKlPerSite));
end

function [selectedJ,audit,aux]=run_audit(method,theta,xy,cfg,master,C,A,logFile,label)
    if method=="exact_kl"
        [selectedJ,audit,aux]=exact_rank_audit(theta,xy,cfg,master,C,A, ...
            logFile,label);
    else
        [selectedJ,audit,aux]=solve_trace_rank_audit(theta,xy,cfg,master,C,A, ...
            logFile,label);
    end
end

function [observed,target,domainXy]=fold_geometry(subset,observedIndex, ...
        targetIndex,cfg,A)
    allCenter=mean(subset.xyKm,1);
    domainXy=(subset.xyKm-allCenter)/cfg.spatial.coordinateScaleKm;
    observed.xyCovariance=domainXy(observedIndex,:);
    target.xyCovariance=domainXy(targetIndex,:);
    center=mean(subset.xyKm(observedIndex,:),1);
    scale=std(subset.xyKm(observedIndex,:),0,1); scale(scale==0)=1;
    observed.design=A.makeDesign(size(subset.U,1), ...
        (subset.xyKm(observedIndex,:)-center)./scale);
    target.design=A.makeDesign(size(subset.U,1), ...
        (subset.xyKm(targetIndex,:)-center)./scale);
end

function write_audit(folder,prefix,selectedJ,audit,aux)
    writetable(audit,fullfile(folder,[prefix '_audit.csv']));
    writetable(aux,fullfile(folder,[prefix '_audit_aux.csv']));
    profile=table((0:numel(selectedJ)-1)',selectedJ,selectedJ.^2, ...
        'VariableNames',{'frequencyIndex','selectedJ','selectedM'});
    writetable(profile,fullfile(folder,[prefix '_rank_profile.csv']));
end

function value=maximum_selected_exact_kl(audit,selectedJ,referenceJ)
    value=0;
    for k=1:numel(selectedJ)
        if selectedJ(k)>=referenceJ, local=0;
        else
            row=audit.frequencyIndex==k-1 & audit.candidateJ==selectedJ(k);
            if ~any(row)
                error('No exact-KL audit row for frequency %d and J=%d.', ...
                    k-1,selectedJ(k));
            end
            local=audit.exactKlPerSite(find(row,1));
        end
        value=max(value,local);
    end
end
