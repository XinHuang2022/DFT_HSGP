function task=run_ghcnh_expanded_method(method,foldID,maxCycles)
%RUN_GHCNH_EXPANDED_METHOD Controlled adaptive-rank fold fit.

    method=string(method);
    if ~ismember(method,["exact_kl","trace_bound","solve_trace"])
        error('Method must be exact_kl, trace_bound, or solve_trace.');
    end
    if nargin<2, foldID=1; end
    C=ghcnh_trace_config(); [cfg,subset,folds,A,E]=ghcnh_expanded_setup();
    if nargin<3 || isempty(maxCycles), maxCycles=E.maxCycles; end
    foldRoot=fullfile(E.resultsDir,sprintf('fold_%02d',foldID));
    methodDir=fullfile(foldRoot,char(method));
    if ~isfolder(methodDir), mkdir(methodDir); end
    completeFile=fullfile(methodDir,'complete_task.mat');
    if isfile(completeFile)
        S=load(completeFile,'task'); task=S.task;
        if task.rankStable, return; end
    end
    logFile=fullfile(methodDir,'progress.log'); totalTimer=tic;
    targetIndex=folds(foldID).targetIndex(:);
    observedIndex=setdiff((1:size(subset.U,2))',targetIndex,'stable');
    Yo=subset.U(:,observedIndex); Yt=subset.U(:,targetIndex);
    [observed,target,domainXy]=fold_geometry(subset,observedIndex, ...
        targetIndex,cfg,A);
    ghcnh_log(logFile,char(method),sprintf( ...
        'fold %d started: nTrain=%d nTest=%d T=%d',foldID, ...
        numel(observedIndex),numel(targetIndex),size(Yo,1)));

    timer=tic;
    master=adaptive_basis_bank(domainXy,cfg.hsgp.boundaryFactor, ...
        C.referenceJ,A,cfg); basisSeconds=toc(timer);
    pilotFile=fullfile(foldRoot,'shared_pilot_fit.mat');
    if isfile(pilotFile)
        S=load(pilotFile,'pilotFit','pilotSeconds');
        pilotFit=S.pilotFit; pilotSeconds=S.pilotSeconds; pilotReused=true;
    else
        pilotBasis=subset_adaptive_basis(master,C.initialJ); timer=tic;
        pilotFit=A.fitFixed(Yo,observed.xyCovariance,observed.design, ...
            cfg,pilotBasis,cfg.start,[],[],cfg.hsgp.optimizer, ...
            cfg.seed.optimizer+1000*foldID);
        pilotSeconds=toc(timer); pilotReused=false;
        save(pilotFile,'pilotFit','pilotSeconds','-v7.3');
    end

    checkpointFile=fullfile(methodDir,'checkpoint.mat');
    rows=cell(maxCycles,1); stable=false; priorMethodSeconds=0;
    resume=false;
    if isfile(checkpointFile)
        S=load(checkpointFile,'checkpoint'); checkpoint=S.checkpoint;
        if ~checkpoint.rankStable && checkpoint.cycle<maxCycles
            resume=true; fit=checkpoint.fit;
            selectedJ=checkpoint.selectedJ;
            previousHistory=checkpoint.history;
            for q=1:height(previousHistory), rows{q}=previousHistory(q,:); end
            fitTotal=sum(previousHistory.fitSeconds);
            initialAuxFile=fullfile(methodDir,'initial_audit_aux.csv');
            if isfile(initialAuxFile)
                initialAux=readtable(initialAuxFile);
                initialAuditSeconds=initialAux.totalSeconds(1);
            else
                initialAuditSeconds=0;
            end
            auditTotal=initialAuditSeconds+sum(previousHistory.auditSeconds);
            priorMethodSeconds=fitTotal+auditTotal;
            startCycle=checkpoint.cycle+1;
            ghcnh_log(logFile,char(method),sprintf( ...
                'fold %d resuming from checkpoint at cycle %d; next cycle=%d', ...
                foldID,checkpoint.cycle,startCycle));
        end
    end
    if ~resume
        timer=tic;
        [selectedJ,audit,aux]=run_audit(method,pilotFit.theta, ...
            observed.xyCovariance,cfg,master,C,A,logFile, ...
            sprintf('%s_fold%d_initial',method,foldID));
        initialAuditSeconds=toc(timer);
        write_audit(methodDir,'initial',selectedJ,audit,aux);
        fit=pilotFit; fitTotal=0; auditTotal=initialAuditSeconds;
        startCycle=1;
    end

    for cycle=startCycle:maxCycles
        ghcnh_log(logFile,char(method),sprintf( ...
            'fold %d cycle %d fit starting; J=%d--%d', ...
            foldID,cycle,min(selectedJ),max(selectedJ)));
        timer=tic;
        fit=adaptive_hsgp_fit(Yo,observed.xyCovariance,observed.design, ...
            cfg,master,selectedJ,fit.theta,cfg.hsgp.optimizer, ...
            cfg.seed.optimizer+1000*foldID+100*cycle,A);
        fitSeconds=toc(timer); fitTotal=fitTotal+fitSeconds;
        timer=tic;
        [requiredJ,audit,aux]=run_audit(method,fit.theta, ...
            observed.xyCovariance,cfg,master,C,A,logFile, ...
            sprintf('%s_fold%d_cycle%d',method,foldID,cycle));
        auditSeconds=toc(timer); auditTotal=auditTotal+auditSeconds;
        enlarged=max(selectedJ,requiredJ); nEnlarged=sum(enlarged>selectedJ);
        stable=nEnlarged==0;
        rows{cycle}=table(cycle,fitSeconds,auditSeconds,fit.objective, ...
            fit.totalFunctionEvaluations,fit.theta.noiseVariance,fit.exitFlag, ...
            fit.convergenceCertified,fit.boundaryHit,min(selectedJ), ...
            median(selectedJ),max(selectedJ),nEnlarged,stable, ...
            'VariableNames',{'cycle','fitSeconds','auditSeconds','objective', ...
            'functionEvaluations','noiseVariance','exitFlag', ...
            'convergenceCertified','boundaryHit','minimumJ','medianJ', ...
            'maximumJ','nFrequenciesEnlarged','rankStable'});
        write_audit(methodDir,sprintf('cycle_%d',cycle),requiredJ, ...
            audit,aux);
        selectedJ=enlarged;
        checkpoint=struct('fit',fit,'selectedJ',selectedJ,'cycle',cycle, ...
            'history',vertcat(rows{1:cycle}),'rankStable',stable);
        save(checkpointFile,'checkpoint','-v7.3');
        ghcnh_log(logFile,char(method),sprintf( ...
            'fold %d cycle %d fit %.2f min audit %.2f s enlarged=%d objective=%.7g', ...
            foldID,cycle,fitSeconds/60,auditSeconds,nEnlarged,fit.objective));
        if stable, break; end
    end

    history=vertcat(rows{1:cycle});
    profile=table((0:numel(selectedJ)-1)',selectedJ,selectedJ.^2, ...
        'VariableNames',{'frequencyIndex','selectedJ','selectedM'});
    if stable
        fit.selectedJ=selectedJ; timer=tic;
        prediction=predict_adaptive_hsgp(Yo,observed.xyCovariance, ...
            target.xyCovariance,observed.design,target.design,fit,cfg,master,A);
        predictionSeconds=toc(timer);
        [metrics,siteMetrics,predictions]=score_ghcnh_predictions( ...
            Yt,prediction,subset,targetIndex,foldID,cfg,method);
    else
        prediction=[]; predictionSeconds=0; metrics=table();
        siteMetrics=table(); predictions=table();
    end
    actualSeconds=toc(totalTimer)+priorMethodSeconds;
    firstUseSeconds=actualSeconds+double(pilotReused)*pilotSeconds;
    timing=table(foldID,basisSeconds,pilotSeconds,initialAuditSeconds, ...
        fitTotal,auditTotal,predictionSeconds,actualSeconds,firstUseSeconds, ...
        pilotReused,'VariableNames',{'foldID','basisSeconds', ...
        'pilotFitSeconds','initialAuditSeconds','adaptiveFitSeconds', ...
        'totalAuditSeconds','predictionSeconds','actualSeconds', ...
        'equivalentFirstUseSeconds','pilotReused'});
    parameters=struct2table(fit.theta,'AsArray',true);
    parameters.foldID=foldID; parameters.objective=fit.objective;
    parameters.exitFlag=fit.exitFlag;
    parameters.convergenceCertified=fit.convergenceCertified;
    parameters.boundaryHit=fit.boundaryHit;
    writetable(history,fullfile(methodDir,'iterations.csv'));
    writetable(profile,fullfile(methodDir,'rank_profile.csv'));
    writetable(timing,fullfile(methodDir,'timing.csv'));
    writetable(parameters,fullfile(methodDir,'parameters.csv'));
    if stable
        writetable(metrics,fullfile(methodDir,'metrics.csv'));
        writetable(siteMetrics,fullfile(methodDir,'site_metrics.csv'));
        writetable(predictions,fullfile(methodDir,'predictions.csv'));
    end
    task=struct('method',method,'foldID',foldID,'metrics',metrics, ...
        'siteMetrics',siteMetrics,'predictions',predictions,'timing',timing, ...
        'parameters',parameters,'rankProfile',profile,'history',history, ...
        'fit',fit,'prediction',prediction,'rankStable',stable, ...
        'observedIndex',observedIndex,'targetIndex',targetIndex, ...
        'panelSize',size(subset.U,2));
    save(completeFile,'task','C','cfg','-v7.3');
    ghcnh_log(logFile,char(method),sprintf( ...
        'fold %d finished stable=%d first-use=%.2f min', ...
        foldID,stable,firstUseSeconds/60));
end

function [selectedJ,audit,aux]=run_audit(method,theta,xy,cfg,master,C,A,logFile,label)
    if method=="exact_kl"
        [selectedJ,audit,aux]=exact_rank_audit(theta,xy,cfg,master,C,A, ...
            logFile,label);
    elseif method=="trace_bound"
        [selectedJ,audit,aux]=trace_rank_audit(theta,xy,cfg,master,C,A, ...
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
