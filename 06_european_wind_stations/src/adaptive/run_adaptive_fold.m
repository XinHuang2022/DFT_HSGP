function task = run_adaptive_fold(C, cfg, subset, fold, A)
%RUN_ADAPTIVE_FOLD Checkpointed Algorithm-1 fit and reconstruction.

    foldDir=fullfile(C.outputRoot,sprintf('fold_%02d',fold.foldID));
    if ~isfolder(foldDir), mkdir(foldDir); end
    finalFile=fullfile(foldDir,'complete_task.mat');
    if C.resume && isfile(finalFile)
        loaded=load(finalFile,'task'); task=loaded.task;
        fprintf('%s fold %d loaded from checkpoint.\n',C.dataset,fold.foldID);
        return;
    end
    logFile=fullfile(foldDir,'progress.log');
    totalTimer=tic;
    log_message(logFile,fold.foldID,'setup','started');
    targetIndex=fold.targetIndex(:);
    observedIndex=setdiff((1:size(subset.U,2))',targetIndex,'stable');
    Yo=subset.U(:,observedIndex); Yt=subset.U(:,targetIndex);
    [observed,target]=fold_geometry(subset,observedIndex,targetIndex,cfg,A);
    [varianceReference,~,~]=A.empiricalVariance(Yo,subset.xyKm(observedIndex,:));
    foldCfg=cfg;
    foldCfg.bounds.noiseVariance(1)=max(cfg.bounds.noiseVariance(1), ...
        cfg.hsgp.nuggetResidualVarianceFraction*varianceReference);
    Jmax=max(C.referenceJ);
    domainXy=(subset.xyKm-mean(subset.xyKm,1))./cfg.spatial.coordinateScaleKm;
    basisTimer=tic;
    master=adaptive_basis_bank(domainXy,cfg.hsgp.boundaryFactor,Jmax,A,foldCfg);
    basisSeconds=toc(basisTimer);
    validate_nested_bank(master,unique([C.initialJ,C.candidateJ,C.referenceJ]));
    log_message(logFile,fold.foldID,'basis',sprintf('Jmax=%d Mmax=%d %.1fs', ...
        Jmax,Jmax^2,basisSeconds));

    seed=cfg.seed.optimizer+1009*fold.foldID;
    pilotFile=fullfile(foldDir,'pilot_fit.mat');
    trustedPilotFile=fullfile(C.trustedPilotDir, ...
        sprintf(C.trustedPilotName,fold.foldID));
    pilotReused=false;
    if C.resume && isfile(pilotFile)
        S=load(pilotFile,'pilotFit','pilotSeconds');
        pilotFit=S.pilotFit; pilotSeconds=S.pilotSeconds;
        pilotReused=true;
        log_message(logFile,fold.foldID,'pilot','loaded checkpoint');
    elseif C.reuseTrustedPilot && isfile(trustedPilotFile)
        trusted=load(trustedPilotFile);
        if C.dataset=="era5"
            pilotFit=trusted.checkpoint.fullFit;
            pilotSeconds=trusted.checkpoint.fitSeconds;
        else
            pilotFit=trusted.task.fullFit;
            pilotSeconds=trusted.task.timing.separableFitSeconds + ...
                trusted.task.timing.fullFitSeconds;
        end
        pilotSource="trusted_checkpoint";
        pilotReused=true;
        save(pilotFile,'pilotFit','pilotSeconds','pilotSource','-v7.3');
        log_message(logFile,fold.foldID,'pilot',sprintf( ...
            'trusted J=%d checkpoint; recorded %.1f min', ...
            C.initialJ,pilotSeconds/60));
    else
        pilotTimer=tic;
        pilotBasis=subset_adaptive_basis(master,C.initialJ);
        separable=A.fitFixed(Yo,observed.xyCovariance,observed.design, ...
            foldCfg,pilotBasis,[],1,[],foldCfg.hsgp.optimizer,seed+1);
        initialTheta=separable.theta;
        initialTheta.epsilon=min(foldCfg.bounds.epsilon(2), ...
            max(foldCfg.bounds.epsilon(1),foldCfg.start.epsilon));
        pilotFit=A.fitFixed(Yo,observed.xyCovariance,observed.design, ...
            foldCfg,pilotBasis,initialTheta,[],[], ...
            foldCfg.hsgp.optimizer,seed+2);
        pilotSeconds=toc(pilotTimer);
        pilotSource="new_fit";
        save(pilotFile,'pilotFit','pilotSeconds','pilotSource','-v7.3');
        log_message(logFile,fold.foldID,'pilot',sprintf( ...
            'objective %.6g, %.1f min',pilotFit.objective,pilotSeconds/60));
    end

    auditTimer=tic;
    [currentJ,initialAudit,initialReference]=adaptive_rank_audit( ...
        pilotFit.theta,observed.xyCovariance,foldCfg,master,C,A, ...
        logFile,sprintf('fold%d_initial',fold.foldID));
    initialAuditSeconds=toc(auditTimer);
    initialAudit.foldID=repmat(fold.foldID,height(initialAudit),1);
    initialAudit.cycle=zeros(height(initialAudit),1);
    initialReference.foldID=repmat(fold.foldID,height(initialReference),1);
    initialReference.cycle=zeros(height(initialReference),1);
    writetable(initialAudit,fullfile(foldDir,'initial_rank_audit.csv'));
    writetable(initialReference,fullfile(foldDir,'initial_reference_stability.csv'));
    writetable(rank_profile_table(currentJ,foldCfg,fold.foldID,0), ...
        fullfile(foldDir,'initial_selected_rank_profile.csv'));

    if C.auditOnly || ~initialReference.referenceStable(end)
        if C.auditOnly && initialReference.referenceStable(end)
            stopStatus="audit_only_reference_stabilized";
            stopMessage='audit-only run completed with stabilized reference';
        elseif C.auditOnly
            stopStatus="audit_only_reference_not_stabilized";
            stopMessage='audit-only run completed without stabilized reference';
        else
            stopStatus="reference_not_stabilized";
            stopMessage='internal reference ladder did not stabilize';
        end
        totalSeconds=toc(totalTimer);
        equivalentFirstUseSeconds=totalSeconds+pilotReused*pilotSeconds;
        rankProfile=rank_profile_table(currentJ,foldCfg,fold.foldID,0);
        timing=table(fold.foldID,basisSeconds,pilotSeconds,initialAuditSeconds, ...
            0,0,0,0,totalSeconds,equivalentFirstUseSeconds,pilotReused,0,false, ...
            stopStatus,'VariableNames',{'foldID', ...
            'basisSeconds','pilotFitSeconds','initialAuditSeconds', ...
            'adaptiveFitSeconds','verificationSeconds','finalCorrectionSeconds', ...
            'predictionSeconds','actualRunSeconds','equivalentFirstUseSeconds', ...
            'pilotLoadedFromTrustedCheckpoint','nAdaptiveCycles', ...
            'profileStabilized','runStatus'});
        parameters=theta_table(pilotFit.theta,fold.foldID,pilotFit.objective, ...
            pilotFit.convergenceCertified,pilotFit.boundaryHit);
        task=struct('dataset',C.dataset,'foldID',fold.foldID,'metrics',table(), ...
            'predictions',table(),'rankProfile',rankProfile,'iterations',table(), ...
            'audits',initialAudit,'references',initialReference,'timing',timing, ...
            'parameters',parameters,'fit',pilotFit,'observedIndex',observedIndex, ...
            'targetIndex',targetIndex,'version',C.version, ...
            'status',stopStatus);
        writetable(rankProfile,fullfile(foldDir,'final_rank_profile.csv'));
        writetable(timing,fullfile(foldDir,'timing.csv'));
        writetable(parameters,fullfile(foldDir,'parameters.csv'));
        save(finalFile,'task','-v7.3');
        log_message(logFile,fold.foldID,'stopped',stopMessage);
        return;
    end

    adaptiveFit=[]; iterationRows=cell(C.maxAdaptiveFits,1);
    verificationAudits=cell(C.maxAdaptiveFits,1);
    verificationReferences=cell(C.maxAdaptiveFits,1);
    fitSecondsTotal=0; verificationSecondsTotal=0; stabilized=false;
    for cycle=1:C.maxAdaptiveFits
        cycleTimer=tic;
        if cycle==1
            optimizer=foldCfg.hsgp.optimizer;
            if all(currentJ==max(C.referenceJ))
                % If the audit retains the largest available reference at
                % every frequency, there is no adaptive dimension reduction.
                % Use a single warm-started bounded refit and record whether
                % it converges, rather than paying for redundant multistarts.
                optimizer.numStarts=1;
                optimizer.maxIterations=min(160,optimizer.maxIterations);
                optimizer.maxFunctionEvaluations=min(800,optimizer.maxFunctionEvaluations);
                optimizer.selectivePolishingEnabled=false;
                log_message(logFile,fold.foldID,'optimizer_fallback', ...
                    'all frequencies use maximum reference; one warm start');
            end
        else
            optimizer=foldCfg.hsgp.optimizer;
            optimizer.numStarts=1;
            optimizer.maxIterations=min(160,optimizer.maxIterations);
            optimizer.maxFunctionEvaluations=min(800,optimizer.maxFunctionEvaluations);
            optimizer.selectivePolishingEnabled=false;
        end
        if cycle==1, startTheta=pilotFit.theta; else, startTheta=adaptiveFit.theta; end
        log_message(logFile,fold.foldID,sprintf('adaptive_fit_%d_start',cycle), ...
            sprintf('starts=%d, iterations<=%d, evaluations<=%d, J range=%d--%d', ...
            optimizer.numStarts,optimizer.maxIterations, ...
            optimizer.maxFunctionEvaluations,min(currentJ),max(currentJ)));
        adaptiveFit=adaptive_hsgp_fit(Yo,observed.xyCovariance, ...
            observed.design,foldCfg,master,currentJ, ...
            startTheta,optimizer, ...
            seed+100*cycle,A);
        fitSeconds=toc(cycleTimer); fitSecondsTotal=fitSecondsTotal+fitSeconds;
        log_message(logFile,fold.foldID,sprintf('adaptive_fit_%d',cycle), ...
            sprintf('objective %.6g, %.1f min',adaptiveFit.objective,fitSeconds/60));
        verifyTimer=tic;
        [requiredJ,audit,reference]=adaptive_rank_audit( ...
            adaptiveFit.theta,observed.xyCovariance,foldCfg,master,C,A, ...
            logFile,sprintf('fold%d_cycle%d',fold.foldID,cycle));
        verifySeconds=toc(verifyTimer);
        verificationSecondsTotal=verificationSecondsTotal+verifySeconds;
        enlarged=max(currentJ,requiredJ);
        nEnlarged=sum(enlarged>currentJ); stabilized=nEnlarged==0;
        iterationRows{cycle}=table(fold.foldID,cycle,fitSeconds, ...
            verifySeconds,adaptiveFit.objective, ...
            adaptiveFit.totalFunctionEvaluations,min(currentJ), ...
            median(currentJ),max(currentJ),nEnlarged,stabilized, ...
            'VariableNames',{'foldID','cycle','fitSeconds', ...
            'verificationSeconds','objective','functionEvaluations', ...
            'minimumJ','medianJ','maximumJ','nFrequenciesEnlarged', ...
            'stabilized'});
        audit.foldID=repmat(fold.foldID,height(audit),1);
        audit.cycle=repmat(cycle,height(audit),1);
        reference.foldID=repmat(fold.foldID,height(reference),1);
        reference.cycle=repmat(cycle,height(reference),1);
        verificationAudits{cycle}=audit; verificationReferences{cycle}=reference;
        checkpoint=struct('cycle',cycle,'fit',adaptiveFit,'currentJ',currentJ, ...
            'requiredJ',requiredJ,'enlargedJ',enlarged,'stabilized',stabilized);
        save(fullfile(foldDir,sprintf('adaptive_cycle_%d.mat',cycle)), ...
            'checkpoint','-v7.3');
        if stabilized, break; end
        currentJ=enlarged;
    end
    nCycles=cycle;
    if ~stabilized
        % The final enlargement is fitted once so prediction never uses a
        % rank profile that was not used for parameter estimation.
        currentJ=enlarged;
        optimizer=foldCfg.hsgp.optimizer; optimizer.numStarts=1;
        optimizer.maxIterations=min(160,optimizer.maxIterations);
        optimizer.maxFunctionEvaluations=min(800,optimizer.maxFunctionEvaluations);
        optimizer.selectivePolishingEnabled=false;
        timer=tic;
        adaptiveFit=adaptive_hsgp_fit(Yo,observed.xyCovariance, ...
            observed.design,foldCfg,master,currentJ,adaptiveFit.theta, ...
            optimizer,seed+999,A);
        finalCorrectionSeconds=toc(timer); fitSecondsTotal=fitSecondsTotal+finalCorrectionSeconds;
        log_message(logFile,fold.foldID,'final_correction',sprintf( ...
            'profile not stabilized; enlarged profile fitted in %.1f min', ...
            finalCorrectionSeconds/60));
    else
        finalCorrectionSeconds=0;
    end

    predictionTimer=tic;
    prediction=predict_adaptive_hsgp(Yo,observed.xyCovariance, ...
        target.xyCovariance,observed.design,target.design,adaptiveFit, ...
        foldCfg,master,A);
    predictionSeconds=toc(predictionTimer);
    [metrics,predictions]=score_predictions(C,subset,Yt,prediction, ...
        adaptiveFit,foldCfg,fold.foldID,targetIndex,observedIndex);
    totalSeconds=toc(totalTimer);
    equivalentFirstUseSeconds=totalSeconds + pilotReused*pilotSeconds;
    rankProfile=rank_profile_table(currentJ,foldCfg,fold.foldID,nCycles);
    iterations=vertcat(iterationRows{1:nCycles});
    allAudits=[initialAudit;vertcat(verificationAudits{1:nCycles})];
    allReferences=[initialReference;vertcat(verificationReferences{1:nCycles})];
    timing=table(fold.foldID,basisSeconds,pilotSeconds,initialAuditSeconds, ...
        fitSecondsTotal,verificationSecondsTotal,finalCorrectionSeconds, ...
        predictionSeconds,totalSeconds,equivalentFirstUseSeconds, ...
        pilotReused,nCycles,stabilized,"complete", ...
        'VariableNames',{'foldID','basisSeconds','pilotFitSeconds', ...
        'initialAuditSeconds','adaptiveFitSeconds','verificationSeconds', ...
        'finalCorrectionSeconds','predictionSeconds','actualRunSeconds', ...
        'equivalentFirstUseSeconds','pilotLoadedFromTrustedCheckpoint', ...
        'nAdaptiveCycles','profileStabilized','runStatus'});
    parameters=theta_table(adaptiveFit.theta,fold.foldID,adaptiveFit.objective, ...
        adaptiveFit.convergenceCertified,adaptiveFit.boundaryHit);
    writetable(metrics,fullfile(foldDir,'prediction_metrics_by_site.csv'));
    writetable(predictions,fullfile(foldDir,'predictions.csv'));
    writetable(rankProfile,fullfile(foldDir,'final_rank_profile.csv'));
    writetable(iterations,fullfile(foldDir,'adaptive_iterations.csv'));
    writetable(allAudits,fullfile(foldDir,'all_rank_audits.csv'));
    writetable(allReferences,fullfile(foldDir,'all_reference_stability.csv'));
    writetable(timing,fullfile(foldDir,'timing.csv'));
    writetable(parameters,fullfile(foldDir,'parameters.csv'));
    task=struct('dataset',C.dataset,'foldID',fold.foldID,'metrics',metrics, ...
        'predictions',predictions,'rankProfile',rankProfile, ...
        'iterations',iterations,'audits',allAudits,'references',allReferences, ...
        'timing',timing,'parameters',parameters,'fit',adaptiveFit, ...
        'observedIndex',observedIndex,'targetIndex',targetIndex, ...
        'version',C.version,'status',"complete");
    save(finalFile,'task','-v7.3');
    log_message(logFile,fold.foldID,'complete',sprintf('%.1f min',totalSeconds/60));
end

function [observed,target]=fold_geometry(subset,observedIndex,targetIndex,cfg,A)
    xo=subset.xyKm(observedIndex,:); xt=subset.xyKm(targetIndex,:);
    center=mean(xo,1); scale=std(xo,0,1); scale(scale<=sqrt(eps))=1;
    observed.xyMean=(xo-center)./scale; target.xyMean=(xt-center)./scale;
    covarianceCenter=mean(subset.xyKm,1);
    observed.xyCovariance=(xo-covarianceCenter)./cfg.spatial.coordinateScaleKm;
    target.xyCovariance=(xt-covarianceCenter)./cfg.spatial.coordinateScaleKm;
    observed.design=A.makeDesign(size(subset.U,1),observed.xyMean);
    target.design=A.makeDesign(size(subset.U,1),target.xyMean);
end

function validate_nested_bank(master,Jvalues)
    for J=Jvalues(:)'
        M=J^2;
        if any(max(master.indices(1:M,:),[],2)>J) || ...
                size(unique(master.indices(1:M,:),'rows'),1)~=M
            error('Nested square-shell basis validation failed at J=%d.',J);
        end
    end
end

function profile=rank_profile_table(selectedJ,cfg,foldID,cycle)
    K=numel(selectedJ); lambda=eq9_fourier_frequencies(cfg.data.expectedTimePoints);
    profile=table(repmat(foldID,K,1),(0:K-1)',lambda(1:K).', ...
        selectedJ(:),selectedJ(:).^2,repmat(cycle,K,1), ...
        'VariableNames',{'foldID','frequencyIndex','lambda', ...
        'selectedJ','selectedM','cycle'});
end

function row=theta_table(theta,foldID,objective,converged,boundaryHit)
    row=table(foldID,objective,theta.signalVariance,theta.alpha,theta.beta, ...
        theta.epsilon,theta.noiseVariance,converged,boundaryHit, ...
        'VariableNames',{'foldID','objective','signalVariance','alpha', ...
        'beta','epsilon','noiseVariance','convergenceCertified','boundaryHit'});
end

function [metrics,predictions]=score_predictions(C,subset,Yt,prediction,fit,cfg, ...
        foldID,targetIndex,observedIndex)
    nTarget=numel(targetIndex); metricRows=cell(nTarget,1); predictionRows=cell(nTarget,1);
    for j=1:nTarget
        value=eq9_prediction_metrics(Yt(:,j),prediction.mean(:,j), ...
            prediction.observationVariance(:,j),cfg.cv.nominalCoverage);
        site=subset.siteTable(targetIndex(j),:);
        nearest=min(sqrt(sum((subset.xyKm(observedIndex,:)- ...
            subset.xyKm(targetIndex(j),:)).^2,2)));
        metricRows{j}=table("adaptive_hsgp",foldID,string(site.stationID), ...
            site.gridRow,site.gridColumn,nearest,value.rmse,value.mae, ...
            value.correlation,value.nlpd,value.crps,value.coverage, ...
            value.meanIntervalWidth,'VariableNames',{'method','foldID', ...
            'targetStation','gridRow','gridColumn','nearestTrainingDistanceKm', ...
            'rmse','mae','correlation','nlpd','crps','coverage95', ...
            'meanIntervalWidth95'});
        T=size(Yt,1);
        predictionRows{j}=table(repmat("adaptive_hsgp",T,1), ...
            repmat(foldID,T,1),repmat(string(site.stationID),T,1), ...
            subset.timeIndex,subset.relativeHours,Yt(:,j),prediction.mean(:,j), ...
            prediction.observationVariance(:,j), ...
            'VariableNames',{'method','foldID','targetStation','timeIndex', ...
            'relativeHours','truthUWindMps','predictionUWindMps', ...
            'observationVariance'});
    end
    metrics=vertcat(metricRows{:}); predictions=vertcat(predictionRows{:});
end

function log_message(file,fold,stage,message)
    stamp=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    line=sprintf('%s\tfold_%02d\t%s\t%s\n',stamp,fold,stage,message);
    fid=fopen(file,'a'); if fid>=0,fprintf(fid,'%s',line);fclose(fid);end
    fprintf('[%s] fold %02d %s: %s\n',stamp,fold,stage,message);
end
