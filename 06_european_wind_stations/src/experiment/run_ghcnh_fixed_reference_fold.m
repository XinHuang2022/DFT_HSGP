function task=run_ghcnh_fixed_reference_fold(foldID)
%RUN_GHCNH_FIXED_REFERENCE_FOLD Fit J_ref at every temporal frequency.
% The comparator reuses the protected GHCNh panel, deterministic folds,
% mean specification, pilot initialization, optimizer settings, and scoring
% conventions of the adaptive-rank experiment. Its reported equivalent
% first-use time includes the shared pilot fit even when loaded from disk.

    if nargin<1, foldID=1; end
    C=ghcnh_trace_config(); [cfg,subset,folds,A,E]=ghcnh_expanded_setup();
    method="fixed_reference";
    foldRoot=fullfile(E.resultsDir,sprintf('fold_%02d',foldID));
    methodDir=fullfile(foldRoot,char(method));
    if ~isfolder(methodDir), mkdir(methodDir); end
    completeFile=fullfile(methodDir,'complete_task.mat');
    if isfile(completeFile)
        S=load(completeFile,'task'); task=S.task; return
    end

    logFile=fullfile(methodDir,'progress.log'); totalTimer=tic;
    targetIndex=folds(foldID).targetIndex(:);
    observedIndex=setdiff((1:size(subset.U,2))',targetIndex,'stable');
    Yo=subset.U(:,observedIndex); Yt=subset.U(:,targetIndex);
    [observed,target,domainXy]=fold_geometry(subset,observedIndex, ...
        targetIndex,cfg,A);
    ghcnh_log(logFile,char(method),sprintf( ...
        'fold %d started: nTrain=%d nTest=%d T=%d J=%d M=%d',foldID, ...
        numel(observedIndex),numel(targetIndex),size(Yo,1), ...
        C.referenceJ,C.referenceJ^2));

    timer=tic;
    master=adaptive_basis_bank(domainXy,cfg.hsgp.boundaryFactor, ...
        C.referenceJ,A,cfg);
    basis=subset_adaptive_basis(master,C.referenceJ);
    basisSeconds=toc(timer);

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

    ghcnh_log(logFile,char(method),sprintf( ...
        'fold %d fixed-reference fit starting from the shared pilot',foldID));
    timer=tic;
    fit=A.fitFixed(Yo,observed.xyCovariance,observed.design,cfg,basis, ...
        pilotFit.theta,[],[],cfg.hsgp.optimizer, ...
        cfg.seed.optimizer+1000*foldID+500);
    fitSeconds=toc(timer);
    ghcnh_log(logFile,char(method),sprintf( ...
        'fold %d fit completed in %.2f min; objective=%.7g exit=%d', ...
        foldID,fitSeconds/60,fit.objective,fit.exitFlag));

    timer=tic;
    prediction=predict_era5_eq9_hsgp(Yo,observed.xyCovariance, ...
        target.xyCovariance,observed.design,target.design,fit,cfg,basis);
    predictionSeconds=toc(timer);
    [metrics,siteMetrics,predictions]=score_ghcnh_predictions( ...
        Yt,prediction,subset,targetIndex,foldID,cfg,method);

    actualSeconds=toc(totalTimer);
    equivalentFirstUseSeconds=actualSeconds+double(pilotReused)*pilotSeconds;
    timing=table(foldID,basisSeconds,pilotSeconds,fitSeconds, ...
        predictionSeconds,actualSeconds,equivalentFirstUseSeconds,pilotReused, ...
        'VariableNames',{'foldID','basisSeconds','pilotFitSeconds', ...
        'fixedFitSeconds','predictionSeconds','actualSeconds', ...
        'equivalentFirstUseSeconds','pilotReused'});
    parameters=struct2table(fit.theta,'AsArray',true);
    parameters.foldID=foldID; parameters.objective=fit.objective;
    parameters.exitFlag=fit.exitFlag;
    parameters.convergenceCertified=fit.convergenceCertified;
    parameters.boundaryHit=fit.boundaryHit;
    writetable(metrics,fullfile(methodDir,'metrics.csv'));
    writetable(siteMetrics,fullfile(methodDir,'site_metrics.csv'));
    writetable(predictions,fullfile(methodDir,'predictions.csv'));
    writetable(timing,fullfile(methodDir,'timing.csv'));
    writetable(parameters,fullfile(methodDir,'parameters.csv'));
    task=struct('method',method,'foldID',foldID,'metrics',metrics, ...
        'siteMetrics',siteMetrics,'predictions',predictions,'timing',timing, ...
        'parameters',parameters,'fit',fit,'prediction',prediction, ...
        'observedIndex',observedIndex,'targetIndex',targetIndex, ...
        'panelSize',size(subset.U,2),'fixedJ',C.referenceJ, ...
        'fixedM',C.referenceJ^2);
    save(completeFile,'task','C','cfg','-v7.3');
    ghcnh_log(logFile,char(method),sprintf( ...
        'fold %d finished: RMSE=%.4f MAE=%.4f first-use=%.2f min', ...
        foldID,metrics.rmse,metrics.mae,equivalentFirstUseSeconds/60));
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
