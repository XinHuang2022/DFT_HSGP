function task=cerra_dense_stage1_task(cfg,referenceFile,M)
%CERRA_DENSE_STAGE1_TASK Fixed-theta full-data HSGP audit for one M.

    taskDir=fullfile(cfg.paths.stage1Dir,'tasks');ensure_directory(taskDir);
    taskFile=fullfile(taskDir,sprintf('M_%04d.mat',M));
    if isfile(taskFile)
        S=load(taskFile,'task');
        if isfield(S,'task')&&S.task.version==string(cfg.version)&&S.task.M==M
            fprintf('Stage 1 M=%d loaded from checkpoint.\n',M);task=S.task;return;
        end
    end
    timerValue=tic;S=load(referenceFile,'reference');R=S.reference;
    basis=make_cerra_hsgp_basis(R.domainXy,cfg.hsgp.boundaryFactor,M,cfg);
    auditPhi=cerra_hsgp_features(R.exactObserved.xyCovariance,basis);
    weights=cerra_hsgp_weights(R.exactFit.theta,R.auditLambda,basis,R.fullCfg);
    nFrequency=numel(R.auditLambda);latentError=zeros(nFrequency,1);
    diagonalRatio=zeros(nFrequency,1);
    for q=1:nFrequency
        approximate=cerra_hsgp_covariance(auditPhi,auditPhi,weights(:,q),0,true);
        exact=R.exactSignal(:,:,q);
        latentError(q)=norm(approximate-exact,'fro')/max(norm(exact,'fro'),eps);
        diagonalRatio(q)=median(diag(approximate)./max(diag(exact),eps));
    end
    fitTimer=tic;
    fixedFit=cerra_fixed_theta_fit("hsgp",R.exactFit.theta,R.Yfull, ...
        R.fullObserved,R.fullCfg,basis);
    likelihoodSeconds=toc(fitTimer);
    predictionTimer=tic;
    prediction=cerra_predict_fit("hsgp",R.Yfull,R.fullObserved, ...
        R.fullTarget,fixedFit,R.fullCfg,basis);
    predictionSeconds=toc(predictionTimer);
    metric=eq9_prediction_metrics(R.Ytarget(:),prediction.mean(:), ...
        prediction.observationVariance(:),cfg.cv.nominalCoverage);
    sorted=sort(latentError);p90=sorted(max(1,ceil(0.9*numel(sorted))));
    passes=metric.rmse<=R.exactMetric.rmse*(1+cfg.hsgp.screenRelativeRmseTolerance);
    summary=table(M,cfg.hsgp.boundaryFactor,basis.maximumIndexByDimension(1), ...
        basis.maximumIndexByDimension(2),basis.maximumRadialFrequency, ...
        median(latentError),p90,max(latentError),median(diagonalRatio), ...
        fixedFit.objective/numel(R.Yfull),metric.rmse,metric.mae, ...
        metric.correlation,metric.nlpd,metric.coverage,passes, ...
        likelihoodSeconds,predictionSeconds,toc(timerValue),cerra_worker_id(), ...
        'VariableNames',{'M','boundaryFactor','maximumIndexX','maximumIndexY', ...
        'maximumRadialFrequency','medianLatentCovarianceError', ...
        'p90LatentCovarianceError','maximumLatentCovarianceError', ...
        'medianDiagonalRatio','objectivePerObservation','pooledRmse','pooledMae', ...
        'correlation','nlpd','coverage95','passesPredictiveScreen', ...
        'likelihoodSeconds','predictionSeconds','elapsedSeconds','workerID'});
    [~,bySite]=cerra_dense_prediction_summary("hsgp_fixed_theta",M, ...
        numel(R.fullIndex),R.Ytarget,prediction,R.targetSiteTable, ...
        fixed_fit_for_summary(fixedFit,likelihoodSeconds),R.fullCfg);
    task=struct('version',string(cfg.version),'M',M,'summary',summary, ...
        'bySite',bySite,'basis',basis,'fixedFit',fixedFit,'prediction',prediction);
    save(taskFile,'task','-v7.3');
    fprintf('Stage 1 M=%d complete: RMSE %.4f, p90 covariance error %.4f.\n', ...
        M,metric.rmse,p90);
end

function fit=fixed_fit_for_summary(fit,seconds)
    fit.elapsedSeconds=seconds;fit.convergenceCertified=true;fit.boundaryHit=false;
end
