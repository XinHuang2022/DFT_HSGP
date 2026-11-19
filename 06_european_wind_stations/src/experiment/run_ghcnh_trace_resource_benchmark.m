function result = run_ghcnh_trace_resource_benchmark(foldID)
%RUN_GHCNH_TRACE_RESOURCE_BENCHMARK Time one objective and Algorithm-2 screen.

    if nargin<1, foldID=1; end
    C=ghcnh_trace_config();
    [cfg,subset,folds,A]=ghcnh_algorithm2_setup();
    fold=folds(foldID); target=fold.targetIndex(:);
    observed=setdiff((1:size(subset.U,2))',target,'stable');
    Yo=subset.U(:,observed); xyAllCentered=(subset.xyKm-mean(subset.xyKm,1));
    xy=xyAllCentered(observed,:)/cfg.spatial.coordinateScaleKm;
    meanCenter=mean(subset.xyKm(observed,:),1);
    meanScale=std(subset.xyKm(observed,:),0,1); meanScale(meanScale==0)=1;
    meanXy=(subset.xyKm(observed,:)-meanCenter)./meanScale;
    design=make_ghcnh_mean_design(size(Yo,1),meanXy);
    label=sprintf('fold%d_benchmark',foldID);
    foldDir=fullfile(C.resultsDir,sprintf('fold_%02d',foldID));
    if ~isfolder(foldDir), mkdir(foldDir); end
    logFile=fullfile(foldDir,'progress.log');

    timer=tic;
    master=adaptive_basis_bank(xyAllCentered/cfg.spatial.coordinateScaleKm, ...
        cfg.hsgp.boundaryFactor,C.referenceJ,A,cfg);
    basisSeconds=toc(timer);
    Phi=A.features(xy,master);
    basisBytes=whos('Phi');
    ghcnh_log(logFile,label,sprintf( ...
        'master basis J=%d M=%d built in %.2fs; training Phi %.1f MB', ...
        C.referenceJ,C.referenceJ^2,basisSeconds,basisBytes.bytes/2^20));

    pilotBasis=subset_adaptive_basis(master,C.initialJ);
    prepared=A.prepare(Yo,xy,design);
    prepared.hsgpPhi=A.features(xy,pilotBasis);
    timer=tic;
    [objective,details]=era5_profiled_whittle_theta_hsgp( ...
        cfg.start,prepared,cfg,pilotBasis);
    objectiveSeconds=toc(timer);
    ghcnh_log(logFile,label,sprintf( ...
        'one J=%d objective evaluation %.2fs; value %.6g', ...
        C.initialJ,objectiveSeconds,objective));

    timer=tic;
    [selectedJ,audit,spotChecks]=trace_rank_audit(cfg.start,xy,cfg, ...
        master,C,A,logFile,label);
    traceSeconds=toc(timer);
    ghcnh_log(logFile,label,sprintf( ...
        'trace screen plus %d exact spot checks %.2fs; J range %d--%d', ...
        height(spotChecks),traceSeconds,min(selectedJ),max(selectedJ)));
    writetable(audit,fullfile(foldDir,'start_parameter_trace_audit.csv'));
    writetable(spotChecks,fullfile(foldDir,'start_parameter_spot_checks.csv'));
    profile=table((0:numel(selectedJ)-1)',selectedJ,selectedJ.^2, ...
        'VariableNames',{'frequencyIndex','selectedJ','selectedM'});
    writetable(profile,fullfile(foldDir,'start_parameter_rank_profile.csv'));
    result=struct('foldID',foldID,'nTraining',numel(observed), ...
        'nTest',numel(target),'basisSeconds',basisSeconds, ...
        'basisMegabytes',basisBytes.bytes/2^20, ...
        'objectiveSeconds',objectiveSeconds,'objective',objective, ...
        'traceSeconds',traceSeconds,'selectedJ',selectedJ, ...
        'spotChecks',spotChecks,'details',details);
    save(fullfile(foldDir,'resource_benchmark.mat'),'result','C','cfg','-v7.3');
end
