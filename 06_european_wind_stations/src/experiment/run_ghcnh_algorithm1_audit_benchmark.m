function result=run_ghcnh_algorithm1_audit_benchmark()
%RUN_GHCNH_ALGORITHM1_AUDIT_BENCHMARK Benchmark exact-KL rank screening.

    C=ghcnh_trace_config();
    [cfg,subset,folds,A]=ghcnh_algorithm2_setup();
    fold=folds(1); target=fold.targetIndex(:);
    observed=setdiff((1:size(subset.U,2))',target,'stable');
    center=mean(subset.xyKm,1);
    domain=(subset.xyKm-center)/cfg.spatial.coordinateScaleKm;
    xy=domain(observed,:);
    master=adaptive_basis_bank(domain,cfg.hsgp.boundaryFactor, ...
        C.referenceJ,A,cfg);
    pilot=load(fullfile(C.resultsDir,'fold_01','algorithm2_pilot_fit.mat'), ...
        'pilotFit');
    logFile=fullfile(C.resultsDir,'fold_01','algorithm1_progress.log');
    K=floor(cfg.data.expectedTimePoints/2)+1;
    indices=unique(round(linspace(1,K,12)));
    [selectedJ,audit,summary]=exact_rank_audit(pilot.pilotFit.theta, ...
        xy,cfg,master,C,A,logFile,'algorithm1_benchmark',indices);
    writetable(audit,fullfile(C.resultsDir,'fold_01', ...
        'algorithm1_benchmark_audit.csv'));
    writetable(summary,fullfile(C.resultsDir,'fold_01', ...
        'algorithm1_benchmark_timing.csv'));
    result=struct('frequencyIndices',indices,'selectedJ',selectedJ, ...
        'audit',audit,'summary',summary);
    save(fullfile(C.resultsDir,'fold_01','algorithm1_benchmark.mat'), ...
        'result','-v7.3');
end
