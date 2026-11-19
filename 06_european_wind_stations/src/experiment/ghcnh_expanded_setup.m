function [cfg,subset,folds,A,E]=ghcnh_expanded_setup()
%GHCNH_EXPANDED_SETUP Controlled setup for the expanded five-fold comparison.

    C=ghcnh_trace_config(); addpath(C.era5SourceRoot,C.adaptiveSourceRoot);
    cfg=era5_northsea_config(); cfg.suiteRoot=fileparts(mfilename('fullpath'));
    cfg.data.expectedTimePoints=C.expectedTimePoints;
    cfg.data.samplingHours=6; cfg.data.responseName='wind_speed';
    cfg.data.responseUnits='m/s'; cfg.spatial.coordinateScaleKm=1000;
    cfg.hsgp.boundaryFactor=2.5; cfg.hsgp.candidateJ=C.candidateJ;
    cfg.hsgp.factorization='auto';
    cfg.start.signalVariance=4; cfg.start.alpha=0.5; cfg.start.beta=0.25;
    cfg.start.epsilon=0.6; cfg.start.noiseVariance=0.25;
    cfg.bounds.signalVariance=[1e-4,200];
    cfg.bounds.noiseVariance=[1e-4,25];
    cfg.hsgp.optimizer.numStarts=1; cfg.hsgp.optimizer.display='off';
    cfg.hsgp.optimizer.maxIterations=160;
    cfg.hsgp.optimizer.maxFunctionEvaluations=600;
    cfg.hsgp.optimizer.tolX=1e-3; cfg.hsgp.optimizer.tolFun=1e-2;
    cfg.hsgp.optimizer.selectivePolishingEnabled=false;

    panelFile=fullfile(C.dataProcessedDir,'ghcnh_complete_panel_expanded.mat');
    if ~isfile(panelFile), subset=prepare_ghcnh_expanded_panel();
    else, S=load(panelFile,'subset'); subset=S.subset; end
    cfg.data.expectedStations=size(subset.U,2);
    folds=make_ghcnh_cv_folds(subset,C.nFolds);
    A.prefix="era5"; A.makeDesign=@make_ghcnh_mean_design;
    A.prepare=@prepare_era5_whittle_data; A.features=@era5_hsgp_features;
    A.weights=@era5_hsgp_weights; A.factor=@era5_hsgp_factor;
    A.solve=@era5_hsgp_solve; A.covariance=@era5_hsgp_covariance;
    A.pack=@era5_pack_constrained_parameters;
    A.unpack=@era5_unpack_constrained_parameters;
    A.evaluateMean=@evaluate_era5_mean; A.fitFixed=@fit_era5_eq9_hsgp;
    E.resultsDir=fullfile(C.resultsDir,sprintf( ...
        'expanded_complete_%04d',size(subset.U,2)));
    if ~isfolder(E.resultsDir), mkdir(E.resultsDir); end
    E.maxCycles=5; E.bootstrapReplicates=5000; E.spatialBlockKm=400;
end
