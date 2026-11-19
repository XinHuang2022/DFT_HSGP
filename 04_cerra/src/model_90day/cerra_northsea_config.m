function cfg = cerra_northsea_config()
%CERRA_NORTHSEA_CONFIG Dense-grid HSGP versus coarse exact experiment.

    root = fileparts(mfilename('fullpath'));
    experimentRoot = fileparts(fileparts(root));
    cfg.version = "0.2.1";
    cfg.suiteRoot = root;
    dataOverride = string(getenv('CERRA_DATA_DIR'));
    if strlength(dataOverride) > 0
        cfg.paths.dataDir = char(dataOverride);
    else
        cfg.paths.dataDir = fullfile(experimentRoot, 'input', 'data');
    end
    resultsOverride = string(getenv('CERRA_REPRO_OUTPUT_DIR'));
    if strlength(resultsOverride) > 0
        cfg.paths.resultsDir = char(resultsOverride);
    else
        cfg.paths.resultsDir = fullfile(experimentRoot, 'results', 'current');
    end
    cfg.paths.stage1Dir = fullfile(cfg.paths.resultsDir, ...
        'stage1_fixed_theta_dense_screen');
    cfg.paths.stage2Dir = fullfile(cfg.paths.resultsDir, ...
        'stage2_selective_hsgp_refit');
    cfg.paths.stage3Dir = fullfile(cfg.paths.resultsDir, ...
        'stage3_coarse_exact_baselines');
    cfg.paths.smokeDir = fullfile(cfg.paths.resultsDir, 'laptop_smoke');

    cfg.data.windFileName = 'u10_full_native_time_by_site.csv';
    cfg.data.siteFileName = 'site_manifest_full_native.csv';
    cfg.data.timeFileName = 'time_manifest_6hourly.csv';
    cfg.data.expectedSourceTimePoints = 360;
    cfg.data.expectedFullStations = 3701;
    cfg.data.periodStartUtc = "2020-12-01T00:00:00Z";
    cfg.data.periodEndUtc = "2021-02-28T18:00:00Z";
    cfg.data.expectedTimePoints = 360;
    cfg.data.samplingHours = 6;
    cfg.data.responseName = 'east_west_u_wind';
    cfg.data.responseUnits = 'm/s';

    cfg.split.nTestSites = 300;
    cfg.split.exactSizes = [250 500 750 1000];
    cfg.split.design = 'deterministic_farthest_point';
    cfg.seed.split = 92017;
    cfg.seed.optimizer = 91401;

    cfg.spatial.coordinateScaleKm = 100;
    cfg.model.dimension = 2;
    cfg.model.nu = 2.5;
    cfg.model.deltaTime = 1;
    cfg.model.aliasTerms = 6;
    cfg.model.normalizationGridSize = 2048;

    cfg.bounds.signalVariance = [1e-3, 200];
    cfg.bounds.alpha = [0.003, 80];
    cfg.bounds.beta = [0.003, 3];
    cfg.bounds.epsilon = [0.02, 1 - 1e-6];
    cfg.bounds.logSpatialScaleContrast = ...
        [-0.5 * log(cfg.bounds.epsilon(2)), ...
         -0.5 * log(cfg.bounds.epsilon(1))];
    cfg.bounds.noiseVariance = [1e-5, 100];
    cfg.parameterization = 'log_spatial_scale_contrast';
    cfg.start = struct('signalVariance',16,'alpha',0.2,'beta',0.25, ...
        'epsilon',0.55,'noiseVariance',0.2,'nu',cfg.model.nu);
    cfg.likelihood.type = 'whittle_reml';
    cfg.mean.coefficientNames = ["intercept";"east";"north"];

    cfg.hsgp.boundaryFactor = 2;
    cfg.hsgp.candidateM = [250 500 750 1000];
    cfg.hsgp.truncation = 'radial_laplacian_eigenvalue';
    cfg.hsgp.factorization = 'feature_qr';
    cfg.hsgp.minimumReciprocalCondition = 1e-12;
    cfg.hsgp.maximumAllowedJitter = 1e-6;
    cfg.hsgp.minimumPredictiveVariance = 1e-10;
    cfg.hsgp.nuggetResidualVarianceFraction = 1e-6;
    cfg.hsgp.auditFrequencyCount = 24;
    cfg.hsgp.auditSpatialSites = 250;
    cfg.hsgp.screenRelativeRmseTolerance = 0.15;

    cfg.stage1.referenceExactSize = 250;
    cfg.stage1.optimizer = optimizer_settings(1,220,1100,true);
    cfg.stage2.nIndependentStarts = 3;
    cfg.stage2.nMultistartCandidates = 2;
    cfg.stage2.singleStartOptimizer = optimizer_settings(1,300,1800,false);
    cfg.stage3.exactSizes = [250 500];
    cfg.stage3.extendedExactSizes = [250 500 750 1000];
    cfg.stage3.nIndependentStarts = 3;
    cfg.stage3.singleStartOptimizer = optimizer_settings(1,300,1800,false);

    cfg.cv.nominalCoverage = 0.95;
    cfg.cv.minimumPredictiveVariance = 1e-10;
    cfg.cv.maximumAllowedImaginary = 1e-8;
    cfg.figures.visible = 'off';
    cfg.parallel.useParallel = true;
    cfg.parallel.autoStartPool = true;
    cfg.parallel.numWorkers = available_workers(6);
end

function n = available_workers(fallback)
    values = [str2double(getenv('CERRA_NUM_WORKERS')), ...
        str2double(getenv('SLURM_CPUS_PER_TASK'))];
    values = values(isfinite(values) & values >= 1);
    if isempty(values), n = fallback; else, n = floor(values(1)); end
end

function value = optimizer_settings(nStarts,maxIterations,maxEvaluations,polish)
    value.display = 'off';
    value.numStarts = nStarts;
    value.maxIterations = maxIterations;
    value.maxFunctionEvaluations = maxEvaluations;
    value.tolX = 3e-4;
    value.tolFun = 3e-4;
    value.randomStartScale = 1.25;
    value.nearBestObjectiveTolerance = 1e-3;
    value.selectivePolishingEnabled = polish;
    value.polishMaxIterations = 800;
    value.polishMaxFunctionEvaluations = 4000;
    value.polishTolX = 1e-4;
    value.polishTolFun = 1e-4;
    value.polishAcceptObjectiveTolerance = 1e-4;
end
