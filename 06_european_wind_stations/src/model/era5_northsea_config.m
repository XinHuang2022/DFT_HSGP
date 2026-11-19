function cfg = era5_northsea_config()
%ERA5_NORTHSEA_CONFIG Reproducible exact/HSGP North Sea experiment.

    root = fileparts(mfilename('fullpath'));
    cfg.version = "0.1.0";
    cfg.suiteRoot = root;
    cfg.paths.dataDir = fullfile(fileparts(root), 'Data');
    cfg.paths.resultsDir = fullfile(root, 'Results');
    cfg.paths.stage1Dir = fullfile(cfg.paths.resultsDir, ...
        'stage1_fixed_theta_J_audit');
    cfg.paths.stage2Dir = fullfile(cfg.paths.resultsDir, ...
        'stage2_hsgp_refit_screen');
    cfg.paths.stage3Dir = fullfile(cfg.paths.resultsDir, ...
        'stage3_paired_cv');

    cfg.data.windFileName = 'u10_time_by_site.csv';
    cfg.data.siteFileName = 'site_manifest_17x17.csv';
    cfg.data.timeFileName = 'time_manifest_6hourly.csv';
    cfg.data.expectedTimePoints = 360;
    cfg.data.expectedStations = 289;
    cfg.data.expectedGridSide = 17;
    cfg.data.samplingHours = 6;
    cfg.data.responseName = 'east_west_u_wind';
    cfg.data.responseUnits = 'm/s';

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
    cfg.start = struct('signalVariance', 16, 'alpha', 0.2, ...
        'beta', 0.25, 'epsilon', 0.55, 'noiseVariance', 0.2, ...
        'nu', cfg.model.nu);
    cfg.likelihood.type = 'whittle_reml';

    cfg.mean.coefficientNames = ["intercept"; "east"; "north"];
    cfg.cv.nFolds = 5;
    cfg.cv.holdoutPattern = 'balanced5';
    cfg.cv.nominalCoverage = 0.95;
    cfg.cv.minimumPredictiveVariance = 1e-10;
    cfg.cv.maximumAllowedImaginary = 1e-8;
    cfg.cv.timeBlocks = 4;

    cfg.seed.optimizer = 90401;
    cfg.optimizer = optimizer_settings(3, 300, 1800);
    cfg.stage1.optimizer = optimizer_settings(3, 300, 1800);
    cfg.stage2.optimizer = optimizer_settings(2, 240, 1300);
    cfg.stage3.optimizer = optimizer_settings(3, 300, 1800);

    cfg.parallel.useParallel = true;
    cfg.parallel.autoStartPool = true;
    cfg.parallel.numWorkers = available_workers(10);

    cfg.hsgp.boundaryFactor = 4;
    cfg.hsgp.candidateJ = [6 8 10 12 14 16 20 24 32 48];
    cfg.hsgp.factorization = 'auto';
    cfg.hsgp.minimumReciprocalCondition = 1e-12;
    cfg.hsgp.maximumAllowedJitter = 1e-6;
    cfg.hsgp.minimumPredictiveVariance = cfg.cv.minimumPredictiveVariance;
    cfg.hsgp.nuggetResidualVarianceFraction = 1e-6;
    cfg.hsgp.reportNuggetBoundaryTolerance = 0.02;
    cfg.hsgp.optimizer = cfg.stage2.optimizer;
    cfg.hsgp.auditFrequencyCount = 24;
    cfg.hsgp.acceptance.maximumP90LatentCovarianceError = 0.05;
    cfg.hsgp.acceptance.maximumLatentCovarianceError = 0.10;
    cfg.hsgp.acceptance.maximumAbsoluteObjectiveDifferencePerObservation = 0.05;
    cfg.hsgp.acceptance.minimumMedianDiagonalRatio = 0.95;
    cfg.hsgp.acceptance.maximumMedianDiagonalRatio = 1.05;

    cfg.screening.foldID = 1;
    cfg.stage2.candidateJ = [8 12 16 20 24 32 48];
    cfg.final.selectedJ = [];
    cfg.final.includeAdjacentJ = true;
    cfg.figures.visible = 'off';
end

function n = available_workers(fallback)
    values = [str2double(getenv('ERA5_NUM_WORKERS')), ...
        str2double(getenv('SLURM_CPUS_PER_TASK'))];
    values = values(isfinite(values) & values >= 1);
    if isempty(values), n = fallback; else, n = floor(values(1)); end
end

function value = optimizer_settings(nStarts, maxIterations, maxEvaluations)
    value.display = 'off';
    value.numStarts = nStarts;
    value.maxIterations = maxIterations;
    value.maxFunctionEvaluations = maxEvaluations;
    value.tolX = 3e-4;
    value.tolFun = 3e-4;
    value.randomStartScale = 1.25;
    value.nearBestObjectiveTolerance = 1e-3;
    value.selectivePolishingEnabled = true;
    value.polishMaxIterations = 800;
    value.polishMaxFunctionEvaluations = 4000;
    value.polishTolX = 1e-4;
    value.polishTolFun = 1e-4;
    value.polishAcceptObjectiveTolerance = 1e-4;
end
