function cfg = pacific_eq9_config(preset)
%PACIFIC_EQ9_CONFIG Configuration for Pacific Equation-9 CV experiments.

    if nargin < 1 || isempty(preset)
        preset = 'audit';
    end
    preset = lower(char(preset));
    validPresets = {'audit', 'smoke3', 'pilot5', ...
        'pilot5_all', 'local5_cv5', 'local7_cv5', ...
        'full17_cv5', 'sparse3'};
    if ~ismember(preset, validPresets)
        error(['Unknown preset "%s". Use audit, smoke3, pilot5, ', ...
            'pilot5_all, local5_cv5, local7_cv5, full17_cv5, ', ...
            'or sparse3.'], preset);
    end

    suiteRoot = fileparts(mfilename('fullpath'));
    experimentRoot = fileparts(fileparts(suiteRoot));
    resultsRoot = getenv('PACIFIC_REPRO_OUTPUT_DIR');
    if isempty(resultsRoot)
        resultsRoot = fullfile(experimentRoot, 'results', 'current');
    end
    cfg.version = '0.2.0';
    cfg.preset = preset;
    cfg.suiteRoot = suiteRoot;
    cfg.paths.dataDir = fullfile(experimentRoot, 'input', 'data');
    cfg.paths.resultsRoot = resultsRoot;
    cfg.paths.outputDir = fullfile(resultsRoot, preset);
    cfg.paths.figureDir = fullfile(cfg.paths.outputDir, 'figures');
    cfg.paths.windFileName = 'NCEP_U.dat';
    cfg.paths.coordinateFileName = 'NCEPlatlon.dat';

    %% Published data contract
    cfg.data.expectedTimePoints = 480;
    cfg.data.expectedStations = 289;
    cfg.data.expectedGridSide = 17;
    cfg.data.samplingHours = 6;
    cfg.data.responseName = 'east_west_u_wind';
    cfg.data.responseUnits = 'm/s';
    % Values beyond this deliberately generous physical limit are treated
    % as possible encoded missing-value sentinels, not as valid wind.
    cfg.data.absolutePlausibilityLimitMps = 150;

    %% Equation-9 covariance
    % Spatial covariance coordinates are measured in units of 100 km.
    % One discrete time step is six hours; beta is therefore per six-hour
    % step and the reported time scale is converted back to hours.
    cfg.spatial.coordinateScaleKm = 100;
    cfg.model.dimension = 2;
    cfg.model.nu = 2.5;
    cfg.model.deltaTime = 1;
    cfg.model.aliasTerms = 6;
    cfg.model.normalizationGridSize = 2048;

    cfg.bounds.signalVariance = [1e-3, 100];
    cfg.bounds.alpha = [0.003, 80];
    cfg.bounds.beta = [0.003, 3.0];
    % The free model approaches, but never numerically equals, the
    % separable limit.  The tight upper bound makes the embedded separable
    % start effectively nested even for the 480-point Pacific record.
    cfg.bounds.epsilon = [0.02, 1 - 1e-6];
    cfg.bounds.logSpatialScaleContrast = ...
        [-0.5 * log(cfg.bounds.epsilon(2)), ...
         -0.5 * log(cfg.bounds.epsilon(1))];
    cfg.bounds.noiseVariance = [1e-5, 50];
    cfg.parameterization = 'log_spatial_scale_contrast';

    % The variance starts are replaced fold-by-fold using the training data.
    cfg.start.signalVariance = 4.0;
    cfg.start.alpha = 0.35;
    cfg.start.beta = 0.25;
    cfg.start.epsilon = 0.55;
    cfg.start.noiseVariance = 0.4;
    cfg.start.nu = cfg.model.nu;
    cfg.likelihood.type = 'whittle_reml';

    %% Mean, validation, and diagnostics
    cfg.mean.coefficientNames = ["intercept"; "east"; "north"];
    cfg.cv.nominalCoverage = 0.95;
    cfg.cv.minimumPredictiveVariance = 1e-10;
    cfg.cv.maximumAllowedImaginary = 1e-8;
    cfg.cv.timeBlocks = 4;
    cfg.cv.savePredictions = true;
    cfg.cv.saveVarianceDecomposition = true;
    cfg.audit.maximumAcfLag = 40;
    cfg.figures.visible = 'off';

    %% Reproducible multistart optimization
    cfg.seed.optimizer = 83101;
    cfg.optimizer.display = 'off';
    cfg.optimizer.tolX = 3e-4;
    cfg.optimizer.tolFun = 3e-4;
    cfg.optimizer.randomStartScale = 1.25;
    cfg.optimizer.nearBestObjectiveTolerance = 1e-3;
    cfg.optimizer.selectivePolishingEnabled = true;
    cfg.optimizer.polishMaxIterations = 800;
    cfg.optimizer.polishMaxFunctionEvaluations = 4000;
    cfg.optimizer.polishTolX = 1e-4;
    cfg.optimizer.polishTolFun = 1e-4;
    cfg.optimizer.polishAcceptObjectiveTolerance = 1e-4;

    %% Outer-fold parallelism only
    cfg.parallel.useParallel = false;
    cfg.parallel.autoStartPool = true;
    cfg.parallel.numWorkers = 5;

    %% Dense full-grid safety gate and HPC metadata
    % The exact full-grid likelihood is intentionally opt-in because each
    % objective evaluation factors hundreds of dense spatial matrices.
    cfg.safety.requireFullGridOptIn = true;
    cfg.safety.fullGridOptIn = false;
    cfg.hpc.recommendedMemoryGb = 64;
    cfg.hpc.maximumWallHours = 168;

    %% Preset definitions
    cfg.runModel = ~strcmp(preset, 'audit');
    switch preset
        case 'audit'
            cfg.subset.name = 'local3';
            cfg.cv.holdoutPattern = 'none';
            cfg.optimizer.numStarts = 1;
            cfg.optimizer.maxIterations = 1;
            cfg.optimizer.maxFunctionEvaluations = 1;
        case 'smoke3'
            cfg.subset.name = 'local3';
            cfg.cv.holdoutPattern = 'center';
            cfg.optimizer.numStarts = 2;
            cfg.optimizer.maxIterations = 180;
            cfg.optimizer.maxFunctionEvaluations = 900;
            cfg.model.normalizationGridSize = 1024;
        case 'pilot5'
            cfg.subset.name = 'local5';
            cfg.cv.holdoutPattern = 'representative5';
            cfg.optimizer.numStarts = 4;
            cfg.optimizer.maxIterations = 350;
            cfg.optimizer.maxFunctionEvaluations = 2000;
        case 'pilot5_all'
            cfg.subset.name = 'local5';
            cfg.cv.holdoutPattern = 'all';
            cfg.optimizer.numStarts = 4;
            cfg.optimizer.maxIterations = 350;
            cfg.optimizer.maxFunctionEvaluations = 2000;
        case 'local5_cv5'
            cfg.subset.name = 'local5';
            cfg.cv.holdoutPattern = 'balanced5';
            cfg.parallel.numWorkers = 5;
            cfg.optimizer.numStarts = 4;
            cfg.optimizer.maxIterations = 350;
            cfg.optimizer.maxFunctionEvaluations = 2000;
        case 'local7_cv5'
            cfg.subset.name = 'local7';
            cfg.cv.holdoutPattern = 'balanced5';
            cfg.parallel.numWorkers = 5;
            cfg.optimizer.numStarts = 4;
            cfg.optimizer.maxIterations = 350;
            cfg.optimizer.maxFunctionEvaluations = 2000;
        case 'full17_cv5'
            cfg.subset.name = 'full17';
            cfg.cv.holdoutPattern = 'balanced5';
            cfg.parallel.numWorkers = 5;
            % Three starts retain the embedded separable start while
            % limiting the cost of this deliberately expensive stage.
            cfg.optimizer.numStarts = 3;
            cfg.optimizer.maxIterations = 300;
            cfg.optimizer.maxFunctionEvaluations = 1800;
        case 'sparse3'
            cfg.subset.name = 'sparse3';
            cfg.cv.holdoutPattern = 'center';
            cfg.optimizer.numStarts = 2;
            cfg.optimizer.maxIterations = 180;
            cfg.optimizer.maxFunctionEvaluations = 900;
            cfg.model.normalizationGridSize = 1024;
    end
end
