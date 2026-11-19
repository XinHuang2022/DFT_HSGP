function cfg = pacific_hsgp_config(basePreset)
%PACIFIC_HSGP_CONFIG Configuration for the Pacific HSGP pilot.

    if nargin < 1 || isempty(basePreset)
        basePreset = 'local7_cv5';
    end
    basePreset = lower(char(basePreset));
    if ~ismember(basePreset, {'local7_cv5', 'full17_cv5'})
        error('HSGP presets are local7_cv5 and full17_cv5.');
    end

    cfg = pacific_scalable_config(basePreset);
    cfg.version = '0.3.1-hsgp';
    cfg.safety.fullGridOptIn = true;
    cfg.parallel.useParallel = true;
    cfg.parallel.numWorkers = 5;

    cfg.paths.hsgpAuditDir = fullfile(cfg.paths.resultsRoot, ...
        [basePreset, '_hsgp_basis_audit']);
    cfg.paths.hsgpCvDir = fullfile(cfg.paths.resultsRoot, ...
        [basePreset, '_hsgp_cv']);

    % Geometry-only HSGP certification.  This development region was fixed
    % from the earlier Pacific exact/subset fits, before inspecting the new
    % HSGP CV target scores.  Every fitted value is checked against it.  A
    % fit outside the box invalidates the selected basis and requires a new,
    % expanded audit rather than silent extrapolation.
    cfg.hsgp.boundaryFactors = [4, 5, 6, 7, 8];
    cfg.hsgp.basisPerDimension = ...
        [8, 10, 12, 14, 16, 18, 20, 22, 24, 28, 32, 36, 40];
    cfg.hsgp.auditAlpha = [0.10, 0.13, 0.17];
    cfg.hsgp.auditBeta = [0.70, 1.00, 1.30];
    cfg.hsgp.auditEpsilon = [0.12, 0.20, 0.32];
    cfg.hsgp.auditSignalVariance = 1;
    cfg.hsgp.diagnosticFrequencies = 12;
    cfg.hsgp.maximumAuditSites = 81;
    cfg.hsgp.acceptance.maximumMedianLatentCovarianceError = 0.05;
    cfg.hsgp.acceptance.maximumLatentCovarianceError = 0.10;
    % Response-free weak-eigenvalue diagnostic.  This is evaluated after
    % adding the audit nugget and whitening by the exact covariance.
    cfg.hsgp.acceptance.maximumWhitenedCovarianceError = 0.20;
    cfg.hsgp.acceptance.maximumAbsoluteObjectiveDifferencePerObservation = ...
        0.05;
    cfg.hsgp.acceptance.minimumMedianDiagonalRatio = 0.95;
    cfg.hsgp.acceptance.maximumMedianDiagonalRatio = 1.05;
    cfg.hsgp.requireCertifiedBasis = true;
    cfg.hsgp.auditNoiseVariance = cfg.bounds.noiseVariance(1);

    % With a very small nugget, M<n leaves n-M directions supported almost
    % entirely by the nugget.  The automatic solver uses a QR/Woodbury
    % factorization only when M<n, and otherwise uses a direct site-space
    % Cholesky.  The selected basis is not required to be low rank in this
    % first scientific pilot.
    cfg.hsgp.factorization = 'auto';
    cfg.hsgp.requireDimensionReduction = false;
    cfg.hsgp.minimumReciprocalCondition = 1e-12;
    cfg.hsgp.maximumAllowedJitter = 1e-6;
    cfg.hsgp.minimumPredictiveVariance = cfg.cv.minimumPredictiveVariance;

    % The earlier nested training-only study supported only a very small
    % scale-adaptive floor.  It remains a lower bound, not a fixed nugget.
    cfg.hsgp.nuggetResidualVarianceFraction = 1e-6;
    cfg.hsgp.reportNuggetBoundaryTolerance = 0.02;
    cfg.hsgp.fixedNoiseVariance = [];

    cfg.hsgp.optimizer = cfg.optimizer;
    if strcmp(basePreset, 'local7_cv5')
        cfg.hsgp.optimizer.numStarts = 2;
        cfg.hsgp.optimizer.maxIterations = 220;
        cfg.hsgp.optimizer.maxFunctionEvaluations = 1100;
        cfg.hsgp.auditExactAtFittedTheta = true;
        cfg.hsgp.predictExactAtFittedTheta = true;
    else
        cfg.hsgp.basisPerDimension = ...
            [16, 20, 24, 28, 32, 40, 48, 56, 64];
        cfg.hsgp.optimizer.numStarts = 3;
        cfg.hsgp.optimizer.maxIterations = 300;
        cfg.hsgp.optimizer.maxFunctionEvaluations = 1800;
        % These are expensive dense audits, not required for HSGP CV.
        cfg.hsgp.auditExactAtFittedTheta = false;
        cfg.hsgp.predictExactAtFittedTheta = false;
    end
    cfg.hsgp.exampleStation = "";
    cfg.hsgp.bootstrapReplicates = 2000;
    cfg.hsgp.bootstrapSeed = 73021;
end
