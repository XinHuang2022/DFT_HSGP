function cfg = pacific_scalable_config(basePreset)
%PACIFIC_SCALABLE_CONFIG Configuration for nugget, subset, and RFF stages.

% The scalable stages retain the v0.2 Equation-9 model and outer folds.
% They change only how covariance parameters are estimated and how the
% spatial covariance matrices are evaluated.

    if nargin < 1 || isempty(basePreset)
        basePreset = 'full17_cv5';
    end
    cfg = pacific_eq9_config(basePreset);
    cfg.version = '0.3.0';
    cfg.safety.fullGridOptIn = false;

    %% Spatially balanced covariance-parameter subset
    cfg.scalable.parameterSubsetSize = 49;
    cfg.scalable.localPairFraction = 0.25;
    cfg.scalable.subsetSeed = 51917;
    cfg.scalable.nSubsetSeedsForSensitivity = 3;
    cfg.scalable.refitMeanOnAllOuterTrainingSites = true;

    %% Nested, training-only nugget study
    % Fractions multiply the empirical residual variance after removing an
    % OLS spatial plane.  This makes the grid invariant to wind units.
    cfg.nugget.innerValidationFoldCount = 5;
    cfg.nugget.innerValidationFold = 1;
    cfg.nugget.residualVarianceFractions = ...
        [1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1];
    cfg.nugget.approximateConfidenceMultiplier = 2;
    cfg.nugget.useDecisionFile = true;
    cfg.nugget.decisionFile = fullfile(cfg.paths.resultsRoot, ...
        [char(basePreset), '_nugget_study'], 'nugget_decision.csv');
    cfg.nugget.optimizer = cfg.optimizer;
    cfg.nugget.optimizer.numStarts = 3;
    cfg.nugget.optimizer.maxIterations = 300;
    cfg.nugget.optimizer.maxFunctionEvaluations = 1800;

    %% Fixed-bank spatial RFF
    cfg.rff.featureCounts = [32, 64, 96, 128, 192, 256];
    cfg.rff.seeds = [12011, 12029, 12041];
    cfg.rff.proposalMixtureWeights = [1/3, 1/3, 1/3];
    cfg.rff.temporalBandEdges = [0, 1/6, 1/3, 1/2, 1];
    cfg.rff.maximumConditionNumber = 1e12;
    cfg.rff.minimumEssFraction = 0.10;
    % Gates compare RFF with exact plug-in calculations at identical theta.
    % They do not use outer-target responses, so the later CV scores remain
    % an honest assessment rather than a feature-count selection score.
    cfg.rff.acceptance.maximumMedianCovarianceError = 0.08;
    cfg.rff.acceptance.maximumPredictionMeanRmseMps = 0.02;
    cfg.rff.acceptance.maximumMedianObservationVarianceRelativeError = ...
        0.10;
    cfg.rff.acceptance.maximumMedianConditionalVarianceRelativeError = ...
        0.10;
    % A Whittle-objective difference of 0.05 nats per scalar observation
    % corresponds to a geometric likelihood-score ratio exp(0.05)=1.051.
    % This gate is essential when the exact nugget is very small: a modest
    % Frobenius covariance error can otherwise hide large errors in inverse
    % covariance directions and in the log determinant.
    cfg.rff.acceptance.maximumAbsoluteObjectiveDifferencePerObservation = ...
        0.05;
    cfg.rff.acceptance.maximumImaginaryPrediction = 1e-8;
    cfg.rff.acceptance.minimumPassingSeedFraction = 2/3;
    % Accuracy and acceleration are separate gates.  The arithmetic ratio
    % below compares one p-space factorization with one dense m-space
    % Cholesky; measured objective timings are reported alongside it.
    cfg.rff.acceptance.maximumEstimatedFlopRatioForAcceleration = 1;
    cfg.rff.fitEnabledOnlyAfterStage0Pass = true;
    cfg.rff.requireAccelerationBeforeFit = true;
    cfg.rff.fitFeatureSeeds = cfg.rff.seeds;
    cfg.rff.fixNuggetAtExactSubsetEstimateDuringFirstFit = true;
    cfg.rff.optimizer = cfg.optimizer;
    cfg.rff.optimizer.numStarts = 3;
    cfg.rff.optimizer.maxIterations = 300;
    cfg.rff.optimizer.maxFunctionEvaluations = 1800;

    cfg.paths.nuggetDir = fullfile(cfg.paths.resultsRoot, ...
        [char(basePreset), '_nugget_study']);
    cfg.paths.subsetPluginDir = fullfile(cfg.paths.resultsRoot, ...
        [char(basePreset), '_subset_plugin']);
    cfg.paths.rffStage0Dir = fullfile(cfg.paths.resultsRoot, ...
        [char(basePreset), '_rff_stage0']);
    cfg.paths.rffFitDir = fullfile(cfg.paths.resultsRoot, ...
        [char(basePreset), '_rff_fit']);
end
