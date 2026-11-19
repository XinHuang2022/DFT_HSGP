function fit = fit_pacific_eq9_model(Y, xyCovariance, design, cfg, ...
        fixedEpsilon, optimizerOverride, seed, fixedNoiseVariance, ...
        initialTheta)
%FIT_PACIFIC_EQ9_MODEL Profiled Whittle REML with optional constraints.
%
% fixedEpsilon and fixedNoiseVariance are empty when the corresponding
% parameter is free.  initialTheta supplies a deterministic warm start; it
% changes only the numerical search, never the likelihood being optimized.

    if nargin < 5
        fixedEpsilon = [];
    end
    if nargin < 6 || isempty(optimizerOverride)
        optimizerOverride = cfg.optimizer;
    end
    if nargin < 7 || isempty(seed)
        seed = cfg.seed.optimizer;
    end
    if nargin < 8
        fixedNoiseVariance = [];
    end
    if nargin < 9
        initialTheta = [];
    end
    trainingVariance = var(Y, 0, 'all');
    if isempty(initialTheta)
        startTheta = cfg.start;
        startTheta.signalVariance = interior_value( ...
            0.8 * trainingVariance, cfg.bounds.signalVariance);
        startTheta.noiseVariance = interior_value( ...
            0.2 * trainingVariance, cfg.bounds.noiseVariance);
    else
        startTheta = initialTheta;
        startTheta.nu = cfg.model.nu;
    end
    if ~isempty(fixedEpsilon)
        startTheta.epsilon = fixedEpsilon;
    end
    if ~isempty(fixedNoiseVariance)
        startTheta.noiseVariance = fixedNoiseVariance;
    end
    prepared = prepare_pacific_whittle_data(Y, xyCovariance, design);
    objective = @(z) pacific_profiled_whittle_reml( ...
        z, prepared, cfg, fixedEpsilon, fixedNoiseVariance);
    baseZ = pacific_pack_constrained_parameters(startTheta, cfg, ...
        fixedEpsilon, fixedNoiseVariance);
    optimization = eq9_multistart_optimize( ...
        objective, baseZ, optimizerOverride, seed);
    [objectiveValue, details] = pacific_profiled_whittle_reml( ...
        optimization.z, prepared, cfg, fixedEpsilon, ...
        fixedNoiseVariance);
    if isfield(details, 'errorMessage')
        error('Best Pacific Equation-9 fit is invalid: %s', ...
            details.errorMessage);
    end
    theta = pacific_unpack_constrained_parameters(optimization.z, cfg, ...
        fixedEpsilon, fixedNoiseVariance);

    fit.theta = theta;
    fit.derived = eq9_derived_parameters(theta, cfg);
    fit.derived.spatialRange10LowFrequencyKm = ...
        fit.derived.spatialRange10LowFrequency * ...
        cfg.spatial.coordinateScaleKm;
    fit.derived.spatialRange10HighFrequencyKm = ...
        fit.derived.spatialRange10HighFrequency * ...
        cfg.spatial.coordinateScaleKm;
    fit.derived.timeScaleHours = cfg.data.samplingHours / theta.beta;
    fit.z = optimization.z;
    fit.eta = details.eta;
    fit.covarianceEta = details.covarianceEta;
    fit.objective = objectiveValue;
    fit.profileNll = details.profileNll;
    fit.remlAdjustment = details.remlAdjustment;
    fit.details = details;
    fit.exitFlag = optimization.exitFlag;
    fit.convergenceCertified = optimization.convergenceCertified;
    fit.elapsedSeconds = optimization.elapsedSeconds;
    fit.secondBestObjectiveGap = optimization.secondBestObjectiveGap;
    fit.nearBestStartCount = optimization.nearBestStartCount;
    fit.polishTriggered = optimization.polishTriggered;
    fit.polishAccepted = optimization.polishAccepted;
    fit.prePolishObjective = optimization.prePolishObjective;
    fit.polishObjective = optimization.polishObjective;
    fit.prePolishExitFlag = optimization.prePolishExitFlag;
    fit.polishExitFlag = optimization.polishExitFlag;
    fit.totalFunctionEvaluations = optimization_function_count( ...
        optimization);
    fit.totalIterations = optimization_iteration_count(optimization);
    fit.fixedEpsilon = fixedEpsilon;
    fit.fixedNoiseVariance = fixedNoiseVariance;
    if isempty(fixedEpsilon) && isempty(fixedNoiseVariance)
        [fit.boundaryHit, fit.minimumBoundaryDistance] = ...
            eq9_parameter_boundary_hit(theta, cfg);
    else
        [fit.boundaryHit, fit.minimumBoundaryDistance] = ...
            constrained_fit_boundary(theta, cfg, fixedEpsilon, ...
            fixedNoiseVariance);
    end
    fit.optimization = optimization;
end

function count = optimization_function_count(optimization)
    count = numel(optimization.objectiveAtStart);
    for j = 1:numel(optimization.outputs)
        if isfield(optimization.outputs{j}, 'funcCount')
            count = count + optimization.outputs{j}.funcCount;
        end
    end
    if optimization.polishTriggered && ...
            isfield(optimization.polishOutput, 'funcCount')
        count = count + optimization.polishOutput.funcCount;
    end
end

function count = optimization_iteration_count(optimization)
    count = 0;
    for j = 1:numel(optimization.outputs)
        if isfield(optimization.outputs{j}, 'iterations')
            count = count + optimization.outputs{j}.iterations;
        end
    end
    if optimization.polishTriggered && ...
            isfield(optimization.polishOutput, 'iterations')
        count = count + optimization.polishOutput.iterations;
    end
end

function value = interior_value(value, bounds)
    value = min(bounds(2) / 1.01, max(bounds(1) * 1.01, value));
end

function [hit, minimumDistance] = constrained_fit_boundary(theta, cfg, ...
        fixedEpsilon, fixedNoiseVariance)
    names = {'signalVariance', 'alpha', 'beta'};
    if isempty(fixedEpsilon)
        names{end + 1} = 'epsilon';
    end
    if isempty(fixedNoiseVariance)
        names{end + 1} = 'noiseVariance';
    end
    distance = zeros(numel(names), 1);
    for j = 1:numel(names)
        if strcmp(names{j}, 'epsilon')
            bounds = cfg.bounds.logSpatialScaleContrast;
            value = -0.5 * log(theta.epsilon);
            fraction = (value - bounds(1)) / (bounds(2) - bounds(1));
        else
            bounds = cfg.bounds.(names{j});
            fraction = (log(theta.(names{j})) - log(bounds(1))) / ...
                (log(bounds(2)) - log(bounds(1)));
        end
        distance(j) = min(fraction, 1 - fraction);
    end
    minimumDistance = min(distance);
    hit = minimumDistance < 0.01;
end
