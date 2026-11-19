function fit = fit_cerra_eq9_hsgp(Y, xyCovariance, design, cfg, ...
        basis, initialTheta, fixedEpsilon, fixedNoiseVariance, ...
        optimizerOverride, seed)
%FIT_CERRA_EQ9_HSGP Fit the deterministic-basis HSGP Whittle REML.

    if nargin < 6 || isempty(initialTheta)
        initialTheta = cfg.start;
        trainingVariance = var(Y, 0, 'all');
        initialTheta.signalVariance = interior_value( ...
            0.8 * trainingVariance, cfg.bounds.signalVariance);
        initialTheta.noiseVariance = interior_value( ...
            0.2 * trainingVariance, cfg.bounds.noiseVariance);
    end
    if nargin < 7, fixedEpsilon = []; end
    if nargin < 8, fixedNoiseVariance = []; end
    if nargin < 9 || isempty(optimizerOverride)
        optimizerOverride = cfg.hsgp.optimizer;
    end
    if nargin < 10 || isempty(seed), seed = cfg.seed.optimizer; end
    if ~isempty(fixedEpsilon), initialTheta.epsilon = fixedEpsilon; end
    if ~isempty(fixedNoiseVariance)
        initialTheta.noiseVariance = fixedNoiseVariance;
    end
    initialTheta.nu = cfg.model.nu;

    prepared = prepare_cerra_whittle_data(Y, xyCovariance, design);
    prepared.hsgpPhi = cerra_hsgp_features(xyCovariance, basis);
    prepared.hsgpFeatureGram = prepared.hsgpPhi' * prepared.hsgpPhi;
    objective = @(z) cerra_profiled_whittle_reml_hsgp( ...
        z, prepared, cfg, basis, fixedEpsilon, fixedNoiseVariance);
    baseZ = cerra_pack_constrained_parameters(initialTheta, cfg, ...
        fixedEpsilon, fixedNoiseVariance);
    optimization = eq9_multistart_optimize( ...
        objective, baseZ, optimizerOverride, seed);
    [objectiveValue, details] = objective(optimization.z);
    if isfield(details, 'errorMessage')
        error('Best CERRA HSGP fit is invalid: %s', details.errorMessage);
    end
    theta = cerra_unpack_constrained_parameters(optimization.z, cfg, ...
        fixedEpsilon, fixedNoiseVariance);
    [boundaryHit, minimumBoundaryDistance] = ...
        eq9_parameter_boundary_hit(theta, cfg);

    fit.theta = theta;
    fit.eta = details.eta;
    fit.covarianceEta = details.covarianceEta;
    fit.objective = objectiveValue;
    fit.profileNll = details.profileNll;
    fit.remlAdjustment = details.remlAdjustment;
    fit.details = details;
    fit.exitFlag = optimization.exitFlag;
    fit.convergenceCertified = optimization.convergenceCertified;
    fit.elapsedSeconds = optimization.elapsedSeconds;
    fit.totalFunctionEvaluations = function_count(optimization);
    fit.totalIterations = iteration_count(optimization);
    fit.optimization = optimization;
    fit.fixedEpsilon = fixedEpsilon;
    fit.fixedNoiseVariance = fixedNoiseVariance;
    fit.boundaryHit = boundaryHit;
    fit.minimumBoundaryDistance = minimumBoundaryDistance;
    fit.basis = basis;
    fit.derived = eq9_derived_parameters(theta, cfg);
    fit.derived.timeScaleHours = cfg.data.samplingHours / theta.beta;
    fit.derived.spatialRange10LowFrequencyKm = ...
        fit.derived.spatialRange10LowFrequency * ...
        cfg.spatial.coordinateScaleKm;
    fit.derived.spatialRange10HighFrequencyKm = ...
        fit.derived.spatialRange10HighFrequency * ...
        cfg.spatial.coordinateScaleKm;
end

function value = interior_value(value, bounds)
    value = min(bounds(2) / (1 + 1e-8), ...
        max(bounds(1) * (1 + 1e-8), value));
end

function count = function_count(optimization)
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

function count = iteration_count(optimization)
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


