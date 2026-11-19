function fit = fit_pacific_eq9_rff(Y, xyCovariance, design, cfg, bank, ...
        initialTheta, fixedEpsilon, fixedNoiseVariance, ...
        optimizerOverride, seed)
%FIT_PACIFIC_EQ9_RFF Optimize the fixed-bank RFF Whittle REML.

    if nargin < 6 || isempty(initialTheta)
        initialTheta = cfg.start;
        trainingVariance = var(Y, 0, 'all');
        initialTheta.signalVariance = min(cfg.bounds.signalVariance(2), ...
            max(cfg.bounds.signalVariance(1), 0.8 * trainingVariance));
        initialTheta.noiseVariance = min(cfg.bounds.noiseVariance(2), ...
            max(cfg.bounds.noiseVariance(1), 0.2 * trainingVariance));
    end
    if nargin < 7
        fixedEpsilon = [];
    end
    if nargin < 8
        fixedNoiseVariance = [];
    end
    if nargin < 9 || isempty(optimizerOverride)
        optimizerOverride = cfg.optimizer;
    end
    if nargin < 10 || isempty(seed)
        seed = cfg.seed.optimizer;
    end
    if ~isempty(fixedEpsilon)
        initialTheta.epsilon = fixedEpsilon;
    end
    if ~isempty(fixedNoiseVariance)
        initialTheta.noiseVariance = fixedNoiseVariance;
    end
    initialTheta.nu = cfg.model.nu;
    prepared = prepare_pacific_whittle_data(Y, xyCovariance, design);
    prepared.rffPhi = pacific_rff_features(xyCovariance, bank);
    prepared.rffFeatureGram = prepared.rffPhi' * prepared.rffPhi;
    objective = @(z) pacific_profiled_whittle_reml_rff( ...
        z, prepared, cfg, bank, fixedEpsilon, fixedNoiseVariance);
    baseZ = pacific_pack_constrained_parameters(initialTheta, cfg, ...
        fixedEpsilon, fixedNoiseVariance);
    optimization = eq9_multistart_optimize( ...
        objective, baseZ, optimizerOverride, seed);
    [objectiveValue, details] = objective_with_details( ...
        optimization.z, prepared, cfg, bank, fixedEpsilon, ...
        fixedNoiseVariance);
    if isfield(details, 'errorMessage')
        error('Best Pacific RFF fit is invalid: %s', details.errorMessage);
    end
    theta = pacific_unpack_constrained_parameters(optimization.z, cfg, ...
        fixedEpsilon, fixedNoiseVariance);
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
    fit.bank = bank;
    fit.derived = eq9_derived_parameters(theta, cfg);
    fit.derived.timeScaleHours = cfg.data.samplingHours / theta.beta;
end

function [value, details] = objective_with_details(z, prepared, cfg, ...
        bank, fixedEpsilon, fixedNoiseVariance)
    [value, details] = pacific_profiled_whittle_reml_rff( ...
        z, prepared, cfg, bank, fixedEpsilon, fixedNoiseVariance);
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
