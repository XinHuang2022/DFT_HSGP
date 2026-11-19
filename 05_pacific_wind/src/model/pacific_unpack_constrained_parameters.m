function theta = pacific_unpack_constrained_parameters(z, cfg, ...
        fixedEpsilon, fixedNoiseVariance)
%PACIFIC_UNPACK_CONSTRAINED_PARAMETERS Unpack free and fixed parameters.

    if nargin < 3
        fixedEpsilon = [];
    end
    if nargin < 4
        fixedNoiseVariance = [];
    end
    expected = 3 + isempty(fixedEpsilon) + isempty(fixedNoiseVariance);
    if numel(z) ~= expected
        error('Constrained parameter vector has length %d; expected %d.', ...
            numel(z), expected);
    end
    z = z(:);
    cursor = 0;
    cursor = cursor + 1;
    theta.signalVariance = bounded(z(cursor), ...
        cfg.bounds.signalVariance, true);
    cursor = cursor + 1;
    theta.alpha = bounded(z(cursor), cfg.bounds.alpha, true);
    cursor = cursor + 1;
    theta.beta = bounded(z(cursor), cfg.bounds.beta, true);
    if isempty(fixedEpsilon)
        cursor = cursor + 1;
        theta.logSpatialScaleContrast = bounded(z(cursor), ...
            cfg.bounds.logSpatialScaleContrast, false);
        theta.epsilon = exp(-2 * theta.logSpatialScaleContrast);
    else
        theta.epsilon = fixedEpsilon;
        theta.logSpatialScaleContrast = -0.5 * log(fixedEpsilon);
    end
    if isempty(fixedNoiseVariance)
        cursor = cursor + 1;
        theta.noiseVariance = bounded(z(cursor), ...
            cfg.bounds.noiseVariance, true);
    else
        theta.noiseVariance = fixedNoiseVariance;
    end
    theta.nu = cfg.model.nu;
end

function value = bounded(z, bounds, useLogScale)
    if z >= 0
        probability = 1 / (1 + exp(-z));
    else
        ez = exp(z);
        probability = ez / (1 + ez);
    end
    if useLogScale
        value = exp(log(bounds(1)) + probability * ...
            (log(bounds(2)) - log(bounds(1))));
    else
        value = bounds(1) + probability * (bounds(2) - bounds(1));
    end
end
