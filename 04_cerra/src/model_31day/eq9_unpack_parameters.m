function theta = eq9_unpack_parameters(z, cfg)
%EQ9_UNPACK_PARAMETERS Map unconstrained coordinates to model parameters.

    if numel(z) ~= 5
        error('The transformed parameter vector must have length five.');
    end
    z = z(:);
    theta.signalVariance = bounded_logistic(z(1), ...
        cfg.bounds.signalVariance, true);
    theta.alpha = bounded_logistic(z(2), cfg.bounds.alpha, true);
    theta.beta = bounded_logistic(z(3), cfg.bounds.beta, true);
    theta.logSpatialScaleContrast = bounded_logistic(z(4), ...
        cfg.bounds.logSpatialScaleContrast, false);
    theta.epsilon = exp(-2 * theta.logSpatialScaleContrast);
    theta.noiseVariance = bounded_logistic(z(5), ...
        cfg.bounds.noiseVariance, true);
    theta.nu = cfg.model.nu;
end

function value = bounded_logistic(z, bounds, useLogScale)
    probability = stable_logistic(z);
    if useLogScale
        lo = log(bounds(1));
        hi = log(bounds(2));
        value = exp(lo + probability * (hi - lo));
    else
        value = bounds(1) + probability * (bounds(2) - bounds(1));
    end
end

function p = stable_logistic(z)
    if z >= 0
        p = 1 / (1 + exp(-z));
    else
        ez = exp(z);
        p = ez / (1 + ez);
    end
end


