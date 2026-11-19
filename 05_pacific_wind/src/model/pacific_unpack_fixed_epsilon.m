function theta = pacific_unpack_fixed_epsilon(z, epsilon, cfg)
%PACIFIC_UNPACK_FIXED_EPSILON Unpack four free covariance parameters.

    if numel(z) ~= 4
        error('Fixed-epsilon transformed parameter vector must have length 4.');
    end
    z = z(:);
    theta.signalVariance = bounded_logistic(z(1), ...
        cfg.bounds.signalVariance);
    theta.alpha = bounded_logistic(z(2), cfg.bounds.alpha);
    theta.beta = bounded_logistic(z(3), cfg.bounds.beta);
    theta.epsilon = epsilon;
    theta.logSpatialScaleContrast = -0.5 * log(epsilon);
    theta.noiseVariance = bounded_logistic(z(4), ...
        cfg.bounds.noiseVariance);
    theta.nu = cfg.model.nu;
end

function value = bounded_logistic(z, bounds)
    if z >= 0
        probability = 1 / (1 + exp(-z));
    else
        ez = exp(z);
        probability = ez / (1 + ez);
    end
    value = exp(log(bounds(1)) + probability * ...
        (log(bounds(2)) - log(bounds(1))));
end
