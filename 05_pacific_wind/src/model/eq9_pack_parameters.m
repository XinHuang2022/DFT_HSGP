function z = eq9_pack_parameters(theta, cfg)
%EQ9_PACK_PARAMETERS Map model parameters to unconstrained coordinates.

    z = zeros(5, 1);
    z(1) = inverse_bounded(theta.signalVariance, ...
        cfg.bounds.signalVariance, true);
    z(2) = inverse_bounded(theta.alpha, cfg.bounds.alpha, true);
    z(3) = inverse_bounded(theta.beta, cfg.bounds.beta, true);
    logSpatialScaleContrast = -0.5 * log(theta.epsilon);
    z(4) = inverse_bounded(logSpatialScaleContrast, ...
        cfg.bounds.logSpatialScaleContrast, false);
    z(5) = inverse_bounded(theta.noiseVariance, ...
        cfg.bounds.noiseVariance, true);
end

function z = inverse_bounded(value, bounds, useLogScale)
    if useLogScale
        value = log(value);
        lo = log(bounds(1));
        hi = log(bounds(2));
    else
        lo = bounds(1);
        hi = bounds(2);
    end
    p = (value - lo) / (hi - lo);
    p = min(1 - 1e-10, max(1e-10, p));
    z = log(p / (1 - p));
end
