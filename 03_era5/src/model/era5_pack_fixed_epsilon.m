function z = era5_pack_fixed_epsilon(theta, cfg)
%ERA5_PACK_FIXED_EPSILON Pack four free parameters for epsilon fixed.

    z = zeros(4, 1);
    z(1) = inverse_bounded(theta.signalVariance, ...
        cfg.bounds.signalVariance);
    z(2) = inverse_bounded(theta.alpha, cfg.bounds.alpha);
    z(3) = inverse_bounded(theta.beta, cfg.bounds.beta);
    z(4) = inverse_bounded(theta.noiseVariance, ...
        cfg.bounds.noiseVariance);
end

function z = inverse_bounded(value, bounds)
    value = log(value);
    lo = log(bounds(1));
    hi = log(bounds(2));
    probability = (value - lo) / (hi - lo);
    probability = min(1 - 1e-10, max(1e-10, probability));
    z = log(probability / (1 - probability));
end

