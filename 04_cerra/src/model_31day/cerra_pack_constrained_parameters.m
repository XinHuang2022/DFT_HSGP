function z = cerra_pack_constrained_parameters(theta, cfg, ...
        fixedEpsilon, fixedNoiseVariance)
%CERRA_PACK_CONSTRAINED_PARAMETERS Pack only parameters that are free.

% The order is signal variance, alpha, beta, optional spatial-scale
% contrast, and optional nugget variance.  Using one common mapping avoids
% silently changing the model when a nugget profile fixes sigma_n^2.

    if nargin < 3
        fixedEpsilon = [];
    end
    if nargin < 4
        fixedNoiseVariance = [];
    end
    values = zeros(3 + isempty(fixedEpsilon) + ...
        isempty(fixedNoiseVariance), 1);
    cursor = 0;
    cursor = cursor + 1;
    values(cursor) = inverse_bounded(theta.signalVariance, ...
        cfg.bounds.signalVariance, true);
    cursor = cursor + 1;
    values(cursor) = inverse_bounded(theta.alpha, ...
        cfg.bounds.alpha, true);
    cursor = cursor + 1;
    values(cursor) = inverse_bounded(theta.beta, ...
        cfg.bounds.beta, true);
    if isempty(fixedEpsilon)
        cursor = cursor + 1;
        delta = -0.5 * log(theta.epsilon);
        values(cursor) = inverse_bounded(delta, ...
            cfg.bounds.logSpatialScaleContrast, false);
    end
    if isempty(fixedNoiseVariance)
        cursor = cursor + 1;
        values(cursor) = inverse_bounded(theta.noiseVariance, ...
            cfg.bounds.noiseVariance, true);
    end
    z = values;
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
    probability = (value - lo) / (hi - lo);
    probability = min(1 - 1e-10, max(1e-10, probability));
    z = log(probability / (1 - probability));
end


