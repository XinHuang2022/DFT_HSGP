function [nll, details] = pacific_profiled_whittle_reml( ...
        z, prepared, cfg, fixedEpsilon, fixedNoiseVariance)
%PACIFIC_PROFILED_WHITTLE_REML Profile spatial-plane mean and evaluate REML.

    if nargin < 4
        fixedEpsilon = [];
    end
    if nargin < 5
        fixedNoiseVariance = [];
    end
    details = struct();
    try
        theta = pacific_unpack_constrained_parameters(z, cfg, ...
            fixedEpsilon, fixedNoiseVariance);
        [nll, details] = pacific_profiled_whittle_theta( ...
            theta, prepared, cfg);
    catch exception
        nll = realmax('double') / 1e100;
        if nargout > 1
            details.errorIdentifier = exception.identifier;
            details.errorMessage = exception.message;
        end
    end
end
