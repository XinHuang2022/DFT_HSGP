function [nll, details] = cerra_profiled_whittle_reml( ...
        z, prepared, cfg, fixedEpsilon, fixedNoiseVariance)
%CERRA_PROFILED_WHITTLE_REML Profile spatial-plane mean and evaluate REML.

    if nargin < 4
        fixedEpsilon = [];
    end
    if nargin < 5
        fixedNoiseVariance = [];
    end
    details = struct();
    try
        theta = cerra_unpack_constrained_parameters(z, cfg, ...
            fixedEpsilon, fixedNoiseVariance);
        [nll, details] = cerra_profiled_whittle_theta( ...
            theta, prepared, cfg);
    catch exception
        nll = realmax('double') / 1e100;
        if nargout > 1
            details.errorIdentifier = exception.identifier;
            details.errorMessage = exception.message;
        end
    end
end


