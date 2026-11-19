function [nll, details] = pacific_profiled_whittle_reml_rff( ...
        z, prepared, cfg, bank, fixedEpsilon, fixedNoiseVariance)
%PACIFIC_PROFILED_WHITTLE_REML_RFF Transformed-parameter RFF objective.

    if nargin < 5
        fixedEpsilon = [];
    end
    if nargin < 6
        fixedNoiseVariance = [];
    end
    details = struct();
    try
        theta = pacific_unpack_constrained_parameters(z, cfg, ...
            fixedEpsilon, fixedNoiseVariance);
        [nll, details] = pacific_profiled_whittle_theta_rff( ...
            theta, prepared, cfg, bank);
    catch exception
        nll = realmax('double') / 1e100;
        if nargout > 1
            details.errorIdentifier = exception.identifier;
            details.errorMessage = exception.message;
        end
    end
end
