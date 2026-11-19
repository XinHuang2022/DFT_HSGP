function [nll, details] = era5_profiled_whittle_reml_hsgp( ...
        z, prepared, cfg, basis, fixedEpsilon, fixedNoiseVariance)
%ERA5_PROFILED_WHITTLE_REML_HSGP Transformed HSGP objective.

    if nargin < 5, fixedEpsilon = []; end
    if nargin < 6, fixedNoiseVariance = []; end
    details = struct();
    try
        theta = era5_unpack_constrained_parameters(z, cfg, ...
            fixedEpsilon, fixedNoiseVariance);
        [nll, details] = era5_profiled_whittle_theta_hsgp( ...
            theta, prepared, cfg, basis);
    catch exception
        nll = realmax('double') / 1e100;
        details.errorIdentifier = exception.identifier;
        details.errorMessage = exception.message;
    end
end

