function covariance = era5_hsgp_covariance(PhiA, PhiB, weight, ...
        nuggetVariance, sameSites)
%ERA5_HSGP_COVARIANCE Form PhiA*diag(weight)*PhiB' plus a nugget.

    if nargin < 4 || isempty(nuggetVariance), nuggetVariance = 0; end
    if nargin < 5, sameSites = false; end
    weight = weight(:);
    if size(PhiA, 2) ~= numel(weight) || ...
            size(PhiB, 2) ~= numel(weight)
        error('HSGP feature and weight dimensions are inconsistent.');
    end
    covariance = (PhiA .* weight.') * PhiB.';
    if sameSites
        covariance = real((covariance + covariance') / 2);
        if nuggetVariance > 0
            covariance = covariance + ...
                nuggetVariance * eye(size(PhiA, 1));
        end
    elseif nuggetVariance ~= 0
        error('A nugget can be added only for a same-site covariance.');
    end
end

