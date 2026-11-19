function [varianceReference, etaOls, residual] = ...
        pacific_empirical_residual_variance(Y, xyKm)
%PACIFIC_EMPIRICAL_RESIDUAL_VARIANCE Scale-free nugget reference variance.

% The spatial plane is estimated by ordinary least squares from station
% time averages.  The returned variance is the mean squared residual over
% all time-site observations and is used only to express candidate nuggets
% as fractions of a fold-specific data scale.

    center = mean(xyKm, 1);
    scale = std(xyKm, 0, 1);
    scale(scale <= sqrt(eps)) = 1;
    X = [ones(size(xyKm, 1), 1), (xyKm - center) ./ scale];
    etaOls = X \ mean(Y, 1).';
    fitted = X * etaOls;
    residual = Y - repmat(fitted.', size(Y, 1), 1);
    varianceReference = mean(residual .^ 2, 'all');
    if ~(isfinite(varianceReference) && varianceReference > 0)
        error('Empirical residual variance is not positive and finite.');
    end
end
