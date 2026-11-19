function [F, details] = eq9_frequency_covariances(theta, xy, lambda, cfg, includeNugget)
%EQ9_FREQUENCY_COVARIANCES Spatial covariance at temporal frequencies.
%
% F(:,:,k) is the covariance of the station DFT vector at lambda(k).

    if nargin < 5
        includeNugget = true;
    end
    if any(~isfinite(xy), 'all') || size(xy, 2) ~= cfg.model.dimension
        error('xy must be a finite m-by-%d coordinate array.', ...
            cfg.model.dimension);
    end
    distances = eq9_pairwise_distances(xy, xy);
    [Hraw, Sraw, kappa] = eq9_raw_half_spectrum( ...
        distances, lambda, theta, cfg);
    normalizer = eq9_spectral_normalizer(theta, cfg);
    F = theta.signalVariance * normalizer * Hraw;

    if includeNugget
        identity = eye(size(xy, 1));
        for k = 1:numel(lambda)
            F(:, :, k) = F(:, :, k) + theta.noiseVariance * identity;
        end
    end

    details.normalizer = normalizer;
    details.rawTemporalPower = Sraw;
    details.kappaByAlias = kappa;
    details.distances = distances;
end



