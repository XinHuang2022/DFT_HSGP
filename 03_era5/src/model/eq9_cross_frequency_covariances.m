function [Fcross, details] = eq9_cross_frequency_covariances( ...
    theta, xyA, xyB, lambda, cfg)
%EQ9_CROSS_FREQUENCY_COVARIANCES Latent cross-spectrum between site sets.
%
% Rows correspond to xyA and columns to xyB.  Measurement nugget is never
% added because independent noise has zero cross-covariance across sites.

    if size(xyA, 2) ~= cfg.model.dimension || ...
            size(xyB, 2) ~= cfg.model.dimension || ...
            any(~isfinite(xyA), 'all') || any(~isfinite(xyB), 'all')
        error('Both coordinate arrays must be finite with the model dimension.');
    end
    distances = eq9_pairwise_distances(xyA, xyB);
    [Hraw, Sraw, kappa] = eq9_raw_half_spectrum( ...
        distances, lambda, theta, cfg);
    normalizer = eq9_spectral_normalizer(theta, cfg);
    Fcross = theta.signalVariance * normalizer * Hraw;
    details.normalizer = normalizer;
    details.rawTemporalPower = Sraw;
    details.kappaByAlias = kappa;
    details.distances = distances;
end

