function [weights, diagnostics] = eq9_hsgp_weights(theta, lambda, basis, cfg)
%EQ9_HSGP_WEIGHTS Equation-9 mixture spectrum at Laplacian frequencies.
%
% eq9_matern_spectral_density integrates to one under the suite's direct
% inverse-Fourier convention. Solin-Sarkka HSGP uses S=(2*pi)^d*rho.

    lambda = lambda(:).';
    [aliasPower, aliasKappa] = eq9_alias_spectrum(theta, lambda, cfg);
    normalizer = eq9_spectral_normalizer(theta, cfg);
    eta = theta.nu - cfg.model.dimension / 2;
    nFrequency = numel(lambda);
    weights = zeros(basis.nFeatures, nFrequency);
    omegaForDensity = zeros(basis.nFeatures, cfg.model.dimension);
    omegaForDensity(:, 1) = basis.radialFrequency;

    for k = 1:nFrequency
        density = eq9_matern_spectral_density(omegaForDensity, ...
            aliasKappa(:, k).', eta);
        weights(:, k) = theta.signalVariance * normalizer * ...
            (2 * pi) ^ cfg.model.dimension * ...
            (density * aliasPower(:, k));
    end
    if any(~isfinite(weights), 'all') || any(weights <= 0, 'all')
        error('HSGP weights must be positive and finite.');
    end

    diagnostics.aliasPower = aliasPower;
    diagnostics.aliasKappa = aliasKappa;
    diagnostics.normalizer = normalizer;
    diagnostics.minimumWeight = min(weights, [], 'all');
    diagnostics.maximumWeight = max(weights, [], 'all');
    diagnostics.maximumToMedianWeight = max(weights, [], 1) ./ ...
        max(realmin('double'), median(weights, 1));
end
