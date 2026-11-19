function [weights, diagnostics] = cerra_hsgp_weights( ...
        theta, lambda, basis, cfg)
%CERRA_HSGP_WEIGHTS Equation-9 spectrum at Laplacian frequencies.
%
% The CERRA Matérn density integrates to one under the direct Fourier
% convention.  The HSGP operator weight is therefore (2*pi)^d times that
% density, summed over the temporal aliases and multiplied by the Equation
% 9 normalizer and signal variance.

    lambda = lambda(:).';
    [aliasPower, aliasKappa] = cerra_eq9_alias_spectrum( ...
        theta, lambda, cfg);
    normalizer = eq9_spectral_normalizer(theta, cfg);
    eta = theta.nu - cfg.model.dimension / 2;
    nFrequency = numel(lambda);
    weights = zeros(basis.nFeatures, nFrequency);
    omegaForDensity = zeros(basis.nFeatures, cfg.model.dimension);
    omegaForDensity(:, 1) = basis.radialFrequency;
    for k = 1:nFrequency
        density = cerra_matern_spectral_density(omegaForDensity, ...
            aliasKappa(:, k).', eta);
        weights(:, k) = theta.signalVariance * normalizer * ...
            (2 * pi) ^ cfg.model.dimension * ...
            (density * aliasPower(:, k));
    end
    if any(~isfinite(weights), 'all') || any(weights <= 0, 'all')
        error('CERRA HSGP weights must be positive and finite.');
    end
    diagnostics.minimumWeight = min(weights, [], 'all');
    diagnostics.maximumWeight = max(weights, [], 'all');
    diagnostics.aliasPower = aliasPower;
    diagnostics.aliasKappa = aliasKappa;
    diagnostics.normalizer = normalizer;
end


