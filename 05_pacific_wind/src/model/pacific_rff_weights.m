function [weights, diagnostics] = pacific_rff_weights( ...
        theta, lambda, bank, cfg)
%PACIFIC_RFF_WEIGHTS Importance weights for the full Equation-9 alias sum.

% The proposal is fixed by bank.thetaReference, but the target numerator is
% recomputed for every theta.  The denominator is always the complete
% mixture q(omega), not the density of the sampled component.

    lambda = lambda(:).';
    [aliasPower, aliasKappa] = pacific_eq9_alias_spectrum( ...
        theta, lambda, cfg);
    normalizer = eq9_spectral_normalizer(theta, cfg);
    p = bank.nFeatures;
    nFrequency = numel(lambda);
    weights = zeros(p, nFrequency);
    eta = theta.nu - cfg.model.dimension / 2;
    for k = 1:nFrequency
        density = pacific_matern_spectral_density( ...
            bank.omega, aliasKappa(:, k).', eta);
        targetDensity = density * aliasPower(:, k);
        weights(:, k) = theta.signalVariance * normalizer * ...
            targetDensity ./ bank.proposalDensity;
    end
    if any(~isfinite(weights), 'all') || any(weights < 0, 'all')
        error('RFF importance weights are invalid.');
    end
    weights = max(realmin('double'), weights);
    sumWeight = sum(weights, 1);
    effectiveSampleSize = sumWeight .^ 2 ./ sum(weights .^ 2, 1);
    diagnostics.effectiveSampleSize = effectiveSampleSize;
    diagnostics.essFraction = effectiveSampleSize / p;
    diagnostics.minimumEssFraction = min(diagnostics.essFraction);
    diagnostics.medianEssFraction = median(diagnostics.essFraction);
    diagnostics.maximumWeightRatio = max(weights, [], 1) ./ ...
        max(realmin('double'), mean(weights, 1));
    diagnostics.aliasPower = aliasPower;
    diagnostics.aliasKappa = aliasKappa;
    diagnostics.normalizer = normalizer;
end
