function [Hraw, Sraw, kappa] = eq9_raw_half_spectrum(distances, lambda, theta, cfg)
%EQ9_RAW_HALF_SPECTRUM Unnormalized discrete-time Equation-9 half-spectrum.
%
% The continuous half-spectrum is summed over temporal aliases. Constants
% independent of frequency and spatial lag are omitted because a separate
% normalizer enforces C(0,0)=signalVariance.

    lambda = lambda(:)';
    nFrequency = numel(lambda);
    nA = size(distances, 1);
    nB = size(distances, 2);
    Hraw = zeros(nA, nB, nFrequency);
    Sraw = zeros(1, nFrequency);
    kappa = zeros(2 * cfg.model.aliasTerms + 1, nFrequency);

    d = cfg.model.dimension;
    eta = theta.nu - d / 2;
    if eta <= 0
        error('Equation 9 requires nu > d/2 for the spatial integral.');
    end

    aliasIndex = 0;
    for q = -cfg.model.aliasTerms:cfg.model.aliasTerms
        aliasIndex = aliasIndex + 1;
        tau = (lambda + 2 * pi * q) / cfg.model.deltaTime;
        A = theta.alpha^2 .* (theta.beta^2 + tau.^2);
        B = theta.beta^2 + theta.epsilon .* tau.^2;
        if any(A <= 0) || any(B <= 0)
            error('Invalid Equation-9 spectral scale.');
        end
        spatialKappa = sqrt(A ./ B);
        spectralPower = A.^(d / 2 - theta.nu) .* B.^(-d / 2) / ...
            cfg.model.deltaTime;
        Sraw = Sraw + spectralPower;
        kappa(aliasIndex, :) = spatialKappa;

        scaledDistance = distances .* reshape(spatialKappa, 1, 1, []);
        spatialCorrelation = matern_correlation(scaledDistance, eta);
        Hraw = Hraw + spatialCorrelation .* ...
            reshape(spectralPower, 1, 1, []);
    end
end
