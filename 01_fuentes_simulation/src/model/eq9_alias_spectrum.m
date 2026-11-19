function [spectralPower, kappa, tau] = eq9_alias_spectrum(theta, lambda, cfg)
%EQ9_ALIAS_SPECTRUM Alias-specific Equation-9 Matérn components.
%
% At each sampled temporal frequency, the spatial half-spectrum is a
% positive mixture of normalized Matérn spectral densities. spectralPower
% contains the unnormalized mixture masses used by eq9_raw_half_spectrum.

    lambda = lambda(:).';
    aliases = -cfg.model.aliasTerms:cfg.model.aliasTerms;
    nAlias = numel(aliases);
    nFrequency = numel(lambda);
    spectralPower = zeros(nAlias, nFrequency);
    kappa = zeros(nAlias, nFrequency);
    tau = zeros(nAlias, nFrequency);
    d = cfg.model.dimension;

    for a = 1:nAlias
        tau(a, :) = (lambda + 2 * pi * aliases(a)) / ...
            cfg.model.deltaTime;
        A = theta.alpha^2 .* (theta.beta^2 + tau(a, :) .^ 2);
        B = theta.beta^2 + theta.epsilon .* tau(a, :) .^ 2;
        if any(A <= 0) || any(B <= 0)
            error('Invalid Equation-9 alias scale.');
        end
        kappa(a, :) = sqrt(A ./ B);
        spectralPower(a, :) = A .^ (d / 2 - theta.nu) .* ...
            B .^ (-d / 2) / cfg.model.deltaTime;
    end
end
