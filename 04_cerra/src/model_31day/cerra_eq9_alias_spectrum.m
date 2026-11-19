function [spectralPower, kappa, tau] = ...
        cerra_eq9_alias_spectrum(theta, lambda, cfg)
%CERRA_EQ9_ALIAS_SPECTRUM Alias-specific Equation-9 spatial spectra.

% For each temporal DFT frequency lambda(k), Equation 9 is a positive
% mixture of normalized Matérn spatial spectra.  This function exposes the
% mixture weights and inverse spatial scales used by the RFF integral.

    lambda = lambda(:).';
    aliasValues = -cfg.model.aliasTerms:cfg.model.aliasTerms;
    nAlias = numel(aliasValues);
    spectralPower = zeros(nAlias, numel(lambda));
    kappa = zeros(nAlias, numel(lambda));
    tau = zeros(nAlias, numel(lambda));
    d = cfg.model.dimension;
    for a = 1:nAlias
        tau(a, :) = (lambda + 2 * pi * aliasValues(a)) / ...
            cfg.model.deltaTime;
        A = theta.alpha^2 .* (theta.beta^2 + tau(a, :) .^ 2);
        B = theta.beta^2 + theta.epsilon .* tau(a, :) .^ 2;
        kappa(a, :) = sqrt(A ./ B);
        spectralPower(a, :) = A .^ (d / 2 - theta.nu) .* ...
            B .^ (-d / 2) / cfg.model.deltaTime;
    end
end


