function density = pacific_matern_spectral_density(omega, kappa, eta)
%PACIFIC_MATERN_SPECTRAL_DENSITY Normalized Matérn frequency density.

% For d-dimensional spatial coordinates and Matérn correlation M_eta,
%
% p_kappa(w) = Gamma(eta+d/2) kappa^(2 eta) /
%              [pi^(d/2) Gamma(eta)
%               (kappa^2 + ||w||^2)^(eta+d/2)],
%
% so integral exp(i w' h) p_kappa(w) dw = M_eta(kappa ||h||).
% Rows of omega are frequencies; kappa may be a scalar or row vector.

    if any(kappa <= 0) || eta <= 0
        error('Matérn kappa and eta must be positive.');
    end
    d = size(omega, 2);
    radiusSquared = sum(omega .^ 2, 2);
    kappa = kappa(:).';
    gamma = eta + d / 2;
    logConstant = gammaln(gamma) - gammaln(eta) - ...
        (d / 2) * log(pi);
    logDensity = logConstant + 2 * eta * log(kappa) - ...
        gamma * log(radiusSquared + kappa .^ 2);
    density = exp(logDensity);
end
