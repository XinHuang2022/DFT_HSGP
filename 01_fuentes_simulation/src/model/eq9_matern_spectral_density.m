function density = eq9_matern_spectral_density(omega, kappa, eta)
%EQ9_MATERN_SPECTRAL_DENSITY Normalized Matérn spatial spectrum.
%
% Rows of omega are d-dimensional frequencies. kappa can be a scalar or a
% row vector. The density integrates to one and its Fourier transform is
% M_eta(kappa*||h||), under the convention used by Equation 9.

    if eta <= 0 || any(kappa <= 0)
        error('Matérn eta and kappa must be positive.');
    end
    d = size(omega, 2);
    radiusSquared = sum(omega .^ 2, 2);
    kappa = kappa(:).';
    exponent = eta + d / 2;
    logConstant = gammaln(exponent) - gammaln(eta) - ...
        (d / 2) * log(pi);
    logDensity = logConstant + 2 * eta * log(kappa) - ...
        exponent * log(radiusSquared + kappa .^ 2);
    density = exp(logDensity);
end
