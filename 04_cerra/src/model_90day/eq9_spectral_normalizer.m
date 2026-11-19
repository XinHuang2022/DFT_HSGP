function normalizer = eq9_spectral_normalizer(theta, cfg)
%EQ9_SPECTRAL_NORMALIZER Scale giving unit latent marginal variance.
%
% With a unitary DFT, C(0,0) is the mean of the discrete spectrum over the
% Fourier frequencies. A dense, T-independent grid is used so changing the
% observed record length does not redefine the covariance parameters.

    nGrid = cfg.model.normalizationGridSize;
    lambda = 2 * pi * (0:nGrid-1) / nGrid;
    lambda(lambda > pi) = lambda(lambda > pi) - 2 * pi;
    [~, Sraw] = eq9_raw_half_spectrum(0, lambda, theta, cfg);
    meanPower = mean(Sraw);
    if ~(isfinite(meanPower) && meanPower > 0)
        error('The Equation-9 spectral normalizer is not finite and positive.');
    end
    normalizer = 1 / meanPower;
end



