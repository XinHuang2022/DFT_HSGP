function Phi = pacific_rff_features(xy, bank)
%PACIFIC_RFF_FEATURES Complex paired Fourier features, scaled by p^(-1/2).

    if size(xy, 2) ~= bank.dimension
        error('Coordinate and RFF dimensions do not agree.');
    end
    Phi = exp(1i * (xy * bank.omega.')) / sqrt(bank.nFeatures);
end
