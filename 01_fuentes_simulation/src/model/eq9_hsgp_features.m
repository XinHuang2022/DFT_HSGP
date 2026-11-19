function Phi = eq9_hsgp_features(xy, basis)
%EQ9_HSGP_FEATURES Evaluate normalized Dirichlet eigenfunctions.

    if size(xy, 2) ~= basis.dimension || any(~isfinite(xy), 'all')
        error('Coordinates do not match the HSGP basis dimension.');
    end
    centered = xy - basis.center;
    if any(abs(centered) > basis.domainHalfWidth + 1e-12, 'all')
        error('At least one site lies outside the fixed HSGP domain.');
    end
    Phi = ones(size(xy, 1), basis.nFeatures);
    for k = 1:basis.dimension
        omega = basis.omegaByDimension(:, k).';
        Phi = Phi .* sin((centered(:, k) + ...
            basis.domainHalfWidth(k)) * omega) / ...
            sqrt(basis.domainHalfWidth(k));
    end
end
