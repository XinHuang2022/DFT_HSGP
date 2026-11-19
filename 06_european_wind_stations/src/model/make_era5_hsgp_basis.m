function basis = make_era5_hsgp_basis(domainXy, boundaryFactor, ...
        basisPerDimension, cfg)
%MAKE_ERA5_HSGP_BASIS Fixed Cartesian Dirichlet Laplacian basis.

    d = cfg.model.dimension;
    if d ~= 2 || size(domainXy, 2) ~= 2 || any(~isfinite(domainXy), 'all')
        error('ERA5 HSGP requires finite two-dimensional coordinates.');
    end
    if ~(isscalar(boundaryFactor) && isfinite(boundaryFactor) && ...
            boundaryFactor > 1)
        error('boundaryFactor must be finite and larger than one.');
    end
    J = basisPerDimension;
    if ~(isscalar(J) && J >= 1 && fix(J) == J)
        error('basisPerDimension must be a positive integer.');
    end

    lower = min(domainXy, [], 1);
    upper = max(domainXy, [], 1);
    center = (lower + upper) / 2;
    supportHalfWidth = (upper - lower) / 2;
    if any(supportHalfWidth <= 0)
        error('The HSGP domain must span both spatial coordinates.');
    end
    L = boundaryFactor * supportHalfWidth;
    [j1, j2] = ndgrid(1:J, 1:J);
    indices = [j1(:), j2(:)];
    omega = pi * indices ./ (2 * L);

    basis.method = 'HSGP';
    basis.boundaryCondition = 'dirichlet';
    basis.dimension = d;
    basis.boundaryFactor = boundaryFactor;
    basis.basisPerDimension = J;
    basis.nFeatures = size(indices, 1);
    basis.center = center;
    basis.supportHalfWidth = supportHalfWidth;
    basis.domainHalfWidth = L;
    basis.indices = indices;
    basis.omegaByDimension = omega;
    basis.laplaceEigenvalue = sum(omega .^ 2, 2);
    basis.radialFrequency = sqrt(basis.laplaceEigenvalue);
    basis.Phi = era5_hsgp_features(domainXy, basis);
end

