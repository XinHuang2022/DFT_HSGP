function basis = make_eq9_hsgp_basis(domainXy, boundaryFactor, ...
        basisPerDimension, cfg)
%MAKE_EQ9_HSGP_BASIS Fixed Cartesian Dirichlet Laplacian eigenbasis.

    d = cfg.model.dimension;
    if d ~= 2 || size(domainXy, 2) ~= d || any(~isfinite(domainXy), 'all')
        error('The v0.3.6 HSGP implementation requires finite 2-D sites.');
    end
    if ~(isfinite(boundaryFactor) && boundaryFactor > 1)
        error('HSGP boundaryFactor must be finite and larger than one.');
    end
    J = basisPerDimension;
    if ~(isscalar(J) && J >= 1 && fix(J) == J)
        error('basisPerDimension must be a positive integer.');
    end

    minimum = min(domainXy, [], 1);
    maximum = max(domainXy, [], 1);
    center = (minimum + maximum) / 2;
    supportHalfWidth = (maximum - minimum) / 2;
    if any(supportHalfWidth <= 0)
        error('HSGP sites must span a nonzero interval in every dimension.');
    end
    L = boundaryFactor * supportHalfWidth;
    [j1, j2] = ndgrid(1:J, 1:J);
    indices = [j1(:), j2(:)];
    omega = pi * indices ./ (2 * L);
    laplaceEigenvalue = sum(omega .^ 2, 2);

    basis.method = 'HSGP';
    basis.boundaryCondition = 'dirichlet';
    basis.ordering = 'tensor_square';
    basis.dimension = d;
    basis.boundaryFactor = boundaryFactor;
    basis.basisPerDimension = J;
    basis.nFeatures = size(indices, 1);
    basis.center = center;
    basis.supportHalfWidth = supportHalfWidth;
    basis.domainHalfWidth = L;
    basis.indices = indices;
    basis.omegaByDimension = omega;
    basis.laplaceEigenvalue = laplaceEigenvalue;
    basis.radialFrequency = sqrt(laplaceEigenvalue);
    basis.Phi = eq9_hsgp_features(domainXy, basis);
end
