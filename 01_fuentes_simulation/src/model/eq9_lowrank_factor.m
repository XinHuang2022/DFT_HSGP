function factor = eq9_lowrank_factor(Phi, weight, noiseVariance)
%EQ9_LOWRANK_FACTOR Stable QR factorization of noise*I+Phi*W*Phi'.

    if ~(isfinite(noiseVariance) && noiseVariance > 0)
        error('Low-rank factorization requires positive noise variance.');
    end
    weight = weight(:);
    if size(Phi, 2) ~= numel(weight) || any(weight <= 0)
        error('Low-rank feature and weight dimensions are inconsistent.');
    end
    Z = Phi .* sqrt(weight).';
    [Q, R] = qr(Z, 0);
    reducedDimension = size(Q, 2);
    small = noiseVariance * eye(reducedDimension) + R * R';
    small = real((small + small') / 2);
    [L, jitter] = chol_with_jitter(small);
    nSites = size(Phi, 1);

    factor.Q = Q;
    factor.R = R;
    factor.L = L;
    factor.noiseVariance = noiseVariance;
    factor.reducedDimension = reducedDimension;
    factor.logDeterminant = (nSites - reducedDimension) * ...
        log(noiseVariance) + 2 * sum(log(diag(L)));
    factor.reciprocalCondition = rcond( ...
        small + jitter * eye(reducedDimension));
    factor.jitter = jitter;
end

function [L, jitterUsed] = chol_with_jitter(A)
    scale = max(1, mean(abs(diag(A))));
    candidates = scale * [0, 1e-13, 1e-11, 1e-9, 1e-7];
    for jitterUsed = candidates
        [L, flag] = chol(A + jitterUsed * eye(size(A)), 'lower');
        if flag == 0
            return;
        end
    end
    error('Reduced HSGP Cholesky factorization failed.');
end
