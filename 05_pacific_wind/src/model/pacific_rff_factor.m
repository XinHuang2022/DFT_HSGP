function factor = pacific_rff_factor(Phi, weight, noiseVariance, ...
        featureGram)
%PACIFIC_RFF_FACTOR Woodbury factor for Phi W Phi* + sigma_n^2 I.

    if ~(isfinite(noiseVariance) && noiseVariance > 0)
        error('RFF Woodbury factor requires a positive nugget variance.');
    end
    weight = weight(:);
    if size(Phi, 2) ~= numel(weight) || any(weight <= 0)
        error('RFF feature and weight dimensions are inconsistent.');
    end
    if nargin < 4 || isempty(featureGram)
        featureGram = Phi' * Phi;
    end
    if ~isequal(size(featureGram), [numel(weight), numel(weight)])
        error('RFF feature Gram matrix has the wrong dimensions.');
    end
    squareRootWeight = sqrt(weight);
    Z = Phi .* squareRootWeight.';
    weightedGram = (squareRootWeight .* featureGram) .* ...
        squareRootWeight.';
    G = eye(size(Z, 2)) + weightedGram / noiseVariance;
    G = (G + G') / 2;
    [L, jitter] = complex_chol(G);
    factor.Z = Z;
    factor.L = L;
    factor.noiseVariance = noiseVariance;
    factor.logDeterminant = size(Phi, 1) * log(noiseVariance) + ...
        2 * sum(log(real(diag(L))));
    factor.reciprocalCondition = rcond(G + jitter * eye(size(G)));
    factor.jitter = jitter;
end

function [L, jitterUsed] = complex_chol(A)
    scale = max(1, mean(abs(diag(A))));
    candidates = scale * [0, 1e-12, 1e-10, 1e-8, 1e-6];
    for jitterUsed = candidates
        [L, flag] = chol(A + jitterUsed * eye(size(A)), 'lower');
        if flag == 0
            return;
        end
    end
    error('RFF feature-space Cholesky factorization failed.');
end
