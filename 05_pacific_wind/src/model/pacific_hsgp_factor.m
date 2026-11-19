function factor = pacific_hsgp_factor(Phi, weight, noiseVariance, cfg)
%PACIFIC_HSGP_FACTOR Automatic stable HSGP covariance factorization.

    if ~(isfinite(noiseVariance) && noiseVariance > 0)
        error('HSGP factorization requires positive noise variance.');
    end
    nSites = size(Phi, 1);
    nFeatures = size(Phi, 2);
    requested = lower(string(cfg.hsgp.factorization));
    useReduced = requested == "feature_qr" || ...
        (requested == "auto" && nFeatures < nSites);
    if useReduced
        core = pacific_rff_factor(Phi, weight, noiseVariance);
        factor.mode = "feature_qr";
        factor.core = core;
        factor.logDeterminant = core.logDeterminant;
        factor.reciprocalCondition = core.reciprocalCondition;
        factor.jitter = core.jitter;
        factor.effectiveDimension = nFeatures;
    else
        K = pacific_hsgp_covariance(Phi, Phi, weight, ...
            noiseVariance, true);
        [L, jitter] = eq9_chol_with_jitter(K);
        factor.mode = "site_cholesky";
        factor.L = L;
        factor.logDeterminant = 2 * sum(log(diag(L)));
        factor.reciprocalCondition = rcond(K + jitter * eye(nSites));
        factor.jitter = jitter;
        factor.effectiveDimension = nSites;
    end
    factor.nSites = nSites;
    factor.nFeatures = nFeatures;
    % The exact numerical rank is cached by the caller when needed.  This
    % inexpensive count is the structural deficit implied by M<n.
    factor.rankDeficientDirectionsBeforeNugget = ...
        max(0, nSites - nFeatures);
end
