function [nll, details] = pacific_profiled_whittle_theta_rff( ...
        theta, prepared, cfg, bank)
%PACIFIC_PROFILED_WHITTLE_THETA_RFF RFF Woodbury Whittle REML.

% No m-by-m covariance is factorized.  Each temporal frequency uses a
% feature-space factorization of I + W^(1/2) Phi*Phi W^(1/2)/sigma_n^2.

    lambdaAll = eq9_fourier_frequencies(prepared.T);
    lambda = lambdaAll(prepared.usedIndices);
    if isfield(prepared, 'rffPhi') && ...
            isfield(prepared, 'rffFeatureGram')
        Phi = prepared.rffPhi;
        featureGram = prepared.rffFeatureGram;
    else
        Phi = pacific_rff_features(prepared.xy, bank);
        featureGram = Phi' * Phi;
    end
    [weights, weightDiagnostics] = pacific_rff_weights( ...
        theta, lambda, bank, cfg);
    pMean = prepared.nCoefficients;
    m = prepared.nStations;
    H = zeros(pMean, pMean);
    g = zeros(pMean, 1);
    baseNll = 0;
    maximumJitter = 0;
    minimumReciprocalCondition = Inf;

    for q = 1:numel(prepared.endpointIndices)
        originalIndex = prepared.endpointIndices(q);
        localIndex = q;
        factor = pacific_rff_factor(Phi, weights(:, localIndex), ...
            theta.noiseVariance, featureGram);
        d = real(prepared.D(originalIndex, :).');
        A = real(reshape(prepared.Afrequency( ...
            originalIndex, :, :), m, pMean));
        solveD = pacific_rff_solve(factor, d);
        solveA = pacific_rff_solve(factor, A);
        H = H + real(A' * solveA);
        g = g + real(A' * solveD);
        baseNll = baseNll + 0.5 * (factor.logDeterminant + ...
            real(d' * solveD) + m * log(2 * pi));
        maximumJitter = max(maximumJitter, factor.jitter);
        minimumReciprocalCondition = min(minimumReciprocalCondition, ...
            factor.reciprocalCondition);
    end

    offset = numel(prepared.endpointIndices);
    for q = 1:numel(prepared.positiveIndices)
        originalIndex = prepared.positiveIndices(q);
        localIndex = offset + q;
        factor = pacific_rff_factor(Phi, weights(:, localIndex), ...
            theta.noiseVariance, featureGram);
        d = prepared.D(originalIndex, :).';
        A = reshape(prepared.Afrequency(originalIndex, :, :), m, pMean);
        solveD = pacific_rff_solve(factor, d);
        solveA = pacific_rff_solve(factor, A);
        H = H + 2 * real(A' * solveA);
        g = g + 2 * real(A' * solveD);
        baseNll = baseNll + factor.logDeterminant + ...
            real(d' * solveD) + m * log(pi);
        maximumJitter = max(maximumJitter, factor.jitter);
        minimumReciprocalCondition = min(minimumReciprocalCondition, ...
            factor.reciprocalCondition);
    end

    H = real((H + H') / 2);
    [LH, meanJitter] = eq9_chol_with_jitter(H);
    eta = LH' \ (LH \ g);
    profileNll = baseNll - 0.5 * real(eta' * g);
    remlAdjustment = sum(log(diag(LH)));
    if strcmpi(cfg.likelihood.type, 'whittle_reml')
        nll = profileNll + remlAdjustment;
    elseif strcmpi(cfg.likelihood.type, 'profile_whittle')
        nll = profileNll;
    else
        error('Unknown likelihood type %s.', cfg.likelihood.type);
    end
    if ~isfinite(nll)
        error('Pacific RFF Whittle objective is nonfinite.');
    end
    details.theta = theta;
    details.eta = eta;
    details.H = H;
    details.g = g;
    details.covarianceEta = LH' \ (LH \ eye(pMean));
    details.profileNll = profileNll;
    details.remlAdjustment = remlAdjustment;
    details.maximumFeatureJitter = maximumJitter;
    details.minimumFeatureReciprocalCondition = ...
        minimumReciprocalCondition;
    details.meanInformationJitter = meanJitter;
    details.minimumEssFraction = ...
        weightDiagnostics.minimumEssFraction;
    details.medianEssFraction = ...
        weightDiagnostics.medianEssFraction;
end
