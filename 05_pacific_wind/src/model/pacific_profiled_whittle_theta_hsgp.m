function [nll, details] = pacific_profiled_whittle_theta_hsgp( ...
        theta, prepared, cfg, basis)
%PACIFIC_PROFILED_WHITTLE_THETA_HSGP Profiled HSGP Whittle REML.

    lambdaAll = eq9_fourier_frequencies(prepared.T);
    lambda = lambdaAll(prepared.usedIndices);
    if isfield(prepared, 'hsgpPhi')
        Phi = prepared.hsgpPhi;
    else
        Phi = pacific_hsgp_features(prepared.xy, basis);
    end
    [weights, weightDiagnostics] = pacific_hsgp_weights( ...
        theta, lambda, basis, cfg);
    p = prepared.nCoefficients;
    m = prepared.nStations;
    H = zeros(p, p);
    g = zeros(p, 1);
    baseNll = 0;
    maximumJitter = 0;
    minimumReciprocalCondition = Inf;
    featureQrCount = 0;
    siteCholeskyCount = 0;

    for q = 1:numel(prepared.endpointIndices)
        originalIndex = prepared.endpointIndices(q);
        factor = pacific_hsgp_factor(Phi, weights(:, q), ...
            theta.noiseVariance, cfg);
        d = real(prepared.D(originalIndex, :).');
        A = real(reshape(prepared.Afrequency( ...
            originalIndex, :, :), m, p));
        solveD = pacific_hsgp_solve(factor, d);
        solveA = pacific_hsgp_solve(factor, A);
        H = H + real(A' * solveA);
        g = g + real(A' * solveD);
        baseNll = baseNll + 0.5 * (factor.logDeterminant + ...
            real(d' * solveD) + m * log(2 * pi));
        [maximumJitter, minimumReciprocalCondition, ...
            featureQrCount, siteCholeskyCount] = update_diagnostics( ...
            factor, maximumJitter, minimumReciprocalCondition, ...
            featureQrCount, siteCholeskyCount);
    end

    offset = numel(prepared.endpointIndices);
    for q = 1:numel(prepared.positiveIndices)
        originalIndex = prepared.positiveIndices(q);
        factor = pacific_hsgp_factor(Phi, weights(:, offset + q), ...
            theta.noiseVariance, cfg);
        d = prepared.D(originalIndex, :).';
        A = reshape(prepared.Afrequency(originalIndex, :, :), m, p);
        solveD = pacific_hsgp_solve(factor, d);
        solveA = pacific_hsgp_solve(factor, A);
        H = H + 2 * real(A' * solveA);
        g = g + 2 * real(A' * solveD);
        baseNll = baseNll + factor.logDeterminant + ...
            real(d' * solveD) + m * log(pi);
        [maximumJitter, minimumReciprocalCondition, ...
            featureQrCount, siteCholeskyCount] = update_diagnostics( ...
            factor, maximumJitter, minimumReciprocalCondition, ...
            featureQrCount, siteCholeskyCount);
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
        error('Pacific HSGP profiled Whittle objective is nonfinite.');
    end

    details.theta = theta;
    details.eta = eta;
    details.H = H;
    details.g = g;
    details.covarianceEta = LH' \ (LH \ eye(p));
    details.profileNll = profileNll;
    details.remlAdjustment = remlAdjustment;
    details.maximumFeatureJitter = maximumJitter;
    details.minimumFeatureReciprocalCondition = ...
        minimumReciprocalCondition;
    details.meanInformationJitter = meanJitter;
    details.minimumWeight = weightDiagnostics.minimumWeight;
    details.maximumWeight = weightDiagnostics.maximumWeight;
    details.featureQrFrequencyCount = featureQrCount;
    details.siteCholeskyFrequencyCount = siteCholeskyCount;
    details.hsgpFeatureRank = rank(Phi);
    details.rankDeficientDirectionsBeforeNugget = ...
        max(0, m - details.hsgpFeatureRank);
end

function [maximumJitter, minimumRcond, qrCount, cholCount] = ...
        update_diagnostics(factor, maximumJitter, minimumRcond, ...
        qrCount, cholCount)
    maximumJitter = max(maximumJitter, factor.jitter);
    minimumRcond = min(minimumRcond, factor.reciprocalCondition);
    if factor.mode == "feature_qr"
        qrCount = qrCount + 1;
    else
        cholCount = cholCount + 1;
    end
end
