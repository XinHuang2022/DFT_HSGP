function [nll, details] = cerra_profiled_whittle_theta_hsgp( ...
        theta, prepared, cfg, basis)
%CERRA_PROFILED_WHITTLE_THETA_HSGP Profiled HSGP Whittle REML.

    lambdaAll = eq9_fourier_frequencies(prepared.T);
    lambda = lambdaAll(prepared.usedIndices);
    if isfield(prepared, 'hsgpPhi')
        Phi = prepared.hsgpPhi;
    else
        Phi = cerra_hsgp_features(prepared.xy, basis);
    end
    if isfield(prepared,'hsgpFeatureGram')
        featureGram=prepared.hsgpFeatureGram;
    else
        featureGram=Phi'*Phi;
    end
    [weights, weightDiagnostics] = cerra_hsgp_weights( ...
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
        factor = cerra_hsgp_factor(Phi, weights(:, q), ...
            theta.noiseVariance, cfg,featureGram);
        d = real(prepared.D(originalIndex, :).');
        A = real(reshape(prepared.Afrequency( ...
            originalIndex, :, :), m, p));
        solveD = cerra_hsgp_solve(factor, d);
        solveA = cerra_hsgp_solve(factor, A);
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
        factor = cerra_hsgp_factor(Phi, weights(:, offset + q), ...
            theta.noiseVariance, cfg,featureGram);
        d = prepared.D(originalIndex, :).';
        A = reshape(prepared.Afrequency(originalIndex, :, :), m, p);
        solveD = cerra_hsgp_solve(factor, d);
        solveA = cerra_hsgp_solve(factor, A);
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
        error('CERRA HSGP profiled Whittle objective is nonfinite.');
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
    % The deterministic basis is constructed with distinct Laplacian modes.
    % Avoid an SVD of the tall feature matrix at every objective evaluation;
    % the preflight checks the cached Gram matrix separately.
    details.hsgpFeatureRank = min(size(Phi));
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


