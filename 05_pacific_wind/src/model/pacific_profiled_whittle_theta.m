function [nll, details] = pacific_profiled_whittle_theta( ...
        theta, prepared, cfg)
%PACIFIC_PROFILED_WHITTLE_THETA Evaluate exact dense REML at fixed theta.

% This is the common exact objective used both by optimization and by the
% subset plug-in stage.  Calling it directly profiles the fixed effects on
% all outer-training sites without reoptimizing covariance parameters.

    lambdaAll = eq9_fourier_frequencies(prepared.T);
    lambda = lambdaAll(prepared.usedIndices);
    F = eq9_frequency_covariances(theta, prepared.xy, ...
        lambda, cfg, true);
    p = prepared.nCoefficients;
    m = prepared.nStations;
    H = zeros(p, p);
    g = zeros(p, 1);
    baseNll = 0;
    maximumJitter = 0;

    for q = 1:numel(prepared.endpointIndices)
        originalIndex = prepared.endpointIndices(q);
        localIndex = q;
        [L, jitter] = eq9_chol_with_jitter(F(:, :, localIndex));
        d = real(prepared.D(originalIndex, :).');
        A = real(reshape(prepared.Afrequency( ...
            originalIndex, :, :), m, p));
        whitenedD = L \ d;
        whitenedA = L \ A;
        H = H + real(whitenedA' * whitenedA);
        g = g + real(whitenedA' * whitenedD);
        baseNll = baseNll + 0.5 * (2 * sum(log(diag(L))) + ...
            real(whitenedD' * whitenedD) + m * log(2 * pi));
        maximumJitter = max(maximumJitter, jitter);
    end

    offset = numel(prepared.endpointIndices);
    for q = 1:numel(prepared.positiveIndices)
        originalIndex = prepared.positiveIndices(q);
        localIndex = offset + q;
        [L, jitter] = eq9_chol_with_jitter(F(:, :, localIndex));
        d = prepared.D(originalIndex, :).';
        A = reshape(prepared.Afrequency(originalIndex, :, :), m, p);
        whitenedD = L \ d;
        whitenedA = L \ A;
        H = H + 2 * real(whitenedA' * whitenedA);
        g = g + 2 * real(whitenedA' * whitenedD);
        baseNll = baseNll + 2 * sum(log(diag(L))) + ...
            real(whitenedD' * whitenedD) + m * log(pi);
        maximumJitter = max(maximumJitter, jitter);
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
        error('Pacific profiled Whittle objective is nonfinite.');
    end
    details.theta = theta;
    details.eta = eta;
    details.H = H;
    details.g = g;
    details.covarianceEta = LH' \ (LH \ eye(p));
    details.profileNll = profileNll;
    details.remlAdjustment = remlAdjustment;
    details.maximumSpectralJitter = maximumJitter;
    details.meanInformationJitter = meanJitter;
    details.nEndpointFrequencies = numel(prepared.endpointIndices);
    details.nComplexFrequencies = numel(prepared.positiveIndices);
end
