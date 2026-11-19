function prediction = predict_era5_eq9(Yobserved, xyObserved, ...
        xyTarget, observedDesign, targetDesign, fit, cfg)
%PREDICT_ERA5_EQ9 Conditional frequency-domain unseen-site prediction.

    [T, nObserved] = size(Yobserved);
    nTarget = size(xyTarget, 1);
    fittedObservedMean = evaluate_era5_mean(observedDesign, fit.eta);
    observedResidual = Yobserved - fittedObservedMean;
    Dobserved = fft(observedResidual, [], 1) / sqrt(T);
    lambda = eq9_fourier_frequencies(T);
    Foo = eq9_frequency_covariances(fit.theta, xyObserved, ...
        lambda, cfg, true);
    Fto = eq9_cross_frequency_covariances(fit.theta, xyTarget, ...
        xyObserved, lambda, cfg);
    Ftt = eq9_frequency_covariances(fit.theta, xyTarget, ...
        lambda, cfg, false);

    p = observedDesign.nCoefficients;
    Dprediction = complex(zeros(T, nTarget));
    sensitivityFrequency = complex(zeros(T, nTarget, p));
    conditionalVarianceFrequency = zeros(T, nTarget);
    maximumJitter = 0;
    for k = 1:T
        [L, jitter] = eq9_chol_with_jitter(Foo(:, :, k));
        maximumJitter = max(maximumJitter, jitter);
        solveD = L' \ (L \ Dobserved(k, :).');
        Dprediction(k, :) = (Fto(:, :, k) * solveD).';
        Aobserved = reshape(observedDesign.Afrequency(k, :, :), ...
            nObserved, p);
        Atarget = reshape(targetDesign.Afrequency(k, :, :), nTarget, p);
        solveA = L' \ (L \ Aobserved);
        sensitivity = Atarget - Fto(:, :, k) * solveA;
        sensitivityFrequency(k, :, :) = ...
            reshape(sensitivity, 1, nTarget, p);
        solveCross = L' \ (L \ Fto(:, :, k)');
        conditionalCovariance = Ftt(:, :, k) - ...
            Fto(:, :, k) * solveCross;
        conditionalCovariance = real((conditionalCovariance + ...
            conditionalCovariance') / 2);
        conditionalVarianceFrequency(k, :) = ...
            max(0, real(diag(conditionalCovariance))).';
    end

    residualComplex = ifft(Dprediction, [], 1) * sqrt(T);
    residualMean = real(residualComplex);
    modeledMean = evaluate_era5_mean(targetDesign, fit.eta);
    totalMean = modeledMean + residualMean;
    conditionalVariance = repmat(mean( ...
        conditionalVarianceFrequency, 1), T, 1);

    sensitivityTimeComplex = ifft(sensitivityFrequency, [], 1) * sqrt(T);
    sensitivityTime = real(sensitivityTimeComplex);
    covarianceEta = real((fit.covarianceEta + fit.covarianceEta') / 2);
    fixedEffectVariance = zeros(T, nTarget);
    for j = 1:nTarget
        C = reshape(sensitivityTime(:, j, :), T, p);
        fixedEffectVariance(:, j) = max(0, real(sum( ...
            (C * covarianceEta) .* C, 2)));
    end
    latentVariance = max(cfg.cv.minimumPredictiveVariance, ...
        conditionalVariance + fixedEffectVariance);
    observationVariance = max(cfg.cv.minimumPredictiveVariance, ...
        latentVariance + fit.theta.noiseVariance);

    prediction.mean = totalMean;
    prediction.meanOnly = modeledMean;
    prediction.residualMean = residualMean;
    prediction.latentVariance = latentVariance;
    prediction.observationVariance = observationVariance;
    prediction.fixedEffectVariance = fixedEffectVariance;
    prediction.conditionalProcessVariance = conditionalVariance;
    prediction.maximumSpectralJitter = maximumJitter;
    prediction.maximumImaginaryPrediction = ...
        max(abs(imag(residualComplex)), [], 'all');
    prediction.maximumImaginarySensitivity = ...
        max(abs(imag(sensitivityTimeComplex)), [], 'all');
end

