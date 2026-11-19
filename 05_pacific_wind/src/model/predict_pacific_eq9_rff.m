function prediction = predict_pacific_eq9_rff(Yobserved, xyObserved, ...
        xyTarget, observedDesign, targetDesign, fit, cfg, bank)
%PREDICT_PACIFIC_EQ9_RFF RFF/Woodbury unseen-location prediction.

    [T, nObserved] = size(Yobserved);
    nTarget = size(xyTarget, 1);
    fittedObservedMean = evaluate_pacific_mean(observedDesign, fit.eta);
    observedResidual = Yobserved - fittedObservedMean;
    Dobserved = fft(observedResidual, [], 1) / sqrt(T);
    lambda = eq9_fourier_frequencies(T);
    PhiObserved = pacific_rff_features(xyObserved, bank);
    PhiTarget = pacific_rff_features(xyTarget, bank);
    observedFeatureGram = PhiObserved' * PhiObserved;
    [weights, weightDiagnostics] = pacific_rff_weights( ...
        fit.theta, lambda, bank, cfg);

    pMean = observedDesign.nCoefficients;
    Dprediction = complex(zeros(T, nTarget));
    sensitivityFrequency = complex(zeros(T, nTarget, pMean));
    conditionalVarianceFrequency = zeros(T, nTarget);
    maximumJitter = 0;
    minimumReciprocalCondition = Inf;
    for k = 1:T
        weight = weights(:, k);
        factor = pacific_rff_factor(PhiObserved, weight, ...
            fit.theta.noiseVariance, observedFeatureGram);
        solveD = pacific_rff_solve(factor, Dobserved(k, :).');
        Dprediction(k, :) = (PhiTarget * ...
            (weight .* (PhiObserved' * solveD))).';
        Aobserved = reshape(observedDesign.Afrequency(k, :, :), ...
            nObserved, pMean);
        Atarget = reshape(targetDesign.Afrequency(k, :, :), ...
            nTarget, pMean);
        solveA = pacific_rff_solve(factor, Aobserved);
        sensitivity = Atarget - PhiTarget * ...
            (weight .* (PhiObserved' * solveA));
        sensitivityFrequency(k, :, :) = ...
            reshape(sensitivity, 1, nTarget, pMean);

        Kot = PhiObserved * (weight .* PhiTarget');
        solveCross = pacific_rff_solve(factor, Kot);
        Ktt = PhiTarget * (weight .* PhiTarget');
        conditionalCovariance = Ktt - Kot' * solveCross;
        conditionalCovariance = real((conditionalCovariance + ...
            conditionalCovariance') / 2);
        conditionalVarianceFrequency(k, :) = ...
            max(0, real(diag(conditionalCovariance))).';
        maximumJitter = max(maximumJitter, factor.jitter);
        minimumReciprocalCondition = min(minimumReciprocalCondition, ...
            factor.reciprocalCondition);
    end

    residualComplex = ifft(Dprediction, [], 1) * sqrt(T);
    residualMean = real(residualComplex);
    modeledMean = evaluate_pacific_mean(targetDesign, fit.eta);
    totalMean = modeledMean + residualMean;
    conditionalVariance = repmat(mean( ...
        conditionalVarianceFrequency, 1), T, 1);
    sensitivityTimeComplex = ifft(sensitivityFrequency, [], 1) * sqrt(T);
    sensitivityTime = real(sensitivityTimeComplex);
    covarianceEta = real((fit.covarianceEta + fit.covarianceEta') / 2);
    fixedEffectVariance = zeros(T, nTarget);
    for j = 1:nTarget
        C = reshape(sensitivityTime(:, j, :), T, pMean);
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
    prediction.maximumFeatureJitter = maximumJitter;
    prediction.minimumFeatureReciprocalCondition = ...
        minimumReciprocalCondition;
    prediction.minimumEssFraction = ...
        weightDiagnostics.minimumEssFraction;
    prediction.medianEssFraction = ...
        weightDiagnostics.medianEssFraction;
    prediction.maximumImaginaryPrediction = ...
        max(abs(imag(residualComplex)), [], 'all');
    prediction.maximumImaginarySensitivity = ...
        max(abs(imag(sensitivityTimeComplex)), [], 'all');
end
