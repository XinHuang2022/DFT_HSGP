function prediction = predict_pacific_eq9_hsgp(Yobserved, xyObserved, ...
        xyTarget, observedDesign, targetDesign, fit, cfg, basis)
%PREDICT_PACIFIC_EQ9_HSGP HSGP whole-history held-out reconstruction.

    [T, nObserved] = size(Yobserved);
    nTarget = size(xyTarget, 1);
    fittedObservedMean = evaluate_pacific_mean(observedDesign, fit.eta);
    observedResidual = Yobserved - fittedObservedMean;
    Dobserved = fft(observedResidual, [], 1) / sqrt(T);
    lambda = eq9_fourier_frequencies(T);
    PhiObserved = pacific_hsgp_features(xyObserved, basis);
    PhiTarget = pacific_hsgp_features(xyTarget, basis);
    weights = pacific_hsgp_weights(fit.theta, lambda, basis, cfg);
    p = observedDesign.nCoefficients;

    Dprediction = complex(zeros(T, nTarget));
    sensitivityFrequency = complex(zeros(T, nTarget, p));
    conditionalVarianceFrequency = zeros(T, nTarget);
    maximumJitter = 0;
    minimumReciprocalCondition = Inf;
    featureQrCount = 0;
    siteCholeskyCount = 0;
    for k = 1:T
        weight = weights(:, k);
        factor = pacific_hsgp_factor(PhiObserved, weight, ...
            fit.theta.noiseVariance, cfg);
        solveD = pacific_hsgp_solve(factor, Dobserved(k, :).');
        Kto = pacific_hsgp_covariance(PhiTarget, PhiObserved, ...
            weight, 0, false);
        Dprediction(k, :) = (Kto * solveD).';

        Aobserved = reshape(observedDesign.Afrequency(k, :, :), ...
            nObserved, p);
        Atarget = reshape(targetDesign.Afrequency(k, :, :), nTarget, p);
        solveA = pacific_hsgp_solve(factor, Aobserved);
        sensitivity = Atarget - Kto * solveA;
        sensitivityFrequency(k, :, :) = ...
            reshape(sensitivity, 1, nTarget, p);

        solveCross = pacific_hsgp_solve(factor, Kto');
        Ktt = pacific_hsgp_covariance(PhiTarget, PhiTarget, ...
            weight, 0, true);
        conditional = Ktt - Kto * solveCross;
        conditional = real((conditional + conditional') / 2);
        conditionalVarianceFrequency(k, :) = ...
            max(0, real(diag(conditional))).';
        maximumJitter = max(maximumJitter, factor.jitter);
        minimumReciprocalCondition = min(minimumReciprocalCondition, ...
            factor.reciprocalCondition);
        if factor.mode == "feature_qr"
            featureQrCount = featureQrCount + 1;
        else
            siteCholeskyCount = siteCholeskyCount + 1;
        end
    end

    residualComplex = ifft(Dprediction, [], 1) * sqrt(T);
    residualMean = real(residualComplex);
    modeledMean = evaluate_pacific_mean(targetDesign, fit.eta);
    totalMean = modeledMean + residualMean;
    conditionalVariance = repmat(mean( ...
        conditionalVarianceFrequency, 1), T, 1);
    sensitivityComplex = ifft(sensitivityFrequency, [], 1) * sqrt(T);
    sensitivityTime = real(sensitivityComplex);
    covarianceEta = real((fit.covarianceEta + fit.covarianceEta') / 2);
    fixedEffectVariance = zeros(T, nTarget);
    for j = 1:nTarget
        C = reshape(sensitivityTime(:, j, :), T, p);
        fixedEffectVariance(:, j) = max(0, real(sum( ...
            (C * covarianceEta) .* C, 2)));
    end
    latentVariance = max(cfg.hsgp.minimumPredictiveVariance, ...
        conditionalVariance + fixedEffectVariance);
    observationVariance = max(cfg.hsgp.minimumPredictiveVariance, ...
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
    prediction.maximumImaginaryPrediction = ...
        max(abs(imag(residualComplex)), [], 'all');
    prediction.maximumImaginarySensitivity = ...
        max(abs(imag(sensitivityComplex)), [], 'all');
    prediction.featureQrFrequencyCount = featureQrCount;
    prediction.siteCholeskyFrequencyCount = siteCholeskyCount;
end
