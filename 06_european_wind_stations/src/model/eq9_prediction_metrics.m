function metrics = eq9_prediction_metrics(truth, meanPrediction, ...
    predictiveVariance, nominalCoverage)
%EQ9_PREDICTION_METRICS Point accuracy and Gaussian predictive calibration.

    truth = truth(:);
    meanPrediction = meanPrediction(:);
    predictiveVariance = predictiveVariance(:);
    keep = isfinite(truth) & isfinite(meanPrediction) & ...
        isfinite(predictiveVariance) & predictiveVariance > 0;
    truth = truth(keep);
    meanPrediction = meanPrediction(keep);
    predictiveVariance = predictiveVariance(keep);
    if isempty(truth)
        error('No finite values are available for prediction metrics.');
    end
    errorValues = meanPrediction - truth;
    metrics.nValues = numel(truth);
    metrics.rmse = sqrt(mean(errorValues.^2));
    metrics.mae = mean(abs(errorValues));
    centeredTruth = truth - mean(truth);
    centeredPrediction = meanPrediction - mean(meanPrediction);
    centeredDenominator = sum(centeredTruth.^2);
    metrics.centeredUnexplainedEnergy = sum( ...
        (centeredPrediction - centeredTruth).^2) / ...
        max(eps, centeredDenominator);
    correlationMatrix = corrcoef(truth, meanPrediction);
    if numel(correlationMatrix) == 4
        metrics.correlation = correlationMatrix(1, 2);
    else
        metrics.correlation = NaN;
    end
    standardDeviation = sqrt(predictiveVariance);
    standardized = (truth - meanPrediction) ./ standardDeviation;
    metrics.nlpd = mean(0.5 * log(2 * pi * predictiveVariance) + ...
        0.5 * standardized.^2);
    normalPdf = exp(-0.5 * standardized.^2) / sqrt(2 * pi);
    normalCdf = 0.5 * (1 + erf(standardized / sqrt(2)));
    crps = standardDeviation .* (standardized .* ...
        (2 * normalCdf - 1) + 2 * normalPdf - 1 / sqrt(pi));
    metrics.crps = mean(crps);
    criticalValue = sqrt(2) * erfinv(nominalCoverage);
    lower = meanPrediction - criticalValue * standardDeviation;
    upper = meanPrediction + criticalValue * standardDeviation;
    metrics.coverage = mean(truth >= lower & truth <= upper);
    metrics.meanIntervalWidth = mean(upper - lower);
    metrics.meanPredictiveSd = mean(standardDeviation);
end
