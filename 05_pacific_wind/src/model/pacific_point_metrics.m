function metrics = pacific_point_metrics(truth, prediction)
%PACIFIC_POINT_METRICS Point-prediction metrics without uncertainty.

    truth = truth(:);
    prediction = prediction(:);
    keep = isfinite(truth) & isfinite(prediction);
    truth = truth(keep);
    prediction = prediction(keep);
    errorValue = prediction - truth;
    metrics.nValues = numel(truth);
    metrics.rmse = sqrt(mean(errorValue .^ 2));
    metrics.mae = mean(abs(errorValue));
    denominator = sum((truth - mean(truth)) .^ 2);
    metrics.centeredUnexplainedEnergy = sum(((prediction - ...
        mean(prediction)) - (truth - mean(truth))) .^ 2) / ...
        max(eps, denominator);
    C = corrcoef(truth, prediction);
    if numel(C) == 4
        metrics.correlation = C(1, 2);
    else
        metrics.correlation = NaN;
    end
end
