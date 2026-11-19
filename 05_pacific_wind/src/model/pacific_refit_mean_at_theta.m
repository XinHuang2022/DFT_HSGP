function fit = pacific_refit_mean_at_theta(Y, xyCovariance, design, ...
        theta, cfg)
%PACIFIC_REFIT_MEAN_AT_THETA Profile fixed effects at fixed covariance.

% This single dense objective evaluation is the bridge between covariance
% pretraining on a spatial subset and prediction conditional on every
% outer-training site.  No covariance parameter is updated here.

    prepared = prepare_pacific_whittle_data(Y, xyCovariance, design);
    [objective, details] = pacific_profiled_whittle_theta( ...
        theta, prepared, cfg);
    fit.theta = theta;
    fit.eta = details.eta;
    fit.covarianceEta = details.covarianceEta;
    fit.objective = objective;
    fit.profileNll = details.profileNll;
    fit.remlAdjustment = details.remlAdjustment;
    fit.details = details;
    fit.derived = eq9_derived_parameters(theta, cfg);
    fit.derived.spatialRange10LowFrequencyKm = ...
        fit.derived.spatialRange10LowFrequency * ...
        cfg.spatial.coordinateScaleKm;
    fit.derived.spatialRange10HighFrequencyKm = ...
        fit.derived.spatialRange10HighFrequency * ...
        cfg.spatial.coordinateScaleKm;
    fit.derived.timeScaleHours = cfg.data.samplingHours / theta.beta;
end
