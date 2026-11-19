function fit = era5_fixed_theta_fit(method, theta, Y, observed, cfg, basis)
%ERA5_FIXED_THETA_FIT Profile only the mean at a supplied covariance theta.

    prepared = prepare_era5_whittle_data(Y, observed.xyCovariance, ...
        observed.design);
    if lower(string(method)) == "exact"
        [objective, details] = era5_profiled_whittle_theta(theta, prepared, cfg);
    else
        prepared.hsgpPhi = era5_hsgp_features(observed.xyCovariance, basis);
        [objective, details] = era5_profiled_whittle_theta_hsgp( ...
            theta, prepared, cfg, basis);
    end
    fit.theta = theta;
    fit.eta = details.eta;
    fit.covarianceEta = details.covarianceEta;
    fit.objective = objective;
    fit.profileNll = details.profileNll;
    fit.remlAdjustment = details.remlAdjustment;
    fit.details = details;
    fit.derived = eq9_derived_parameters(theta, cfg);
end
