function fit = cerra_fixed_theta_fit(method, theta, Y, observed, cfg, basis)
%CERRA_FIXED_THETA_FIT Profile only the mean at a supplied covariance theta.

    prepared = prepare_cerra_whittle_data(Y, observed.xyCovariance, ...
        observed.design);
    if lower(string(method)) == "exact"
        [objective, details] = cerra_profiled_whittle_theta(theta, prepared, cfg);
    else
        prepared.hsgpPhi = cerra_hsgp_features(observed.xyCovariance, basis);
        prepared.hsgpFeatureGram = prepared.hsgpPhi' * prepared.hsgpPhi;
        [objective, details] = cerra_profiled_whittle_theta_hsgp( ...
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

