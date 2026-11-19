function derived = eq9_derived_parameters(theta, cfg)
%EQ9_DERIVED_PARAMETERS Interpretable Equation-9 scale summaries.
%
% kappaLow  = alpha is the spatial inverse scale at temporal frequency 0.
% kappaHigh = alpha/sqrt(epsilon) is its asymptotic high-frequency value.
% delta      = log(kappaHigh/kappaLow) = -0.5*log(epsilon).

    eta = theta.nu - cfg.model.dimension / 2;
    targetCorrelation = 0.10;
    rootFunction = @(x) matern_correlation(x, eta) - targetCorrelation;
    maternDistance10 = fzero(rootFunction, [1e-8, 50]);

    derived.logSpatialScaleContrast = -0.5 * log(theta.epsilon);
    derived.spatialRangeRatioHighToLow = sqrt(theta.epsilon);
    derived.kappaLowFrequency = theta.alpha;
    derived.kappaHighFrequency = theta.alpha / sqrt(theta.epsilon);
    derived.spatialRange10LowFrequency = maternDistance10 / ...
        derived.kappaLowFrequency;
    derived.spatialRange10HighFrequency = maternDistance10 / ...
        derived.kappaHighFrequency;
    derived.inverseBetaTimeScale = 1 / theta.beta;
end



