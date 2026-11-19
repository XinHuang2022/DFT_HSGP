function [foldCfg, varianceReference] = era5_fold_config(cfg, Yobserved, xyObservedKm, optimizer)
%ERA5_FOLD_CONFIG Apply a training-only, scale-adaptive nugget floor.

    [varianceReference, ~, ~] = era5_empirical_residual_variance( ...
        Yobserved, xyObservedKm);
    foldCfg = cfg;
    floorValue = max(cfg.bounds.noiseVariance(1), ...
        cfg.hsgp.nuggetResidualVarianceFraction * varianceReference);
    foldCfg.bounds.noiseVariance(1) = floorValue;
    foldCfg.optimizer = optimizer;
    foldCfg.hsgp.optimizer = optimizer;
end
