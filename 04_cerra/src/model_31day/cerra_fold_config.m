function [foldCfg, varianceReference] = cerra_fold_config(cfg, Yobserved, xyObservedKm, optimizer)
%CERRA_FOLD_CONFIG Apply a training-only, scale-adaptive nugget floor.

    [varianceReference, ~, ~] = cerra_empirical_residual_variance( ...
        Yobserved, xyObservedKm);
    foldCfg = cfg;
    floorValue = max(cfg.bounds.noiseVariance(1), ...
        cfg.hsgp.nuggetResidualVarianceFraction * varianceReference);
    foldCfg.bounds.noiseVariance(1) = floorValue;
    foldCfg.optimizer = optimizer;
    foldCfg.hsgp.optimizer = optimizer;
end

