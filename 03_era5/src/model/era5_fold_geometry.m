function [observed, target] = era5_fold_geometry(subset, observedIndex, targetIndex, cfg)
%ERA5_FOLD_GEOMETRY Mean and covariance coordinates for a held-out fold.

    xyObservedKm = subset.xyKm(observedIndex, :);
    xyTargetKm = subset.xyKm(targetIndex, :);
    meanCenter = mean(xyObservedKm, 1);
    meanScale = std(xyObservedKm, 0, 1);
    meanScale(meanScale <= sqrt(eps)) = 1;
    observed.xyMean = (xyObservedKm - meanCenter) ./ meanScale;
    target.xyMean = (xyTargetKm - meanCenter) ./ meanScale;
    covarianceCenter = mean(subset.xyKm, 1);
    observed.xyCovariance = (xyObservedKm - covarianceCenter) ./ ...
        cfg.spatial.coordinateScaleKm;
    target.xyCovariance = (xyTargetKm - covarianceCenter) ./ ...
        cfg.spatial.coordinateScaleKm;
    observed.xyKm = xyObservedKm;
    target.xyKm = xyTargetKm;
    T = size(subset.U, 1);
    observed.design = make_era5_mean_design(T, observed.xyMean);
    target.design = make_era5_mean_design(T, target.xyMean);
end
