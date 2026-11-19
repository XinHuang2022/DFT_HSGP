function report=run_era5_preflight()
%RUN_ERA5_PREFLIGHT Fast data, geometry, exact, and HSGP interface checks.
    cfg=era5_northsea_config();cfg.model.normalizationGridSize=256;
    data=load_era5_wind_data(cfg);subset=make_era5_subset(data,'central5');
    [folds,manifest,summary]=make_era5_cv_folds(subset,'balanced5');
    targetIndex=folds(1).targetIndex;observedIndex=setdiff((1:25)',targetIndex);
    [observed,target]=era5_fold_geometry(subset,observedIndex,targetIndex,cfg);
    [foldCfg,~]=era5_fold_config(cfg,subset.U(:,observedIndex), ...
        observed.xyKm,cfg.stage2.optimizer);
    exact=era5_fixed_theta_fit("exact",cfg.start,subset.U(:,observedIndex), ...
        observed,foldCfg,[]);
    domain=(subset.xyKm-mean(subset.xyKm,1))/cfg.spatial.coordinateScaleKm;
    basis=make_era5_hsgp_basis(domain,4,6,cfg);
    hsgp=era5_fixed_theta_fit("hsgp",cfg.start,subset.U(:,observedIndex), ...
        observed,foldCfg,basis);
    ep=era5_predict_fit("exact",subset.U(:,observedIndex),observed,target,exact,foldCfg,[]);
    hp=era5_predict_fit("hsgp",subset.U(:,observedIndex),observed,target,hsgp,foldCfg,basis);
    names=["data are finite 360-by-289";"five folds cover central test grid"; ...
        "exact objective finite";"HSGP objective finite"; ...
        "exact predictions finite";"HSGP predictions finite"];
    passed=[isequal(size(data.U),[360 289])&&all(isfinite(data.U),'all'); ...
        numel(folds)==5&&all(manifest.foldID>0)&&height(summary)==5; ...
        isfinite(exact.objective);isfinite(hsgp.objective); ...
        all(isfinite(ep.mean),'all');all(isfinite(hp.mean),'all')];
    report=table(names,passed,'VariableNames',{'check','passed'});
    disp(report);if ~all(passed),error('ERA5 preflight failed.');end
    fprintf('ERA5 preflight passed: %d/%d checks.\n',sum(passed),numel(passed));
end
