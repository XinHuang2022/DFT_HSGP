function report = preflight_era5_package(experimentRoot)
%PREFLIGHT_ERA5_PACKAGE Validate local data, folds, and model interfaces.

    cfg = era5_northsea_config();
    contract = load_era5_wind_data(cfg);
    subset = make_era5_subset(contract, 'full17');
    [folds, manifest] = make_era5_cv_folds(subset, 'balanced5');

    retainedPath = fullfile(experimentRoot, 'input', 'retained', ...
        'fixed_rank', 'stage3_paired_cv', 'cv_fold_manifest.csv');
    retained = readtable(retainedPath, 'VariableNamingRule', 'preserve');
    generatedFold = zeros(height(subset.siteTable), 1);
    for f = 1:numel(folds)
        generatedFold(folds(f).targetIndex) = folds(f).foldID;
    end
    sameStationOrder = isequal(string(retained.stationID), ...
        string(subset.siteTable.stationID));
    sameFolds = isequal(retained.foldID, generatedFold);
    interfaceReport = run_era5_preflight();

    check = ["finite 360-by-289 extracted matrix"; ...
        "retained station order"; "retained five-fold assignment"; ...
        "fixed-model interface checks"];
    passed = [isequal(size(contract.U), [360 289]) && ...
        all(isfinite(contract.U), 'all'); sameStationOrder; sameFolds; ...
        all(interfaceReport.passed)];
    detail = [string(sprintf('%d-by-%d', size(contract.U))); ...
        string(sprintf('%d sites', height(subset.siteTable))); ...
        string(sprintf('%d folds, %d assigned sites', numel(folds), ...
        sum(manifest.foldID > 0))); ...
        string(sprintf('%d interface checks', height(interfaceReport)))];
    report = table(check, passed, detail);
    disp(report);
    if ~all(passed), error('ERA5 reproducibility-package preflight failed.'); end
end
