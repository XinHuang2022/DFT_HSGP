function report = preflight_pacific_package(experimentRoot)
%PREFLIGHT_PACIFIC_PACKAGE Validate local data, folds, and interfaces.

    cfg = pacific_hsgp_config('full17_cv5');
    contract = load_pacific_wind_data(cfg);
    subset = make_pacific_subset(contract, 'full17');
    [folds, manifest] = make_pacific_cv_folds(subset, 'balanced5');

    retainedPath = fullfile(experimentRoot, 'input', 'retained', ...
        'fixed_rank', 'cv_fold_manifest.csv');
    retained = readtable(retainedPath, 'VariableNamingRule', 'preserve');
    generatedFold = zeros(height(subset.siteTable), 1);
    for f = 1:numel(folds)
        generatedFold(folds(f).targetIndex) = folds(f).foldID;
    end
    sameStationOrder = isequal(string(retained.stationID), ...
        string(subset.siteTable.stationID));
    sameFolds = isequal(retained.foldID, generatedFold);
    interfaceReport = run_pacific_hsgp_preflight('full17_cv5');

    check = ["finite 480-by-289 wind matrix"; ...
        "retained station order"; "retained five-fold assignment"; ...
        "fixed and HSGP interface checks"];
    passed = [isequal(size(contract.U), [480 289]) && ...
        all(isfinite(contract.U), 'all'); sameStationOrder; sameFolds; ...
        all(interfaceReport.base.passed) && all(interfaceReport.hsgp.passed)];
    detail = [string(sprintf('%d-by-%d', size(contract.U))); ...
        string(sprintf('%d sites', height(subset.siteTable))); ...
        string(sprintf('%d folds, %d assigned sites', numel(folds), ...
        sum(manifest.foldID > 0))); ...
        string(sprintf('%d checks', height(interfaceReport.base) + ...
        height(interfaceReport.hsgp)))];
    report = table(check, passed, detail);
    disp(report);
    if ~all(passed)
        error('Pacific reproducibility-package preflight failed.');
    end
end
