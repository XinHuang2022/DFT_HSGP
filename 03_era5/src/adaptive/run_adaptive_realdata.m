function result = run_adaptive_realdata(dataset, override)
%RUN_ADAPTIVE_REALDATA Run Algorithm 1 on all prescribed spatial folds.

    if nargin < 2, override = struct(); end
    C = adaptive_experiment_config(dataset);
    names = fieldnames(override);
    for j = 1:numel(names), C.(names{j}) = override.(names{j}); end
    if ~isfolder(C.sourceRoot)
        error('Trusted source directory not found: %s', C.sourceRoot);
    end
    addpath(C.sourceRoot);
    ensure_local_directory(C.outputRoot);
    copy_fixed_comparator(C);

    [cfg, subset, folds, adapter] = prepare_dataset(C);
    if C.smoke
        cfg.hsgp.optimizer.numStarts=1;
        cfg.hsgp.optimizer.maxIterations=2;
        cfg.hsgp.optimizer.maxFunctionEvaluations=8;
        cfg.hsgp.optimizer.selectivePolishingEnabled=false;
    end
    folds = folds(ismember([folds.foldID], C.runFolds));
    writetable(subset.siteTable, fullfile(C.outputRoot, ...
        'selected_subset_manifest.csv'));
    write_fold_manifest(folds, subset, C.outputRoot);
    save(fullfile(C.outputRoot, 'configuration.mat'), 'C', 'cfg', '-v7.3');

    runTimer = tic;
    tasks = cell(numel(folds), 1);
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool(min(C.nWorkers, numel(folds))); end
    fprintf('\n%s adaptive HSGP: %d folds on %d workers.\n', ...
        upper(char(C.dataset)), numel(folds), pool.NumWorkers);
    parfor f = 1:numel(folds)
        tasks{f} = run_adaptive_fold(C, cfg, subset, folds(f), adapter);
    end
    result = aggregate_adaptive_results(C, cfg, subset, folds, tasks, ...
        toc(runTimer));
    fprintf('\n%s adaptive HSGP complete in %.1f minutes.\nResults: %s\n', ...
        upper(char(C.dataset)), result.runSummary.elapsedSeconds / 60, ...
        C.outputRoot);
end

function [cfg, subset, folds, A] = prepare_dataset(C)
    if C.dataset == "era5"
        cfg = era5_northsea_config();
        cfg.suiteRoot = C.sourceRoot;
        cfg.paths.dataDir = C.dataDir;
        contract = load_era5_wind_data(cfg);
        subset = make_era5_subset(contract, 'full17');
        [folds, ~, ~] = make_era5_cv_folds(subset, 'balanced5');
        A.prefix = "era5";
        A.makeDesign = @make_era5_mean_design;
        A.prepare = @prepare_era5_whittle_data;
        A.features = @era5_hsgp_features;
        A.weights = @era5_hsgp_weights;
        A.factor = @era5_hsgp_factor;
        A.solve = @era5_hsgp_solve;
        A.covariance = @era5_hsgp_covariance;
        A.pack = @era5_pack_constrained_parameters;
        A.unpack = @era5_unpack_constrained_parameters;
        A.evaluateMean = @evaluate_era5_mean;
        A.fitFixed = @fit_era5_eq9_hsgp;
        A.empiricalVariance = @era5_empirical_residual_variance;
    else
        cfg = pacific_hsgp_config('full17_cv5');
        cfg.suiteRoot = C.sourceRoot;
        cfg.paths.dataDir = fullfile(C.sourceRoot, 'Data');
        % Match the boundary factor certified and used by the trusted J=48
        % Pacific fit.  The source config stores the audit ladder rather
        % than the selected scalar.
        cfg.hsgp.boundaryFactor = 4;
        contract = load_pacific_wind_data(cfg);
        subset = make_pacific_subset(contract, 'full17');
        [folds, ~, ~] = make_pacific_cv_folds(subset, 'balanced5');
        A.prefix = "pacific";
        A.makeDesign = @make_pacific_mean_design;
        A.prepare = @prepare_pacific_whittle_data;
        A.features = @pacific_hsgp_features;
        A.weights = @pacific_hsgp_weights;
        A.factor = @pacific_hsgp_factor;
        A.solve = @pacific_hsgp_solve;
        A.covariance = @pacific_hsgp_covariance;
        A.pack = @pacific_pack_constrained_parameters;
        A.unpack = @pacific_unpack_constrained_parameters;
        A.evaluateMean = @evaluate_pacific_mean;
        A.fitFixed = @fit_pacific_eq9_hsgp;
        A.empiricalVariance = @pacific_empirical_residual_variance;
    end
end

function write_fold_manifest(folds, subset, outputRoot)
    station = string(subset.siteTable.stationID);
    foldID = zeros(numel(station), 1);
    for f = 1:numel(folds), foldID(folds(f).targetIndex) = folds(f).foldID; end
    manifest = table(station, foldID, 'VariableNames', ...
        {'stationID','foldID'});
    writetable(manifest, fullfile(outputRoot, 'cv_fold_manifest.csv'));
end

function copy_fixed_comparator(C)
    files = {'method_summary.csv','task_timing.csv','fit_diagnostics.csv', ...
        'run_summary.csv'};
    destination = fullfile(C.outputRoot, 'fixed_rank_comparator');
    ensure_local_directory(destination);
    for j = 1:numel(files)
        source = fullfile(C.fixedResultRoot, files{j});
        if isfile(source)
            copyfile(source, fullfile(destination, files{j}));
        end
    end
end

function ensure_local_directory(pathValue)
    if ~isfolder(pathValue), mkdir(pathValue); end
end
