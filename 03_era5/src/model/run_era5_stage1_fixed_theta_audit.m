function result = run_era5_stage1_fixed_theta_audit(cfg)
%RUN_ERA5_STAGE1_FIXED_THETA_AUDIT Fit exact once, then audit all J values.

    if nargin < 1, cfg = era5_northsea_config(); end
    ensure_directory(cfg.paths.stage1Dir);
    ensure_directory(fullfile(cfg.paths.stage1Dir, 'tasks'));
    ensure_directory(fullfile(cfg.paths.stage1Dir, 'figures'));
    contract = load_era5_wind_data(cfg);
    subset = make_era5_subset(contract, 'full17');
    [folds, foldManifest, foldSummary] = make_era5_cv_folds( ...
        subset, cfg.cv.holdoutPattern);
    fold = folds(cfg.screening.foldID);
    targetIndex = fold.targetIndex;
    observedIndex = setdiff((1:size(subset.U,2))', targetIndex, 'stable');
    [observed, target] = era5_fold_geometry(subset, observedIndex, ...
        targetIndex, cfg);
    Yobserved = subset.U(:, observedIndex); Ytarget = subset.U(:, targetIndex);
    [foldCfg, varianceReference] = era5_fold_config(cfg, Yobserved, ...
        observed.xyKm, cfg.stage1.optimizer);
    referenceFile = fullfile(cfg.paths.stage1Dir, 'exact_reference.mat');
    if isfile(referenceFile)
        S = load(referenceFile, 'reference'); reference = S.reference;
        fprintf('Stage 1 exact fit loaded from checkpoint.\n');
    else
        fprintf('Stage 1 exact fit starting on screening fold %d.\n', fold.foldID);
        seed = cfg.seed.optimizer + 1009 * fold.foldID;
        separableFile = fullfile(cfg.paths.stage1Dir, 'exact_separable_fit.mat');
        if isfile(separableFile)
            Sfit = load(separableFile, 'separableFit');
            separableFit = Sfit.separableFit;
            fprintf('Stage 1 exact separable fit loaded from checkpoint.\n');
        else
            separableFit = fit_era5_eq9_model(Yobserved, ...
                observed.xyCovariance, observed.design, foldCfg, 1, ...
                cfg.stage1.optimizer, seed + 1, [], []);
            save(separableFile, 'separableFit', '-v7.3');
        end
        fullFile = fullfile(cfg.paths.stage1Dir, 'exact_full_fit.mat');
        if isfile(fullFile)
            Sfit = load(fullFile, 'exactFit'); exactFit = Sfit.exactFit;
            fprintf('Stage 1 exact full fit loaded from checkpoint.\n');
        else
            initial = separableFit.theta;
            initial.epsilon = min(foldCfg.bounds.epsilon(2), ...
                max(foldCfg.bounds.epsilon(1), foldCfg.start.epsilon));
            fullOptimizer = cfg.stage1.optimizer;
            fullOptimizer.additionalStartZ = eq9_pack_parameters(initial, foldCfg);
            exactFit = fit_era5_eq9_model(Yobserved, ...
                observed.xyCovariance, observed.design, foldCfg, [], ...
                fullOptimizer, seed + 2, [], initial);
            save(fullFile, 'exactFit', '-v7.3');
        end
        exactPrediction = era5_predict_fit("exact", Yobserved, observed, ...
            target, exactFit, foldCfg, []);
        exactPoint = era5_point_metrics(Ytarget(:), exactPrediction.mean(:));
        lambdaAll = eq9_fourier_frequencies(size(Yobserved,1));
        selected = unique(round(linspace(1, numel(lambdaAll), ...
            cfg.hsgp.auditFrequencyCount)));
        auditLambda = lambdaAll(selected);
        exactSignal = eq9_frequency_covariances(exactFit.theta, ...
            observed.xyCovariance, auditLambda, foldCfg, false);
        fullDomainXy = (subset.xyKm - mean(subset.xyKm,1)) ./ ...
            cfg.spatial.coordinateScaleKm;
        reference = struct('version',string(cfg.version), 'subset',subset, ...
            'fold',fold, 'foldCfg',foldCfg, 'observedIndex',observedIndex, ...
            'targetIndex',targetIndex, 'observed',observed, 'target',target, ...
            'Yobserved',Yobserved, 'Ytarget',Ytarget, ...
            'separableFit',separableFit, 'exactFit',exactFit, ...
            'exactPrediction',exactPrediction, 'exactRmse',exactPoint.rmse, ...
            'auditLambda',auditLambda, 'exactSignal',exactSignal, ...
            'fullDomainXy',fullDomainXy, 'varianceReference',varianceReference);
        save(referenceFile, 'reference', '-v7.3');
        fprintf('Stage 1 exact reference saved.\n');
    end
    writetable(subset.siteTable, fullfile(cfg.paths.stage1Dir, ...
        'selected_subset_manifest.csv'));
    writetable(foldManifest, fullfile(cfg.paths.stage1Dir, 'cv_fold_manifest.csv'));
    writetable(foldSummary, fullfile(cfg.paths.stage1Dir, 'cv_fold_summary.csv'));
    Jvalues = cfg.hsgp.candidateJ(:);
    tasks = cell(numel(Jvalues),1);
    [parallelUsed,poolWorkers] = era5_open_pool(cfg,numel(Jvalues));
    if parallelUsed
        parfor q = 1:numel(Jvalues)
            tasks{q} = era5_stage1_j_task(cfg,referenceFile,Jvalues(q));
        end
    else
        for q = 1:numel(Jvalues)
            tasks{q} = era5_stage1_j_task(cfg,referenceFile,Jvalues(q));
        end
    end
    summaryParts = cellfun(@(x)x.summary,tasks,'UniformOutput',false);
    frequencyParts = cellfun(@(x)x.frequency,tasks,'UniformOutput',false);
    summary = vertcat(summaryParts{:});
    frequency = vertcat(frequencyParts{:});
    passIndex = find(summary.passesAudit,1,'first');
    if isempty(passIndex), recommendedJ = NaN; else, recommendedJ = summary.J(passIndex); end
    writetable(summary, fullfile(cfg.paths.stage1Dir,'hsgp_J_audit_summary.csv'));
    writetable(frequency, fullfile(cfg.paths.stage1Dir,'hsgp_J_frequency_diagnostics.csv'));
    runSummary = table(string(cfg.version),height(subset.siteTable),fold.foldID, ...
        numel(Jvalues),recommendedJ,parallelUsed,poolWorkers, ...
        'VariableNames',{'version','nSites','screeningFoldID','nJValues', ...
        'recommendedJ','parallelUsed','poolWorkers'});
    writetable(runSummary,fullfile(cfg.paths.stage1Dir,'run_summary.csv'));
    make_era5_stage1_figures(summary,cfg);
    save(fullfile(cfg.paths.stage1Dir,'complete_result.mat'),'cfg','reference', ...
        'summary','frequency','recommendedJ','runSummary','-v7.3');
    result = struct('summary',summary,'recommendedJ',recommendedJ, ...
        'runSummary',runSummary);
end
