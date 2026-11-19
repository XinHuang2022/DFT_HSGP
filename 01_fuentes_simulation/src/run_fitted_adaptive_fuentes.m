function result = run_fitted_adaptive_fuentes(cfg)
%RUN_FITTED_ADAPTIVE_FUENTES Select frequency ranks at fitted parameters.
%
% Analytic full-rank Fuentes estimates are read from the trusted paired
% recovery experiment. For each replicate, the KL/site=0.01 rank profile is
% selected using fitting coordinates and the replicate-specific estimate.
% Held-site responses are not used in rank selection.

    if nargin < 1 || isempty(cfg)
        cfg = fitted_adaptive_config();
    end
    validate_inputs(cfg);
    addpath(cfg.sourceSuiteRoot);
    addpath(cfg.oracleSourceDirectory);

    loaded = load(cfg.oracleResultFile, 'result');
    oracle = loaded.result;
    sourceCfg = oracle.sourceCfg;
    sites = oracle.sites;
    allXy = sites.xy;
    configurations = oracle.configurations;
    folds = oracle.foldManifest;
    observedMask = folds.foldID ~= sourceCfg.fixedVariance.predictionFold;
    targetMask = ~observedMask;
    xyObserved = sites.xy(observedMask, :);
    xyTarget = sites.xy(targetMask, :);
    T = sourceCfg.main.T;
    nPositive = floor(T / 2) + 1;
    lambdaAll = eq9_fourier_frequencies(T);
    lambda = lambdaAll(1:nPositive);
    frequencyIndex = (0:nPositive-1).';

    parameters = readtable(cfg.fittedParameterFile, ...
        'TextType', 'string', 'VariableNamingRule', 'preserve');
    parameters = parameters(parameters.fitRegime == cfg.fitRegime & ...
        parameters.engine == cfg.fitEngine & ...
        ismember(parameters.configurationID, ...
        configurations.configurationID), :);
    trustedPrediction = readtable(cfg.trustedPredictionFile, ...
        'TextType', 'string', 'VariableNamingRule', 'preserve');

    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    outputDir = fullfile(cfg.outputRoot, ['fitted_adaptive_', stamp]);
    figureDir = fullfile(outputDir, 'figures');
    ensure_directory(outputDir);
    ensure_directory(figureDir);
    copyfile(cfg.fittedParameterFile, fullfile(outputDir, ...
        'trusted_paired_parameter_recovery.csv'));
    copyfile(cfg.trustedPredictionFile, fullfile(outputDir, ...
        'trusted_paired_prediction_rmse.csv'));
    writetable(folds, fullfile(outputDir, 'holdout_manifest.csv'));

    if cfg.useParallel
        pool = gcp('nocreate');
        if isempty(pool)
            parpool('local', cfg.nWorkers);
        elseif pool.NumWorkers ~= cfg.nWorkers
            warning('Using existing pool with %d workers, requested %d.', ...
                pool.NumWorkers, cfg.nWorkers);
        end
    end

    auditCells = cell(height(configurations), 1);
    profileCells = cell(height(configurations), 1);
    predictionCells = cell(height(configurations), 1);
    validationCells = cell(height(configurations), 1);

    fprintf(['Fitted-parameter Fuentes experiment: %d configurations, ', ...
        '%d replicates, %d frequencies, KL/site=%g.\n'], ...
        height(configurations), cfg.nReplicates, nPositive, ...
        cfg.klPerFittingSiteTolerance);

    for g = 1:height(configurations)
        configuration = configurations(g, :);
        configurationID = configuration.configurationID;
        thetaTruth = theta_from_configuration(configuration, sourceCfg);
        fitRows = parameters(parameters.configurationID == ...
            configurationID, :);
        fitRows = sortrows(fitRows, 'replicateID');
        if height(fitRows) ~= cfg.nReplicates || ...
                any(fitRows.replicateID ~= (1:cfg.nReplicates).')
            error('Expected fitted estimates for replicates 1:%d in config %d.', ...
                cfg.nReplicates, configurationID);
        end
        if any(~fitRows.convergenceCertified)
            error('Configuration %d contains uncertified parameter fits.', ...
                configurationID);
        end

        common = common_profile(oracle, configurationID, nPositive);
        oracleProfile = oracle_kl_profile(oracle, configurationID, ...
            cfg.klPerFittingSiteTolerance);
        boundaryFactor = common.boundaryFactor(1);
        [bank, validation] = make_nested_geometry_bank(allXy, ...
            xyObserved, xyTarget, boundaryFactor, cfg.candidateJ, sourceCfg);
        validation.configurationID = repmat(configurationID, ...
            height(validation), 1);
        validationCells{g} = movevars(validation, 'configurationID', ...
            'Before', 1);

        audits = cell(cfg.nReplicates, 1);
        profiles = cell(cfg.nReplicates, 1);
        predictions = cell(cfg.nReplicates, 1);
        if cfg.useParallel
            parfor r = 1:cfg.nReplicates
                [audits{r}, profiles{r}, predictions{r}] = ...
                    run_one_replicate(r, configurationID, thetaTruth, ...
                    fitRows(r, :), allXy, observedMask, targetMask, ...
                    xyObserved, xyTarget, lambda, frequencyIndex, ...
                    common, oracleProfile, bank, sourceCfg, cfg);
            end
        else
            for r = 1:cfg.nReplicates
                [audits{r}, profiles{r}, predictions{r}] = ...
                    run_one_replicate(r, configurationID, thetaTruth, ...
                    fitRows(r, :), allXy, observedMask, targetMask, ...
                    xyObserved, xyTarget, lambda, frequencyIndex, ...
                    common, oracleProfile, bank, sourceCfg, cfg);
            end
        end
        auditCells{g} = vertcat(audits{:});
        profileCells{g} = vertcat(profiles{:});
        predictionCells{g} = vertcat(predictions{:});
        fprintf('Configuration %d complete.\n', configurationID);
    end

    result.schemaVersion = 1;
    result.createdAt = char(datetime('now', ...
        'Format', 'yyyy-MM-dd HH:mm:ss'));
    result.cfg = cfg;
    result.sourceCfg = sourceCfg;
    result.outputDir = outputDir;
    result.configurations = configurations;
    result.sites = sites;
    result.foldManifest = folds;
    result.fittedParameters = parameters;
    result.frequencyAudit = vertcat(auditCells{:});
    result.fittedRankProfiles = vertcat(profileCells{:});
    result.profileFrequencySummary = summarize_frequency_profiles( ...
        result.fittedRankProfiles);
    result.profileOverallSummary = summarize_profiles( ...
        result.fittedRankProfiles);
    result.nestingValidation = vertcat(validationCells{:});
    result.predictionRmse = vertcat(predictionCells{:});
    result.predictionSummary = summarize_predictions( ...
        result.predictionRmse, cfg);
    result.trustedRegression = compare_trusted_predictions( ...
        result.predictionRmse, trustedPrediction);

    if any(result.trustedRegression.maximumAbsoluteRmseDifference > 1e-10)
        error('Trusted analytic/common prediction regression check failed.');
    end

    writetable(result.frequencyAudit, fullfile(outputDir, ...
        'fitted_frequency_audit.csv'));
    writetable(result.fittedRankProfiles, fullfile(outputDir, ...
        'fitted_rank_profiles.csv'));
    writetable(result.profileFrequencySummary, fullfile(outputDir, ...
        'rank_profile_frequency_summary.csv'));
    writetable(result.profileOverallSummary, fullfile(outputDir, ...
        'rank_profile_overall_summary.csv'));
    writetable(result.nestingValidation, fullfile(outputDir, ...
        'nesting_validation.csv'));
    writetable(result.predictionRmse, fullfile(outputDir, ...
        'prediction_rmse.csv'));
    writetable(result.predictionSummary, fullfile(outputDir, ...
        'prediction_rmse_summary.csv'));
    writetable(result.trustedRegression, fullfile(outputDir, ...
        'trusted_prediction_regression.csv'));
    save(fullfile(outputDir, 'result.mat'), 'result', '-v7.3');
    make_figures(result, figureDir, cfg);
    write_metadata(result, fullfile(outputDir, 'run_metadata.txt'));
    fprintf('Fitted-parameter experiment complete: %s\n', outputDir);
end

function validate_inputs(cfg)
    requiredFiles = {cfg.oracleResultFile, cfg.fittedParameterFile, ...
        cfg.trustedPredictionFile};
    if ~isfolder(cfg.sourceSuiteRoot)
        error('Trusted source suite not found: %s', cfg.sourceSuiteRoot);
    end
    for j = 1:numel(requiredFiles)
        if ~isfile(requiredFiles{j})
            error('Required input not found: %s', requiredFiles{j});
        end
    end
    if any(diff(cfg.candidateJ) <= 0) || ...
            any(cfg.candidateJ ~= fix(cfg.candidateJ))
        error('candidateJ must contain increasing integers.');
    end
    if cfg.klPerFittingSiteTolerance <= 0
        error('KL tolerance must be positive.');
    end
end

function theta = theta_from_configuration(row, sourceCfg)
    theta = sourceCfg.truth;
    fields = {'signalVariance', 'alpha', 'beta', 'epsilon', ...
        'noiseVariance'};
    for j = 1:numel(fields)
        theta.(fields{j}) = row.(fields{j});
    end
end

function theta = theta_from_fit(row, sourceCfg)
    theta = sourceCfg.truth;
    fields = {'signalVariance', 'alpha', 'beta', 'epsilon', ...
        'noiseVariance'};
    for j = 1:numel(fields)
        theta.(fields{j}) = row.(fields{j});
    end
end

function profile = common_profile(oracle, configurationID, nFrequency)
    source = oracle.klRankProfiles( ...
        oracle.klRankProfiles.configurationID == configurationID, :);
    commonJ = source.commonBasisPerDimension(1);
    commonM = source.commonM(1);
    boundaryFactor = oracle.frequencyAudit.boundaryFactor(find( ...
        oracle.frequencyAudit.configurationID == configurationID, 1));
    profile = table((0:nFrequency-1).', repmat(commonJ, nFrequency, 1), ...
        repmat(commonM, nFrequency, 1), repmat(boundaryFactor, ...
        nFrequency, 1), 'VariableNames', {'frequencyIndex', ...
        'basisPerDimension', 'M', 'boundaryFactor'});
end

function profile = oracle_kl_profile(oracle, configurationID, tolerance)
    source = oracle.klRankProfiles( ...
        oracle.klRankProfiles.configurationID == configurationID & ...
        abs(oracle.klRankProfiles.klPerFittingSiteTolerance - tolerance) ...
        < 1e-12, :);
    profile = table(source.frequencyIndex, ...
        source.selectedBasisPerDimension, source.selectedM, ...
        'VariableNames', {'frequencyIndex', 'basisPerDimension', 'M'});
end

function [bank, validation] = make_nested_geometry_bank(domainXy, ...
        xyObserved, xyTarget, boundaryFactor, candidateJ, sourceCfg)
    maximumJ = max(candidateJ);
    maximumBasis = make_eq9_hsgp_basis( ...
        domainXy, boundaryFactor, maximumJ, sourceCfg);
    PhiObservedMaximum = eq9_hsgp_features(xyObserved, maximumBasis);
    PhiTargetMaximum = eq9_hsgp_features(xyTarget, maximumBasis);
    entries = repmat(struct(), numel(candidateJ), 1);
    rows = cell(numel(candidateJ), 1);
    for q = 1:numel(candidateJ)
        J = candidateJ(q);
        basis = make_eq9_hsgp_basis(domainXy, boundaryFactor, J, sourceCfg);
        [found, map] = ismember(basis.indices, maximumBasis.indices, 'rows');
        if ~all(found)
            error('J=%d is not nested in maximum J=%d.', J, maximumJ);
        end
        PhiObservedDirect = eq9_hsgp_features(xyObserved, basis);
        PhiTargetDirect = eq9_hsgp_features(xyTarget, basis);
        PhiObserved = PhiObservedMaximum(:, map);
        PhiTarget = PhiTargetMaximum(:, map);
        observedDifference = max(abs(PhiObservedDirect - PhiObserved), ...
            [], 'all');
        targetDifference = max(abs(PhiTargetDirect - PhiTarget), [], 'all');
        rows{q} = table(J, J^2, maximumJ, maximumJ^2, ...
            observedDifference, targetDifference, ...
            'VariableNames', {'basisPerDimension', 'M', ...
            'maximumBasisPerDimension', 'maximumM', ...
            'maximumObservedFeatureDifference', ...
            'maximumTargetFeatureDifference'});
        entries(q).J = J;
        entries(q).M = J^2;
        entries(q).basis = basis;
        entries(q).map = map;
        entries(q).PhiObserved = PhiObserved;
        entries(q).PhiTarget = PhiTarget;
    end
    validation = vertcat(rows{:});
    if any(validation.maximumObservedFeatureDifference > 1e-12) || ...
            any(validation.maximumTargetFeatureDifference > 1e-12)
        error('Nested feature validation failed.');
    end
    bank.entries = entries;
    bank.maximumBasis = maximumBasis;
end

function [audit, profile, prediction] = run_one_replicate( ...
        replicateID, configurationID, thetaTruth, fitRow, allXy, ...
        observedMask, targetMask, xyObserved, xyTarget, lambda, ...
        frequencyIndex, common, oracleProfile, bank, sourceCfg, cfg)
    thetaFit = theta_from_fit(fitRow, sourceCfg);
    [audit, fittedProfile, weightBank] = select_fitted_profile( ...
        replicateID, configurationID, thetaFit, xyObserved, lambda, ...
        frequencyIndex, common, oracleProfile, bank, sourceCfg, cfg);
    profile = fittedProfile;

    exact = build_exact_operators(thetaFit, xyObserved, xyTarget, ...
        lambda, sourceCfg);
    commonOperators = build_hsgp_operators(thetaFit, ...
        common.basisPerDimension, bank, weightBank, lambda);
    oracleOperators = build_hsgp_operators(thetaFit, ...
        oracleProfile.basisPerDimension, bank, weightBank, lambda);
    fittedOperators = build_hsgp_operators(thetaFit, ...
        fittedProfile.selectedBasisPerDimension, bank, weightBank, ...
        lambda);

    seed = sourceCfg.recovery.predictionSeeds(replicateID);
    simulation = simulate_eq9_circulant_components( ...
        thetaTruth, allXy, sourceCfg.main.T, sourceCfg, seed);
    Yobserved = simulation.Y(:, observedMask);
    latentTarget = simulation.latent(:, targetMask);
    noisyTarget = simulation.Y(:, targetMask);
    methods = ["analytic_reference_fit"; "common_rank_hsgp_fit"; ...
        "oracle_profile_hsgp_fit"; "fitted_profile_hsgp_fit"];
    operators = {exact; commonOperators; oracleOperators; fittedOperators};
    rows = cell(numel(methods), 1);
    exactMean = reconstruct(Yobserved, exact);
    for q = 1:numel(methods)
        if q == 1
            predicted = exactMean;
        else
            predicted = reconstruct(Yobserved, operators{q});
        end
        rows{q} = table(replicateID, seed, fitRow.estimationSeed, ...
            configurationID, methods(q), ...
            sqrt(mean((predicted - latentTarget).^2, 'all')), ...
            sqrt(mean((predicted - noisyTarget).^2, 'all')), ...
            sqrt(mean((predicted - exactMean).^2, 'all')), ...
            'VariableNames', {'replicateID', 'predictionSeed', ...
            'estimationSeed', 'configurationID', 'method', ...
            'latentRmse', 'noisyRmse', ...
            'predictionMeanRmseVersusAnalyticReference'});
    end
    prediction = vertcat(rows{:});
end

function [audit, profile, weightBank] = select_fitted_profile( ...
        replicateID, configurationID, thetaFit, xyObserved, lambda, ...
        frequencyIndex, common, oracleProfile, bank, sourceCfg, cfg)
    nCandidate = numel(bank.entries);
    nFrequency = numel(lambda);
    nObserved = size(xyObserved, 1);
    Fexact = eq9_frequency_covariances( ...
        thetaFit, xyObserved, lambda, sourceCfg, false);
    weightBank = cell(nCandidate, 1);
    errors = zeros(nCandidate, nFrequency);
    jitterReference = zeros(nCandidate, nFrequency);
    jitterHsgp = zeros(nCandidate, nFrequency);
    auditCells = cell(nCandidate, 1);
    for q = 1:nCandidate
        entry = bank.entries(q);
        weights = eq9_hsgp_weights(thetaFit, lambda, ...
            entry.basis, sourceCfg);
        weightBank{q} = weights;
        rows = cell(nFrequency, 1);
        for k = 1:nFrequency
            P = Fexact(:, :, k) + thetaFit.noiseVariance * eye(nObserved);
            Qsignal = eq9_hsgp_covariance(entry.PhiObserved, ...
                entry.PhiObserved, weights(:, k), 0, true);
            Q = Qsignal + thetaFit.noiseVariance * eye(nObserved);
            endpoint = (k == 1) || (k == nFrequency);
            factor = 1 - 0.5 * endpoint;
            [kl, jitterReference(q, k), jitterHsgp(q, k)] = ...
                gaussian_kl_covariance(P, Q, factor);
            errors(q, k) = kl / nObserved;
            rows{k} = table(replicateID, configurationID, ...
                frequencyIndex(k), lambda(k), endpoint, entry.J, ...
                entry.M, errors(q, k), jitterReference(q, k), ...
                jitterHsgp(q, k), ...
                'VariableNames', {'replicateID', 'configurationID', ...
                'frequencyIndex', 'lambda', 'endpointFrequency', ...
                'basisPerDimension', 'M', 'klPerFittingSite', ...
                'referenceCholeskyJitter', 'hsgpCholeskyJitter'});
        end
        auditCells{q} = vertcat(rows{:});
    end
    audit = vertcat(auditCells{:});

    candidateJ = [bank.entries.J].';
    selectedJ = zeros(nFrequency, 1);
    selectedError = zeros(nFrequency, 1);
    persistentError = zeros(nFrequency, 1);
    thresholdMet = false(nFrequency, 1);
    for k = 1:nFrequency
        tailMaximum = flipud(cummax(flipud(errors(:, k))));
        position = find(tailMaximum <= ...
            cfg.klPerFittingSiteTolerance, 1, 'first');
        thresholdMet(k) = ~isempty(position);
        if isempty(position)
            position = nCandidate;
        end
        selectedJ(k) = candidateJ(position);
        selectedError(k) = errors(position, k);
        persistentError(k) = tailMaximum(position);
    end
    profile = table(repmat(replicateID, nFrequency, 1), ...
        repmat(configurationID, nFrequency, 1), frequencyIndex, lambda(:), ...
        selectedJ, selectedJ.^2, thresholdMet, selectedError, ...
        persistentError, common.basisPerDimension, common.M, ...
        oracleProfile.basisPerDimension, oracleProfile.M, ...
        'VariableNames', {'replicateID', 'configurationID', ...
        'frequencyIndex', 'lambda', 'selectedBasisPerDimension', ...
        'selectedM', 'thresholdMet', 'selectedKlPerFittingSite', ...
        'persistentTailMaximumKlPerFittingSite', ...
        'commonBasisPerDimension', 'commonM', ...
        'oracleBasisPerDimension', 'oracleM'});
end

function operators = build_exact_operators(theta, xyObserved, xyTarget, ...
        lambda, sourceCfg)
    Foo = eq9_frequency_covariances(theta, xyObserved, ...
        lambda, sourceCfg, true);
    Fto = eq9_cross_frequency_covariances( ...
        theta, xyTarget, xyObserved, lambda, sourceCfg);
    operators = zeros(size(xyTarget, 1), size(xyObserved, 1), numel(lambda));
    for k = 1:numel(lambda)
        [L, ~] = eq9_chol_with_jitter(Foo(:, :, k));
        operators(:, :, k) = (L' \ (L \ Fto(:, :, k)')).';
    end
end

function operators = build_hsgp_operators(theta, profileJ, bank, ...
        weightBank, lambda)
    nTarget = size(bank.entries(1).PhiTarget, 1);
    nObserved = size(bank.entries(1).PhiObserved, 1);
    operators = zeros(nTarget, nObserved, numel(lambda));
    for k = 1:numel(lambda)
        q = find([bank.entries.J] == profileJ(k), 1);
        entry = bank.entries(q);
        weights = weightBank{q}(:, k);
        factor = eq9_lowrank_factor( ...
            entry.PhiObserved, weights, theta.noiseVariance);
        Kto = eq9_hsgp_covariance(entry.PhiTarget, ...
            entry.PhiObserved, weights, 0, false);
        operators(:, :, k) = eq9_lowrank_solve(factor, Kto').';
    end
end

function value = reconstruct(Yobserved, operators)
    T = size(Yobserved, 1);
    nFrequency = size(operators, 3);
    D = fft(Yobserved, [], 1) / sqrt(T);
    Dmean = complex(zeros(T, size(operators, 1)));
    for k = 1:nFrequency
        Dmean(k, :) = (operators(:, :, k) * D(k, :).').';
    end
    for k = 2:nFrequency-1
        Dmean(T - k + 2, :) = conj(Dmean(k, :));
    end
    value = real(ifft(Dmean, [], 1) * sqrt(T));
end

function summary = summarize_frequency_profiles(profiles)
    configurations = unique(profiles.configurationID, 'stable');
    nFrequency = numel(unique(profiles.frequencyIndex));
    rows = cell(numel(configurations) * nFrequency, 1);
    row = 0;
    for g = 1:numel(configurations)
        for k = 0:nFrequency-1
            subset = profiles(profiles.configurationID == ...
                configurations(g) & profiles.frequencyIndex == k, :);
            values = subset.selectedM;
            row = row + 1;
            rows{row} = table(configurations(g), k, subset.lambda(1), ...
                median(values), local_quantile(values, 0.25), ...
                local_quantile(values, 0.75), min(values), max(values), ...
                subset.commonM(1), subset.oracleM(1), ...
                mean(values == subset.oracleM), ...
                mean(values > subset.oracleM), ...
                mean(values < subset.oracleM), ...
                'VariableNames', {'configurationID', 'frequencyIndex', ...
                'lambda', 'medianFittedM', 'q25FittedM', 'q75FittedM', ...
                'minimumFittedM', 'maximumFittedM', 'commonM', ...
                'oracleM', 'proportionEqualOracle', ...
                'proportionAboveOracle', 'proportionBelowOracle'});
        end
    end
    summary = vertcat(rows{:});
end

function summary = summarize_profiles(profiles)
    configurations = unique(profiles.configurationID, 'stable');
    rows = cell(numel(configurations), 1);
    for g = 1:numel(configurations)
        subset = profiles(profiles.configurationID == configurations(g), :);
        rows{g} = table(configurations(g), height(subset), ...
            mean(subset.selectedM), median(subset.selectedM), ...
            min(subset.selectedM), max(subset.selectedM), ...
            mean(subset.selectedM == subset.oracleM), ...
            mean(subset.selectedM > subset.oracleM), ...
            mean(subset.selectedM < subset.oracleM), ...
            sum(~subset.thresholdMet), ...
            'VariableNames', {'configurationID', ...
            'nReplicateFrequencyPairs', 'meanFittedM', 'medianFittedM', ...
            'minimumFittedM', 'maximumFittedM', ...
            'proportionEqualOracle', 'proportionAboveOracle', ...
            'proportionBelowOracle', 'numberUnresolved'});
    end
    summary = vertcat(rows{:});
end

function summary = summarize_predictions(rows, cfg)
    configurations = unique(rows.configurationID, 'stable');
    methods = unique(rows.method, 'stable');
    cells = cell(numel(configurations) * numel(methods), 1);
    rng(cfg.bootstrapSeed, 'twister');
    row = 0;
    for g = 1:numel(configurations)
        reference = rows(rows.configurationID == configurations(g) & ...
            rows.method == "analytic_reference_fit", :);
        common = rows(rows.configurationID == configurations(g) & ...
            rows.method == "common_rank_hsgp_fit", :);
        for q = 1:numel(methods)
            subset = rows(rows.configurationID == configurations(g) & ...
                rows.method == methods(q), :);
            deltaReference = subset.latentRmse - reference.latentRmse;
            deltaCommon = subset.latentRmse - common.latentRmse;
            [low, high] = bootstrap_mean_ci( ...
                subset.latentRmse, cfg.bootstrapReplicates);
            [referenceLow, referenceHigh] = bootstrap_mean_ci( ...
                deltaReference, cfg.bootstrapReplicates);
            [commonLow, commonHigh] = bootstrap_mean_ci( ...
                deltaCommon, cfg.bootstrapReplicates);
            row = row + 1;
            cells{row} = table(configurations(g), methods(q), ...
                height(subset), mean(subset.latentRmse), low, high, ...
                mean(subset.noisyRmse), mean(deltaReference), ...
                referenceLow, referenceHigh, mean(deltaCommon), ...
                commonLow, commonHigh, ...
                'VariableNames', {'configurationID', 'method', ...
                'nReplicates', 'meanLatentRmse', 'bootstrap95Low', ...
                'bootstrap95High', 'meanNoisyRmse', ...
                'meanPairedDifferenceFromAnalyticReference', ...
                'pairedReferenceDifferenceBootstrap95Low', ...
                'pairedReferenceDifferenceBootstrap95High', ...
                'meanPairedDifferenceFromCommon', ...
                'pairedCommonDifferenceBootstrap95Low', ...
                'pairedCommonDifferenceBootstrap95High'});
        end
    end
    summary = vertcat(cells{:});
end

function regression = compare_trusted_predictions(rows, trusted)
    configurations = unique(rows.configurationID, 'stable');
    methods = ["analytic_reference_fit"; "common_rank_hsgp_fit"];
    trustedEngines = ["exact"; "hsgp"];
    cells = cell(numel(configurations) * numel(methods), 1);
    row = 0;
    for g = 1:numel(configurations)
        for q = 1:numel(methods)
            current = rows(rows.configurationID == configurations(g) & ...
                rows.method == methods(q), :);
            target = trusted(trusted.configurationID == configurations(g) & ...
                trusted.thetaSource == "free_exact" & ...
                trusted.predictionEngine == trustedEngines(q), :);
            current = sortrows(current, 'replicateID');
            target = sortrows(target, 'replicateID');
            if height(current) ~= height(target) || ...
                    any(current.replicateID ~= target.replicateID)
                error('Trusted prediction rows do not align.');
            end
            row = row + 1;
            cells{row} = table(configurations(g), methods(q), ...
                max(abs(current.latentRmse - target.latentRmse)), ...
                max(abs(current.noisyRmse - target.noisyRmse)), ...
                'VariableNames', {'configurationID', 'method', ...
                'maximumAbsoluteRmseDifference', ...
                'maximumAbsoluteNoisyRmseDifference'});
        end
    end
    regression = vertcat(cells{:});
end

function [low, high] = bootstrap_mean_ci(values, nBootstrap)
    values = values(:);
    indices = randi(numel(values), numel(values), nBootstrap);
    means = sort(mean(values(indices), 1));
    low = interpolated_quantile(means, 0.025);
    high = interpolated_quantile(means, 0.975);
end

function value = interpolated_quantile(sortedValues, probability)
    position = 1 + (numel(sortedValues) - 1) * probability;
    lower = floor(position);
    upper = ceil(position);
    value = sortedValues(lower) + (position - lower) * ...
        (sortedValues(upper) - sortedValues(lower));
end

function value = local_quantile(values, probability)
    value = interpolated_quantile(sort(values(:)), probability);
end

function make_figures(result, figureDir, cfg)
    configurations = unique(result.profileFrequencySummary.configurationID, ...
        'stable');
    fig = figure('Visible', cfg.figureVisible, 'Color', 'w', ...
        'Position', [100, 100, 1250, 680]);
    tiledlayout(numel(configurations), 1, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    for g = 1:numel(configurations)
        subset = result.profileFrequencySummary( ...
            result.profileFrequencySummary.configurationID == ...
            configurations(g), :);
        nexttile;
        x = subset.frequencyIndex;
        fill([x; flipud(x)], [subset.q25FittedM; ...
            flipud(subset.q75FittedM)], [0.9290, 0.6940, 0.1250], ...
            'FaceAlpha', 0.22, 'EdgeColor', 'none');
        hold on;
        stairs(x, subset.commonM, '-', 'Color', [0, 0.4470, 0.7410], ...
            'LineWidth', 1.8);
        stairs(x, subset.oracleM, '--', 'Color', [0.20, 0.20, 0.20], ...
            'LineWidth', 1.8);
        stairs(x, subset.medianFittedM, '-', ...
            'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2.0);
        xlabel('DFT frequency index, k');
        ylabel('retained rank, M_k');
        title(sprintf('Configuration %s', configuration_label(configurations(g))));
        legend({'Fitted-profile IQR', 'Fixed common rank', ...
            'Generating-parameter profile', 'Fitted-profile median'}, ...
            'Location', 'northeast');
        grid on;
        box on;
    end
    save_figure(fig, figureDir, 'fitted_vs_oracle_rank_profiles');

    fig = figure('Visible', cfg.figureVisible, 'Color', 'w', ...
        'Position', [100, 100, 1250, 620]);
    tiledlayout(1, numel(configurations), 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    for g = 1:numel(configurations)
        subset = result.fittedRankProfiles( ...
            result.fittedRankProfiles.configurationID == ...
            configurations(g), :);
        subset = sortrows(subset, {'replicateID', 'frequencyIndex'});
        replicateID = unique(subset.replicateID, 'stable');
        nFrequency = numel(unique(subset.frequencyIndex));
        values = reshape(subset.selectedM, nFrequency, ...
            numel(replicateID)).';
        nexttile;
        imagesc(0:size(values, 2)-1, replicateID, values);
        axis xy;
        xlabel('DFT frequency index, k');
        ylabel('fitted replicate');
        title(sprintf('Configuration %s', configuration_label(configurations(g))));
        colorbar;
    end
    save_figure(fig, figureDir, 'fitted_rank_profile_heatmaps');

    fig = figure('Visible', cfg.figureVisible, 'Color', 'w', ...
        'Position', [100, 100, 1150, 500]);
    tiledlayout(1, numel(configurations), 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    for g = 1:numel(configurations)
        subset = result.predictionSummary( ...
            result.predictionSummary.configurationID == ...
            configurations(g), :);
        nexttile;
        x = 1:height(subset);
        errorbar(x, subset.meanLatentRmse, ...
            subset.meanLatentRmse - subset.bootstrap95Low, ...
            subset.bootstrap95High - subset.meanLatentRmse, ...
            'o', 'LineWidth', 1.5, 'MarkerFaceColor', [0, 0.4470, 0.7410]);
        xticks(x);
        xticklabels(pretty_labels(subset.method));
        xtickangle(18);
        ylabel('held-site latent RMSE');
        title(sprintf('Configuration %s', configuration_label(configurations(g))));
        grid on;
    end
    save_figure(fig, figureDir, 'fitted_adaptive_prediction_rmse');
end

function label = configuration_label(configurationID)
    switch configurationID
        case 5
            label = 'A';
        case 21
            label = 'B';
        otherwise
            label = string(configurationID);
    end
end

function labels = pretty_labels(methods)
    labels = cell(numel(methods), 1);
    for q = 1:numel(methods)
        switch string(methods(q))
            case "analytic_reference_fit"
                labels{q} = 'analytic fit';
            case "common_rank_hsgp_fit"
                labels{q} = 'common rank';
            case "oracle_profile_hsgp_fit"
                labels{q} = 'oracle profile';
            case "fitted_profile_hsgp_fit"
                labels{q} = 'fitted profile';
            otherwise
                labels{q} = char(methods(q));
        end
    end
end

function save_figure(fig, directory, name)
    exportgraphics(fig, fullfile(directory, [name, '.pdf']), ...
        'ContentType', 'vector');
    exportgraphics(fig, fullfile(directory, [name, '.png']), ...
        'Resolution', 300);
    savefig(fig, fullfile(directory, [name, '.fig']));
    close(fig);
end

function write_metadata(result, filename)
    fid = fopen(filename, 'w');
    if fid < 0
        error('Could not create metadata file.');
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, 'version=%s\n', result.cfg.version);
    fprintf(fid, 'created_at=%s\n', result.createdAt);
    fprintf(fid, 'fit_regime=%s\n', result.cfg.fitRegime);
    fprintf(fid, 'fit_engine=%s\n', result.cfg.fitEngine);
    fprintf(fid, 'kl_per_fitting_site_tolerance=%.17g\n', ...
        result.cfg.klPerFittingSiteTolerance);
    fprintf(fid, 'candidate_J=%s\n', mat2str(result.cfg.candidateJ));
    fprintf(fid, 'replicates=%d\n', result.cfg.nReplicates);
    fprintf(fid, 'parameter_file_sha256=%s\n', ...
        result.cfg.parameterFileSha256);
    fprintf(fid, 'prediction_file_sha256=%s\n', ...
        result.cfg.predictionFileSha256);
    fprintf(fid, ['selection_data=fitting coordinates and analytic ', ...
        'full-rank parameter estimate; held-site responses excluded\n']);
    fprintf(fid, ['evaluation_data=trusted independent prediction seeds ', ...
        'and fixed whole-site holdout\n']);
end

function ensure_directory(pathValue)
    if ~isfolder(pathValue)
        mkdir(pathValue);
    end
end
