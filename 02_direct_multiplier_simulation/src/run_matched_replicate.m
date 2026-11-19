function result = run_matched_replicate(C)
%RUN_MATCHED_REPLICATE Compare four methods on one matched realization.
%
% The numerical Hankel covariance is intentionally excluded from parameter
% initialization, reference stabilization, KL rank selection, and adaptive
% refitting. It is retained only for simulation and external evaluation.

    if nargin < 1 || isempty(C)
        C = matched_config('smoke');
    end
    validate_config(C);
    addpath(C.baseSourceRoot);

    totalTimer = tic;
    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    outputDir = fullfile(C.outputRoot, ['algorithm_b_', C.profile, '_', stamp]);
    figureDir = fullfile(outputDir, 'figures');
    ensure_directory(outputDir);
    ensure_directory(figureDir);

    Dsgn = design(C);
    observed = Dsgn.train(:);
    target = Dsgn.test(:);
    xo = Dsgn.xy(observed, :);
    xt = Dsgn.xy(target, :);
    T = C.T;
    K = floor(T / 2) + 1;
    lambdaAll = freqs(T);
    lambda = lambdaAll(1:K);
    frequencyIndex = (0:K-1).';
    manifest = table(Dsgn.site, Dsgn.xy(:,1), Dsgn.xy(:,2), ...
        ismember(Dsgn.site, observed), ismember(Dsgn.site, target), ...
        'VariableNames', {'site','x','y','fitting','heldOut'});
    writetable(manifest, fullfile(outputDir, 'site_manifest.csv'));

    referenceTableTimer = tic;
    R = rtable(C);
    referenceTableSeconds = toc(referenceTableTimer);

    basisTimer = tic;
    maximumM = max([C.candidateM(:); C.referenceM(:); C.Minit]);
    master = mkbasis(Dsgn.xy, C.c, maximumM);
    master.PhiObserved = features(xo, master);
    master.PhiTarget = features(xt, master);
    validate_nested_bases(Dsgn.xy, xo, master, ...
        unique([C.Minit,C.candidateM,C.referenceM]));
    basisBankSeconds = toc(basisTimer);

    simulationTimer = tic;
    simulation = sim_ref(C.truth, Dsgn.xy, T, R, C.simulationSeed);
    Yobserved = simulation.Y(:, observed);
    latentTarget = simulation.X(:, target);
    noisyTarget = simulation.Y(:, target);
    Dpositive = fft(Yobserved, [], 1) / sqrt(T);
    Dpositive = Dpositive(1:K, :);
    simulationSeconds = toc(simulationTimer);

    initialBasis = subset_basis(master, C.Minit, Dsgn.xy);
    initialFitTimer = tic;
    initialFit = fit_model(Yobserved, xo, C, R, 'hsgp', ...
        initialBasis, C.optimizerSeed);
    initialFitStageSeconds = toc(initialFitTimer);

    initialReferenceTimer = tic;
    [initialReferenceStability, initialReferenceChoice] = ...
        stabilize_reference(initialFit.theta, master, lambda, ...
        frequencyIndex, C, "initial_fit");
    initialReferenceSeconds = toc(initialReferenceTimer);
    require_stable_reference(initialReferenceChoice, 'initial fit');

    initialAuditTimer = tic;
    [initialAudit, initialProfile] = rank_audit_and_select( ...
        initialFit.theta, master, lambda, frequencyIndex, ...
        initialReferenceChoice.selectedReferenceM, C, ...
        "initial_fit", 0);
    initialAuditSeconds = toc(initialAuditTimer);
    require_resolved_profile(initialProfile, 'initial fit');

    currentM = initialProfile.selectedM;
    adaptiveStartTheta = initialFit.theta;
    adaptiveFitSeconds = 0;
    postFitReferenceSeconds = 0;
    postFitAuditSeconds = 0;
    iterationCells = cell(C.maxAdaptiveIterations, 1);
    postReferenceCells = cell(C.maxAdaptiveIterations, 1);
    postAuditCells = cell(C.maxAdaptiveIterations, 1);
    postProfileCells = cell(C.maxAdaptiveIterations, 1);
    adaptiveFit = struct();
    finalReferenceChoice = initialReferenceChoice;
    stabilized = false;
    nAdaptiveFits = 0;
    for cycle = 1:C.maxAdaptiveIterations
        nAdaptiveFits = cycle;
        if cycle == 1
            numberStarts = C.starts;
        else
            numberStarts = C.correctionStarts;
        end
        fitTimer = tic;
        adaptiveFit = fit_variable_rank(Dpositive, adaptiveStartTheta, ...
            currentM, master, lambda, C, C.optimizerSeed + 100 * cycle, ...
            numberStarts);
        fitSeconds = toc(fitTimer);
        adaptiveFitSeconds = adaptiveFitSeconds + fitSeconds;

        referenceTimer = tic;
        [cycleReference, cycleReferenceChoice] = stabilize_reference( ...
            adaptiveFit.theta, master, lambda, frequencyIndex, C, ...
            "adaptive_refit_cycle_" + cycle);
        referenceSeconds = toc(referenceTimer);
        postFitReferenceSeconds = postFitReferenceSeconds + referenceSeconds;
        require_stable_reference(cycleReferenceChoice, ...
            sprintf('adaptive cycle %d', cycle));

        auditTimer = tic;
        [cycleAudit, requiredProfile] = rank_audit_and_select( ...
            adaptiveFit.theta, master, lambda, frequencyIndex, ...
            cycleReferenceChoice.selectedReferenceM, C, ...
            "adaptive_refit_verification", cycle);
        auditSeconds = toc(auditTimer);
        postFitAuditSeconds = postFitAuditSeconds + auditSeconds;
        require_resolved_profile(requiredProfile, ...
            sprintf('adaptive cycle %d', cycle));

        enlargedM = max(currentM, requiredProfile.selectedM);
        numberEnlarged = sum(enlargedM > currentM);
        stableNow = numberEnlarged == 0;
        iterationCells{cycle} = table(cycle, fitSeconds, ...
            referenceSeconds, auditSeconds, adaptiveFit.nll, ...
            adaptiveFit.exitFlag, fit_eval_count(adaptiveFit), ...
            numberStarts, ...
            min(currentM), median(currentM), mean(currentM), max(currentM), ...
            numberEnlarged, stableNow, ...
            'VariableNames', {'cycle','fitSeconds', ...
            'referenceStabilitySeconds','auditSeconds','objective', ...
            'exitFlag','objectiveEvaluations','numberStarts','minimumMUsed', ...
            'medianMUsed','meanMUsed','maximumMUsed', ...
            'numberFrequenciesEnlarged','profileStabilized'});
        postReferenceCells{cycle} = cycleReference;
        postAuditCells{cycle} = cycleAudit;
        postProfileCells{cycle} = requiredProfile;
        finalReferenceChoice = cycleReferenceChoice;
        if stableNow
            stabilized = true;
            break;
        end
        currentM = enlargedM;
        adaptiveStartTheta = adaptiveFit.theta;
    end
    if ~stabilized
        error('Algorithm-B profile did not stabilize after %d fits.', ...
            C.maxAdaptiveIterations);
    end

    finalAudit = postAuditCells{nAdaptiveFits};
    finalProfile = profile_at_used_ranks(finalAudit, currentM, ...
        frequencyIndex, lambda, C, "adaptive_final_used", ...
        nAdaptiveFits);
    require_resolved_profile(finalProfile, 'final adaptive fit');

    fixedReferenceM = finalReferenceChoice.selectedReferenceM;
    fixedReferenceTimer = tic;
    fixedReferenceBasis = subset_basis(master, fixedReferenceM, Dsgn.xy);
    fixedReferenceFit = fit_model(Yobserved, xo, C, R, 'hsgp', ...
        fixedReferenceBasis, C.optimizerSeed);
    fixedReferenceFitStageSeconds = toc(fixedReferenceTimer);

    fixedReferenceValidationTimer = tic;
    [fixedReferenceStability, fixedReferenceChoice] = ...
        stabilize_reference(fixedReferenceFit.theta, master, lambda, ...
        frequencyIndex, C, "fixed_reference_fit");
    fixedReferenceValidationSeconds = toc(fixedReferenceValidationTimer);
    require_stable_reference(fixedReferenceChoice, 'fixed reference fit');
    if fixedReferenceChoice.selectedReferenceM ~= fixedReferenceM
        error(['The fixed-reference fit requires M_ref=%d rather than ', ...
            'the prespecified stabilized value %d.'], ...
            fixedReferenceChoice.selectedReferenceM, fixedReferenceM);
    end

    numericalFitTimer = tic;
    numericalFit = fit_model(Yobserved, xo, C, R, 'exact', [], ...
        C.optimizerSeed);
    numericalFitStageSeconds = toc(numericalFitTimer);

    predictionRows = cell(5,1);
    predictionSeconds = zeros(5,1);
    predictionMethods = ["numerical_inversion_fit"; ...
        "fixed_screened_rank_fit"; "fixed_reference_rank_fit"; ...
        "algorithm_B_adaptive_fit"; ...
        "fixed_reference_at_adaptive_theta"];
    matchedReferenceFit = struct('theta',adaptiveFit.theta);
    predictionFits = {numericalFit, initialFit, fixedReferenceFit, ...
        adaptiveFit, ...
        matchedReferenceFit};
    predictionProfiles = {[], repmat(C.Minit,K,1), ...
        repmat(fixedReferenceM,K,1), ...
        finalProfile.selectedM, repmat(fixedReferenceM,K,1)};
    predictionEngines = ["exact";"hsgp";"hsgp";"hsgp";"hsgp"];
    predictions = cell(5,1);
    for q = 1:5
        timer = tic;
        predictions{q} = predict_positive(Yobserved, xo, xt, ...
            predictionFits{q}.theta, predictionEngines(q), ...
            predictionProfiles{q}, master, lambda, R, C);
        predictionSeconds(q) = toc(timer);
        predictionRows{q} = table(predictionMethods(q), ...
            sqrt(mean((predictions{q} - latentTarget).^2, 'all')), ...
            sqrt(mean((predictions{q} - noisyTarget).^2, 'all')), ...
            predictionSeconds(q), ...
            'VariableNames', {'method','latentRmse','noisyRmse', ...
            'predictionSeconds'});
    end
    predictionMetrics = vertcat(predictionRows{:});

    externalAuditTimer = tic;
    externalAudit = external_numerical_audit(adaptiveFit.theta, master, ...
        xo, lambda, frequencyIndex, finalReferenceChoice.selectedReferenceM, ...
        finalProfile.selectedM, C.Minit, R, C);
    externalAuditSeconds = toc(externalAuditTimer);

    parameters = [parameter_row("truth", C.truth, NaN, NaN, NaN); ...
        parameter_row("initialization_fit", initialFit.theta, ...
            initialFit.nll, initialFit.seconds, fit_eval_count(initialFit)); ...
        parameter_row("fixed_reference_rank_fit", ...
            fixedReferenceFit.theta, fixedReferenceFit.nll, ...
            fixedReferenceFit.seconds, fit_eval_count(fixedReferenceFit)); ...
        parameter_row("numerical_inversion_fit", numericalFit.theta, ...
            numericalFit.nll, numericalFit.seconds, ...
            fit_eval_count(numericalFit)); ...
        parameter_row("algorithm_B_adaptive_fit", adaptiveFit.theta, ...
            adaptiveFit.nll, adaptiveFit.seconds, ...
            fit_eval_count(adaptiveFit))];
    fitStatus = table(["fixed_screened_rank_fit"; ...
        "fixed_reference_rank_fit"; "numerical_inversion_fit"; ...
        "algorithm_B_adaptive_fit"], ...
        [initialFit.exitFlag; fixedReferenceFit.exitFlag; ...
        numericalFit.exitFlag; adaptiveFit.exitFlag], ...
        'VariableNames', {'method','exitFlag'});

    objectiveTimings = benchmark_objectives(Yobserved, Dpositive, xo, ...
        initialFit.theta, fixedReferenceFit.theta, adaptiveFit.theta, ...
        numericalFit.theta, initialBasis, fixedReferenceBasis, ...
        finalProfile.selectedM, master, lambda, R, C);

    adaptivePredictionSeconds = predictionMetrics.predictionSeconds( ...
        predictionMetrics.method == "algorithm_B_adaptive_fit");
    screenedFixedPredictionSeconds = predictionMetrics.predictionSeconds( ...
        predictionMetrics.method == "fixed_screened_rank_fit");
    fixedReferencePredictionSeconds = predictionMetrics.predictionSeconds( ...
        predictionMetrics.method == "fixed_reference_rank_fit");
    numericalPredictionSeconds = predictionMetrics.predictionSeconds( ...
        predictionMetrics.method == "numerical_inversion_fit");
    algorithmBFirstUseSeconds = basisBankSeconds + ...
        initialFitStageSeconds + initialReferenceSeconds + ...
        initialAuditSeconds + adaptiveFitSeconds + ...
        postFitReferenceSeconds + postFitAuditSeconds + ...
        adaptivePredictionSeconds;
    % Common reported boundary: parameter fitting, any method-specific rank
    % selection, and held-site reconstruction. Shared data generation,
    % numerical-reference table construction, and basis-bank construction are
    % recorded separately and excluded from all four compared totals.
    adaptiveComparableTotalSeconds = initialFitStageSeconds + ...
        initialReferenceSeconds + initialAuditSeconds + ...
        adaptiveFitSeconds + postFitReferenceSeconds + ...
        postFitAuditSeconds + adaptivePredictionSeconds;
    screenedFixedTotalSeconds = initialFitStageSeconds + ...
        screenedFixedPredictionSeconds;
    fixedReferenceTotalSeconds = fixedReferenceFitStageSeconds + ...
        fixedReferencePredictionSeconds;
    numericalTotalSeconds = numericalFitStageSeconds + ...
        numericalPredictionSeconds;
    totalSeconds = toc(totalTimer);
    stageTimings = table(referenceTableSeconds, simulationSeconds, ...
        basisBankSeconds, initialFitStageSeconds, initialFit.seconds, ...
        initialReferenceSeconds, initialAuditSeconds, ...
        adaptiveFitSeconds, postFitReferenceSeconds, postFitAuditSeconds, ...
        fixedReferenceFitStageSeconds, fixedReferenceFit.seconds, ...
        fixedReferenceValidationSeconds, ...
        numericalFitStageSeconds, numericalFit.seconds, ...
        adaptivePredictionSeconds, screenedFixedPredictionSeconds, ...
        fixedReferencePredictionSeconds, ...
        numericalPredictionSeconds, externalAuditSeconds, ...
        algorithmBFirstUseSeconds, adaptiveComparableTotalSeconds, ...
        screenedFixedTotalSeconds, fixedReferenceTotalSeconds, ...
        numericalTotalSeconds, totalSeconds, ...
        'VariableNames', {'numericalReferenceTableSeconds', ...
        'simulationSeconds','basisBankSeconds', ...
        'initialFitStageWallSeconds','initialFitOptimizationSeconds', ...
        'initialReferenceStabilitySeconds','initialRankAuditSeconds', ...
        'adaptiveFitSeconds','postFitReferenceStabilitySeconds', ...
        'postFitRankAuditSeconds','fixedReferenceFitStageWallSeconds', ...
        'fixedReferenceFitOptimizationSeconds', ...
        'fixedReferenceValidationSeconds', ...
        'numericalFitStageWallSeconds','numericalFitOptimizationSeconds', ...
        'adaptivePredictionSeconds','screenedFixedPredictionSeconds', ...
        'fixedReferencePredictionSeconds', ...
        'numericalPredictionSeconds','externalNumericalAuditSeconds', ...
        'algorithmBFirstUseSeconds','adaptiveComparableTotalSeconds', ...
        'screenedFixedTotalSeconds','fixedReferenceTotalSeconds', ...
        'numericalReferenceTotalSeconds', ...
        'entireScriptSeconds'});

    allReferenceStability = [initialReferenceStability; ...
        vertcat(postReferenceCells{1:nAdaptiveFits}); ...
        fixedReferenceStability];
    allAudits = [initialAudit; vertcat(postAuditCells{1:nAdaptiveFits})];
    allProfiles = [initialProfile; vertcat(postProfileCells{1:nAdaptiveFits}); ...
        finalProfile];
    adaptiveIterations = vertcat(iterationCells{1:nAdaptiveFits});
    rankSummary = summarize_profile(finalProfile, C.Minit, fixedReferenceM, ...
        size(xo,1));
    referenceChoice = [initialReferenceChoice; ...
        vertcat_reference_choices(postReferenceCells, postProfileCells, ...
        nAdaptiveFits, adaptiveFit, master, lambda, frequencyIndex, C); ...
        fixedReferenceChoice];
    % The helper above reconstructs choice summaries from stored stability
    % rows so every adaptive cycle is represented in one machine-readable file.

    result.schemaVersion = 2;
    result.createdAt = char(datetime('now', ...
        'Format', 'yyyy-MM-dd HH:mm:ss'));
    result.C = C;
    result.outputDir = outputDir;
    result.design = Dsgn;
    result.referenceChoice = referenceChoice;
    result.referenceStability = allReferenceStability;
    result.frequencyAudit = allAudits;
    result.rankProfiles = allProfiles;
    result.finalProfile = finalProfile;
    result.rankSummary = rankSummary;
    result.adaptiveIterations = adaptiveIterations;
    result.parameters = parameters;
    result.fitStatus = fitStatus;
    result.predictionMetrics = predictionMetrics;
    result.externalAudit = externalAudit;
    result.objectiveTimings = objectiveTimings;
    result.stageTimings = stageTimings;

    writetable(referenceChoice, fullfile(outputDir, 'reference_choice.csv'));
    writetable(allReferenceStability, fullfile(outputDir, ...
        'reference_stability.csv'));
    writetable(allAudits, fullfile(outputDir, 'frequency_audit.csv'));
    writetable(allProfiles, fullfile(outputDir, 'rank_profiles.csv'));
    writetable(finalProfile, fullfile(outputDir, 'final_rank_profile.csv'));
    writetable(rankSummary, fullfile(outputDir, 'rank_summary.csv'));
    writetable(adaptiveIterations, fullfile(outputDir, ...
        'adaptive_iterations.csv'));
    writetable(parameters, fullfile(outputDir, 'parameters.csv'));
    writetable(fitStatus, fullfile(outputDir, 'fit_status.csv'));
    writetable(predictionMetrics, fullfile(outputDir, ...
        'prediction_metrics.csv'));
    writetable(externalAudit, fullfile(outputDir, ...
        'external_numerical_audit.csv'));
    writetable(objectiveTimings, fullfile(outputDir, ...
        'objective_timings.csv'));
    writetable(stageTimings, fullfile(outputDir, 'stage_timings.csv'));
    if ~isfield(C,'writeReplicateFigures') || C.writeReplicateFigures
        make_figures(result, figureDir);
    end
    write_summary(result, fullfile(outputDir, 'SUMMARY.md'));
    write_metadata(result, fullfile(outputDir, 'run_metadata.txt'));
    if ~isfield(C,'saveReplicateMat') || C.saveReplicateMat
        save(fullfile(outputDir, 'result.mat'), 'result', '-v7.3');
    end
    fprintf('Algorithm-B experiment complete: %s\n', outputDir);
end

function validate_config(C)
    if ~isfolder(C.baseSourceRoot)
        error('Established source suite is absent: %s', C.baseSourceRoot);
    end
    if any(diff(C.candidateM) <= 0) || any(diff(C.referenceM) <= 0)
        error('Candidate and reference ranks must be strictly increasing.');
    end
    if max(C.candidateM) >= min(C.referenceM)
        error('Every candidate rank must be below every reference rank.');
    end
    if C.Minit > max(C.candidateM)
        error('Initial rank must be represented by the master basis.');
    end
    if C.klTolerance <= 0 || C.referenceTolerance <= 0
        error('KL tolerances must be positive.');
    end
    validAuditRules = {'exact_kl','trace_certificate'};
    if ~isfield(C,'auditRule') || ...
            ~any(strcmpi(C.auditRule,validAuditRules))
        error('auditRule must be exact_kl or trace_certificate.');
    end
    if C.maxAdaptiveIterations < 1 || C.correctionStarts < 1 || C.starts < 1
        error('Iteration and multistart counts must be positive integers.');
    end
end

function validate_nested_bases(domainXy, observedXy, master, ranks)
    for M = ranks(:).'
        candidate = mkbasis(domainXy, master.c, M);
        if max(abs(candidate.radial - master.radial(1:M))) > 1e-12 || ...
                any(candidate.ind ~= master.ind(1:M,:), 'all')
            error('Basis rank M=%d is not nested in the master basis.', M);
        end
        direct = features(observedXy, candidate);
        if max(abs(direct - master.PhiObserved(:,1:M)), [], 'all') > 1e-12
            error('Nested feature validation failed at M=%d.', M);
        end
    end
end

function B = subset_basis(master, M, domainXy)
    B = mkbasis(domainXy, master.c, M);
end

function [rows, choice] = stabilize_reference(theta, master, lambda, ...
        frequencyIndex, C, stage)
    n = size(master.PhiObserved,1);
    endpoint = false(numel(lambda),1);
    endpoint([1,end]) = true;
    rowCells = cell((numel(C.referenceM)-1) * numel(lambda), 1);
    pairMax = nan(numel(C.referenceM)-1,1);
    pairMedian = nan(numel(C.referenceM)-1,1);
    cursor = 0;
    selectedPair = [];
    for q = 1:numel(C.referenceM)-1
        M = C.referenceM(q);
        Mplus = C.referenceM(q+1);
        errors = zeros(numel(lambda),1);
        for k = 1:numel(lambda)
            P = hsgp_covariance(theta, master, M, lambda(k), true);
            Q = hsgp_covariance(theta, master, Mplus, lambda(k), true);
            FP = covariance_factor(P);
            FQ = covariance_factor(Q);
            factor = 1 - 0.5 * endpoint(k);
            forward = covariance_kl_factored(P,FP,Q,FQ,factor) / n;
            reverse = covariance_kl_factored(Q,FQ,P,FP,factor) / n;
            errors(k) = max(forward,reverse);
            cursor = cursor + 1;
            rowCells{cursor} = table(stage, frequencyIndex(k), ...
                lambda(k), endpoint(k), M, Mplus, forward, reverse, ...
                errors(k), 'VariableNames', {'stage','frequencyIndex', ...
                'lambda','endpointFrequency','referenceM', ...
                'largerReferenceM','forwardKlPerSite', ...
                'reverseKlPerSite','symmetricMaximumKlPerSite'});
        end
        pairMax(q) = max(errors);
        pairMedian(q) = median(errors);
        if pairMax(q) <= C.referenceTolerance
            selectedPair = q;
            break;
        end
    end
    rows = vertcat(rowCells{1:cursor});
    stable = ~isempty(selectedPair);
    if isempty(selectedPair)
        selectedPair = numel(C.referenceM)-1;
    end
    choice = table(stage, C.referenceM(selectedPair), ...
        C.referenceM(selectedPair+1), pairMax(selectedPair), ...
        pairMedian(selectedPair), C.referenceTolerance, stable, ...
        'VariableNames', {'stage','selectedReferenceM', ...
        'checkingReferenceM','maximumSymmetricKlPerSite', ...
        'medianSymmetricKlPerSite','referenceTolerance', ...
        'referenceStable'});
end

function [audit, profile] = rank_audit_and_select(theta, master, lambda, ...
        frequencyIndex, referenceM, C, stage, cycle)
    switch lower(C.auditRule)
        case 'exact_kl'
            [audit, profile] = exact_kl_audit_and_select(theta, master, ...
                lambda, frequencyIndex, referenceM, C, stage, cycle);
        case 'trace_certificate'
            [audit, profile] = trace_audit_and_select(theta, master, ...
                lambda, frequencyIndex, referenceM, C, stage, cycle);
        otherwise
            error('Unknown rank-audit rule: %s',C.auditRule);
    end
end

function [audit, profile] = exact_kl_audit_and_select(theta, master, lambda, ...
        frequencyIndex, referenceM, C, stage, cycle)
    nCandidate = numel(C.candidateM);
    nFrequency = numel(lambda);
    n = size(master.PhiObserved,1);
    endpoint = false(nFrequency,1);
    endpoint([1,end]) = true;
    errors = zeros(nCandidate,nFrequency);
    rowCells = cell(nCandidate*nFrequency,1);
    cursor = 0;
    for k = 1:nFrequency
        P = hsgp_covariance(theta, master, referenceM, lambda(k), true);
        FP = covariance_factor(P);
        factor = 1 - 0.5 * endpoint(k);
        for q = 1:nCandidate
            M = C.candidateM(q);
            Q = hsgp_covariance(theta, master, M, lambda(k), true);
            FQ = covariance_factor(Q);
            errors(q,k) = covariance_kl_factored(P,FP,Q,FQ,factor)/n;
            cursor = cursor+1;
            rowCells{cursor} = table(stage, cycle, frequencyIndex(k), ...
                lambda(k), endpoint(k), referenceM, M, errors(q,k), ...
                'VariableNames', {'stage','cycle','frequencyIndex', ...
                'lambda','endpointFrequency','referenceM','candidateM', ...
                'klPerFittingSite'});
        end
    end
    audit = vertcat(rowCells{1:cursor});
    selectedM = zeros(nFrequency,1);
    selectedError = zeros(nFrequency,1);
    tailError = zeros(nFrequency,1);
    thresholdMet = false(nFrequency,1);
    for k = 1:nFrequency
        tailMaximum = flipud(cummax(flipud(errors(:,k))));
        position = find(tailMaximum <= C.klTolerance,1,'first');
        thresholdMet(k) = ~isempty(position);
        if isempty(position)
            position = nCandidate;
        end
        selectedM(k) = C.candidateM(position);
        selectedError(k) = errors(position,k);
        tailError(k) = tailMaximum(position);
    end
    profile = table(repmat(stage,nFrequency,1), ...
        repmat(cycle,nFrequency,1), frequencyIndex, lambda(:), ...
        repmat(referenceM,nFrequency,1), selectedM, thresholdMet, ...
        selectedError, tailError, ...
        'VariableNames', {'stage','cycle','frequencyIndex','lambda', ...
        'referenceM','selectedM','thresholdMet','selectedKlPerSite', ...
        'persistentTailMaximumKlPerSite'});
end

function [audit, profile] = trace_audit_and_select(theta, master, lambda, ...
        frequencyIndex, referenceM, C, stage, cycle)
%TRACE_AUDIT_AND_SELECT Apply the covariance-scaled trace certificate.
%
% For candidate covariance Q and positive-semidefinite omitted tail R,
% t=tr(Q^{-1}R) is evaluated with the Cholesky factor of Q. The selection
% score c_k/n * {t-log(1+t)} is an upper bound on the normalized Gaussian
% KL discrepancy. The internal reference rank is appended as a guaranteed
% zero-tail fallback but is selected only when no smaller candidate passes.

    candidates = unique([C.candidateM(:).', referenceM], 'sorted');
    candidates = candidates(candidates <= referenceM);
    nCandidate = numel(candidates);
    nFrequency = numel(lambda);
    n = size(master.PhiObserved,1);
    endpoint = false(nFrequency,1);
    endpoint([1,end]) = true;
    bounds = zeros(nCandidate,nFrequency);
    rowCells = cell(nCandidate*nFrequency,1);
    cursor = 0;
    for k = 1:nFrequency
        factor = 1 - 0.5 * endpoint(k);
        allWeights = mult(theta,master.radial(1:referenceM),lambda(k));
        for q = 1:nCandidate
            M = candidates(q);
            Q = hsgp_covariance(theta, master, M, lambda(k), true);
            FQ = covariance_factor(Q);
            if M < referenceM
                PhiTail = master.PhiObserved(:,M+1:referenceM);
                solvedTail = FQ.L\PhiTail;
                tailWeights = allWeights(M+1:referenceM);
                traceAction = real(sum(tailWeights(:).' .* ...
                    sum(abs(solvedTail).^2,1)));
            else
                traceAction = 0;
            end
            traceAction = max(traceAction,0);
            bounds(q,k) = factor * ...
                (traceAction-log1p(traceAction)) / n;
            cursor = cursor+1;
            rowCells{cursor} = table(stage, cycle, frequencyIndex(k), ...
                lambda(k), endpoint(k), referenceM, M, bounds(q,k), ...
                'VariableNames', {'stage','cycle','frequencyIndex', ...
                'lambda','endpointFrequency','referenceM','candidateM', ...
                'klPerFittingSite'});
        end
    end
    audit = vertcat(rowCells{1:cursor});
    selectedM = zeros(nFrequency,1);
    selectedBound = zeros(nFrequency,1);
    tailBound = zeros(nFrequency,1);
    thresholdMet = false(nFrequency,1);
    for k = 1:nFrequency
        tailMaximum = flipud(cummax(flipud(bounds(:,k))));
        position = find(tailMaximum <= C.klTolerance,1,'first');
        thresholdMet(k) = ~isempty(position);
        if isempty(position)
            position = nCandidate;
        end
        selectedM(k) = candidates(position);
        selectedBound(k) = bounds(position,k);
        tailBound(k) = tailMaximum(position);
    end
    profile = table(repmat(stage,nFrequency,1), ...
        repmat(cycle,nFrequency,1), frequencyIndex, lambda(:), ...
        repmat(referenceM,nFrequency,1), selectedM, thresholdMet, ...
        selectedBound, tailBound, ...
        'VariableNames', {'stage','cycle','frequencyIndex','lambda', ...
        'referenceM','selectedM','thresholdMet','selectedKlPerSite', ...
        'persistentTailMaximumKlPerSite'});
end

function profile = profile_at_used_ranks(audit, usedM, frequencyIndex, ...
        lambda, C, stage, cycle)
    nFrequency = numel(lambda);
    selectedError = zeros(nFrequency,1);
    tailError = zeros(nFrequency,1);
    thresholdMet = false(nFrequency,1);
    referenceM = zeros(nFrequency,1);
    for k = 1:nFrequency
        subset = audit(audit.frequencyIndex == frequencyIndex(k),:);
        [candidateM,order] = sort(subset.candidateM);
        errors = subset.klPerFittingSite(order);
        position = find(candidateM == usedM(k),1);
        if isempty(position)
            error('Used rank M=%d is absent from the audit.', usedM(k));
        end
        tails = flipud(cummax(flipud(errors)));
        selectedError(k) = errors(position);
        tailError(k) = tails(position);
        thresholdMet(k) = tails(position) <= C.klTolerance;
        referenceM(k) = subset.referenceM(1);
    end
    profile = table(repmat(stage,nFrequency,1), ...
        repmat(cycle,nFrequency,1), frequencyIndex, lambda(:), ...
        referenceM, usedM(:), thresholdMet, selectedError, tailError, ...
        'VariableNames', {'stage','cycle','frequencyIndex','lambda', ...
        'referenceM','selectedM','thresholdMet','selectedKlPerSite', ...
        'persistentTailMaximumKlPerSite'});
end

function require_stable_reference(choice, location)
    if ~choice.referenceStable
        error(['High-rank reference failed the prespecified stability ', ...
            'criterion at %s.'], location);
    end
end

function require_resolved_profile(profile, location)
    if any(~profile.thresholdMet)
        bad = profile.frequencyIndex(~profile.thresholdMet);
        error('Candidate ranks unresolved at %s; frequencies: %s', ...
            location, mat2str(bad(:).'));
    end
end

function K = hsgp_covariance(theta, master, M, lambda, sameSites)
    P = master.PhiObserved(:,1:M);
    w = mult(theta,master.radial(1:M),lambda);
    K = (P.*w')*P';
    if sameSites
        K = real((K+K')/2)+theta(5)*eye(size(P,1));
    end
end

function value = covariance_kl(P,Q,factor)
    P = real((P+P')/2);
    Q = real((Q+Q')/2);
    FP = covariance_factor(P);
    FQ = covariance_factor(Q);
    value = covariance_kl_factored(P,FP,Q,FQ,factor);
end

function F = covariance_factor(K)
    K = real((K+K')/2);
    [F.L,~] = cholj(K);
    F.logdet = 2*sum(log(diag(F.L)));
end

function value = covariance_kl_factored(P,FP,Q,FQ,factor)
    n = size(P,1);
    traceTerm = trace(FQ.L'\(FQ.L\P));
    value = factor*real(traceTerm-n+FQ.logdet-FP.logdet);
    if value < 0 && value > -1e-8
        value = 0;
    elseif value < -1e-8
        error('Materially negative covariance KL divergence.');
    end
end

function fit = fit_variable_rank(Dpositive, initialTheta, profileM, ...
        master, lambda, C, seed, numberStarts)
    base = to_z(initialTheta,C);
    rng(seed,'twister');
    starts = zeros(numel(base),numberStarts);
    starts(:,1) = base;
    for q = 2:numberStarts
        starts(:,q) = base+C.opt.startScale*randn(size(base));
    end
    opts = optimset('Display','off','TolX',C.opt.tolX, ...
        'TolFun',C.opt.tolFun,'MaxIter',C.opt.maxIter, ...
        'MaxFunEvals',C.opt.maxEval);
    variableCache = make_variable_cache(profileM,master);
    objective = @(z) variable_nll_z(z,Dpositive,profileM,master,lambda,C, ...
        variableCache);
    solution = zeros(size(starts));
    value = inf(numberStarts,1);
    flag = zeros(numberStarts,1);
    output = cell(numberStarts,1);
    timer = tic;
    for q = 1:numberStarts
        [solution(:,q),value(q),flag(q),output{q}] = ...
            fminsearch(objective,starts(:,q),opts);
    end
    [fit.nll,best] = min(value);
    fit.z = solution(:,best);
    fit.theta = from_z(fit.z,C);
    fit.exitFlag = flag(best);
    fit.output = output;
    fit.seconds = toc(timer);
    fit.starts = starts;
    fit.numberStarts = numberStarts;
    fit.solutions = solution;
    fit.values = value;
    fit.flags = flag;
    fit.engine = 'adaptive_hsgp';
    fit.profileM = profileM;
    fit.boundaryDistance = min((fit.theta-C.lower)./(C.upper-C.lower), ...
        1-(fit.theta-C.lower)./(C.upper-C.lower));
end

function value = variable_nll_z(z,Dpositive,profileM,master,lambda,C,cache)
    theta = from_z(z,C);
    value = variable_nll_theta(theta,Dpositive,profileM,master,lambda,C,cache);
end

function value = variable_nll_theta(theta,Dpositive,profileM,master,lambda,C,cache)
    n = size(Dpositive,2);
    K = numel(lambda);
    value = 0;
    if nargin < 7 || isempty(cache)
        cache = make_variable_cache(profileM,master);
    end
    try
        for k = 1:K
            entry = cache.entries(cache.frequencyEntry(k));
            w = mult(theta,entry.radial,lambda(k));
            F = lrfactor(entry.P,w,theta(5),entry.G);
            if k == 1 || k == K
                d = real(Dpositive(k,:)).';
                value = value+0.5*(F.logdet+real(d'*lrsolve(F,d))+ ...
                    n*log(2*pi));
            else
                d = Dpositive(k,:).';
                value = value+F.logdet+real(d'*lrsolve(F,d))+n*log(pi);
            end
        end
        if ~isfinite(value)
            error('Nonfinite variable-rank objective.');
        end
    catch
        value = realmax/1e100;
    end
end

function cache = make_variable_cache(profileM,master)
    uniqueM = unique(profileM(:),'stable');
    entries = repmat(struct('M',0,'P',[],'G',[],'radial',[]), ...
        numel(uniqueM),1);
    for q = 1:numel(uniqueM)
        M = uniqueM(q);
        P = master.PhiObserved(:,1:M);
        entries(q).M = M;
        entries(q).P = P;
        entries(q).G = P'*P;
        entries(q).radial = master.radial(1:M);
    end
    [found,map] = ismember(profileM(:),uniqueM);
    if ~all(found)
        error('Could not map all frequency-specific ranks into the cache.');
    end
    cache.entries = entries;
    cache.frequencyEntry = map;
end

function prediction = predict_positive(Yo,xo,xt,theta,engine,profileM, ...
        master,lambda,R,C)
    T = size(Yo,1);
    K = numel(lambda);
    Do = fft(Yo,[],1)/sqrt(T);
    Dmean = complex(zeros(T,size(xt,1)));
    if engine == "hsgp"
        hsgpCache = make_variable_cache(profileM,master);
    end
    for k = 1:K
        if engine == "exact"
            Koo = cov_ref(theta,xo,xo,lambda(k),R,true);
            Kto = cov_ref(theta,xt,xo,lambda(k),R,false);
            [L,~] = cholj(Koo);
            solved = L'\(L\Do(k,:).');
        else
            entry = hsgpCache.entries(hsgpCache.frequencyEntry(k));
            M = entry.M;
            Po = entry.P;
            Pt = master.PhiTarget(:,1:M);
            w = mult(theta,entry.radial,lambda(k));
            factor = lrfactor(Po,w,theta(5),entry.G);
            Kto = (Pt.*w')*Po';
            solved = lrsolve(factor,Do(k,:).');
        end
        Dmean(k,:) = (Kto*solved).';
    end
    if mod(T,2)==0
        interior = 2:K-1;
    else
        interior = 2:K;
    end
    for k = interior
        Dmean(T-k+2,:) = conj(Dmean(k,:));
    end
    prediction = real(ifft(Dmean,[],1)*sqrt(T));
end

function audit = external_numerical_audit(theta,master,xo,lambda, ...
        frequencyIndex,referenceM,adaptiveM,initialM,R,C)
    methods = ["high_rank_internal_reference"; ...
        "initialization_rank"; "algorithm_B_adaptive_rank"];
    n = size(xo,1);
    rows = cell(numel(methods)*numel(lambda),1);
    cursor = 0;
    for k = 1:numel(lambda)
        P = cov_ref(theta,xo,xo,lambda(k),R,true);
        factor = 1-0.5*((k==1)||(k==numel(lambda)));
        ranks = [referenceM,initialM,adaptiveM(k)];
        for q = 1:numel(methods)
            Q = hsgp_covariance(theta,master,ranks(q),lambda(k),true);
            cursor = cursor+1;
            rows{cursor} = table(methods(q),frequencyIndex(k),lambda(k), ...
                ranks(q),covariance_kl(P,Q,factor)/n, ...
                norm(P-Q,'fro')/norm(P,'fro'), ...
                'VariableNames', {'method','frequencyIndex','lambda','M', ...
                'klPerSiteVersusNumericalReference','relativeFrobeniusError'});
        end
    end
    audit = vertcat(rows{:});
end

function row = parameter_row(method,theta,objective,seconds,evaluations)
    row = table(method,objective,seconds,evaluations,theta(1),theta(2), ...
        theta(3),theta(4),theta(5), ...
        'VariableNames', {'method','objective','fitSeconds', ...
        'objectiveEvaluations','vf','rho','kappa0','gamma','vn'});
end

function count = fit_eval_count(fit)
    count = NaN;
    if ~isfield(fit,'output')
        return;
    end
    out = fit.output;
    if isstruct(out)
        out = {out};
    end
    values = zeros(numel(out),1);
    for q = 1:numel(out)
        if isfield(out{q},'funcCount')
            values(q) = out{q}.funcCount;
        else
            values(q) = NaN;
        end
    end
    count = sum(values,'omitnan');
end

function summary = summarize_profile(profile,initialM,referenceM,n)
    M = profile.selectedM;
    summary = table(numel(M),n,min(M),median(M),mean(M),max(M), ...
        sum(M<n),sum(M>=n),sum(double(M).^3), ...
        numel(M)*double(initialM)^3,numel(M)*double(referenceM)^3, ...
        sum(double(M).^3)/(numel(M)*double(initialM)^3), ...
        sum(double(M).^3)/(numel(M)*double(referenceM)^3), ...
        max(profile.persistentTailMaximumKlPerSite), ...
        'VariableNames', {'numberFrequencies','numberFittingSites', ...
        'minimumM','medianM','meanM','maximumM', ...
        'numberFeatureSpaceFrequencies','numberObservationSpaceFrequencies', ...
        'adaptiveCubicProxy','initializationCubicProxy', ...
        'fixedReferenceCubicProxy','adaptiveToInitializationCubicProxyRatio', ...
        'adaptiveToFixedReferenceCubicProxyRatio', ...
        'maximumPersistentKlPerSite'});
end

function rows = benchmark_objectives(Y,Dpositive,xo,thetaInitial, ...
        thetaReference,thetaAdaptive,thetaNumerical,initialBasis, ...
        referenceBasis,profileM,master,lambda,R,C)
    names = ["initialization_fit";"fixed_reference_rank_fit"; ...
        "algorithm_B_adaptive_fit";"numerical_inversion_fit"];
    initialCache.P = features(xo,initialBasis);
    initialCache.G = initialCache.P'*initialCache.P;
    referenceCache.P = features(xo,referenceBasis);
    referenceCache.G = referenceCache.P'*referenceCache.P;
    adaptiveCache = make_variable_cache(profileM,master);
    functions = {@() nll_model(to_z(thetaInitial,C),Y,xo,C,R, ...
        'hsgp',initialBasis,initialCache), ...
        @() nll_model(to_z(thetaReference,C),Y,xo,C,R, ...
        'hsgp',referenceBasis,referenceCache), ...
        @() variable_nll_theta(thetaAdaptive,Dpositive,profileM, ...
        master,lambda,C,adaptiveCache), ...
        @() nll_model(to_z(thetaNumerical,C),Y,xo,C,R, ...
        'exact',[],[])};
    rowCells = cell(numel(names),1);
    for q = 1:numel(names)
        objective = functions{q}();
        seconds = zeros(C.nObjectiveTimingReplicates,1);
        for r = 1:C.nObjectiveTimingReplicates
            timer = tic;
            functions{q}();
            seconds(r) = toc(timer);
        end
        rowCells{q} = table(names(q),objective, ...
            C.nObjectiveTimingReplicates,mean(seconds),median(seconds), ...
            min(seconds),max(seconds), ...
            'VariableNames', {'method','objective','nReplicates', ...
            'meanSeconds','medianSeconds','minimumSeconds','maximumSeconds'});
    end
    rows = vertcat(rowCells{:});
end

function choices = vertcat_reference_choices(referenceCells,profileCells, ...
        nCycles,adaptiveFit,master,lambda,frequencyIndex,C)
    choices = cell(nCycles,1);
    theta = adaptiveFit.theta;
    for cycle = 1:nCycles
        stage = "adaptive_refit_cycle_" + cycle;
        rows = referenceCells{cycle};
        pair = unique(rows(:,{'referenceM','largerReferenceM'}),'rows','stable');
        maxima = zeros(height(pair),1);
        medians = zeros(height(pair),1);
        for q = 1:height(pair)
            subset = rows(rows.referenceM==pair.referenceM(q) & ...
                rows.largerReferenceM==pair.largerReferenceM(q),:);
            maxima(q) = max(subset.symmetricMaximumKlPerSite);
            medians(q) = median(subset.symmetricMaximumKlPerSite);
        end
        selected = find(maxima<=C.referenceTolerance,1,'first');
        stable = ~isempty(selected);
        if isempty(selected), selected=height(pair); end
        choices{cycle} = table(stage,pair.referenceM(selected), ...
            pair.largerReferenceM(selected),maxima(selected), ...
            medians(selected),C.referenceTolerance,stable, ...
            'VariableNames', {'stage','selectedReferenceM', ...
            'checkingReferenceM','maximumSymmetricKlPerSite', ...
            'medianSymmetricKlPerSite','referenceTolerance', ...
            'referenceStable'});
    end
    choices = vertcat(choices{:});
    %#ok<NASGU> Inputs retained in the signature to make the reconstruction
    % of choice summaries explicit and future-proof if the stored rows change.
end

function make_figures(result,figureDir)
    C = result.C;
    profile = result.finalProfile;
    figure('Visible',C.figureVisible,'Color','w', ...
        'Position',[100,100,1050,430]);
    stairs(profile.frequencyIndex,profile.selectedM,'LineWidth',1.8, ...
        'Color',[0.00,0.45,0.74]);
    hold on;
    fixedReferenceM = result.referenceChoice.selectedReferenceM(1);
    yline(fixedReferenceM,'--','fixed M_{ref}', ...
        'Color',[0.85,0.33,0.10],'LineWidth',1.4);
    ylim([0,1.08*fixedReferenceM]);
    xlabel('DFT frequency index k');
    ylabel('retained rank M_k');
    title('Algorithm-B fitted and verified rank profile');
    grid on;
    save_figure(gcf,figureDir,'algorithm_b_rank_profile');

    figure('Visible',C.figureVisible,'Color','w', ...
        'Position',[100,100,1050,430]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile;
    bar(result.predictionMetrics.latentRmse,'FaceColor',[0.00,0.45,0.74]);
    xticks(1:height(result.predictionMetrics));
    xticklabels(pretty_method(result.predictionMetrics.method));
    xtickangle(20);
    ylabel('held-site latent RMSE');
    grid on;
    nexttile;
    times = [result.stageTimings.numericalReferenceTotalSeconds; ...
        result.stageTimings.screenedFixedTotalSeconds; ...
        result.stageTimings.fixedReferenceTotalSeconds; ...
        result.stageTimings.adaptiveComparableTotalSeconds];
    bar(times,'FaceColor',[0.85,0.33,0.10]);
    xticks(1:4);
    xticklabels({'numerical inversion','fixed M=121', ...
        'fixed M_{ref}','adaptive M_k'});
    xtickangle(20);
    ylabel('wall time (s)');
    grid on;
    save_figure(gcf,figureDir,'algorithm_b_accuracy_and_time');
end

function labels = pretty_method(methods)
    labels = cell(numel(methods),1);
    for q = 1:numel(methods)
        labels{q} = strrep(char(methods(q)),'_',' ');
    end
end

function save_figure(fig,directory,name)
    exportgraphics(fig,fullfile(directory,[name,'.pdf']), ...
        'ContentType','vector');
    exportgraphics(fig,fullfile(directory,[name,'.png']),'Resolution',300);
    savefig(fig,fullfile(directory,[name,'.fig']));
    close(fig);
end

function write_summary(result,filename)
    C = result.C;
    p = result.predictionMetrics;
    r = result.rankSummary;
    t = result.stageTimings;
    e = result.externalAudit;
    internal = e(e.method=="high_rank_internal_reference",:);
    adaptive = e(e.method=="algorithm_B_adaptive_rank",:);
    fid = fopen(filename,'w');
    if fid<0, error('Could not create summary file.'); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'# Algorithm-B %s result\n\n',C.profile);
    fprintf(fid,'- Fitting sites: %d; held sites: %d; T: %d.\n', ...
        r.numberFittingSites,numel(result.design.test),C.T);
    fprintf(fid,'- KL tolerance: %.4g; reference tolerance: %.4g.\n', ...
        C.klTolerance,C.referenceTolerance);
    fprintf(fid,'- Final ranks: min %.0f, median %.0f, mean %.2f, max %.0f.\n', ...
        r.minimumM,r.medianM,r.meanM,r.maximumM);
    fprintf(fid,'- Fixed reference rank: %.0f.\n', ...
        result.referenceChoice.selectedReferenceM(1));
    fprintf(fid,'- Adaptive/fixed-reference cubic proxy: %.4f.\n', ...
        r.adaptiveToFixedReferenceCubicProxyRatio);
    fprintf(fid,'- Maximum persistent KL per fitting site: %.6g.\n', ...
        r.maximumPersistentKlPerSite);
    fprintf(fid,'- Algorithm-B first-use wall time: %.3f s.\n', ...
        t.algorithmBFirstUseSeconds);
    fprintf(fid,['- Adaptive matched-boundary wall time (excluding shared ', ...
        'basis setup): %.3f s.\n'],t.adaptiveComparableTotalSeconds);
    fprintf(fid,'- Fixed screened-rank fit-plus-prediction wall time: %.3f s.\n', ...
        t.screenedFixedTotalSeconds);
    fprintf(fid,'- Fixed-reference fit-plus-prediction wall time: %.3f s.\n', ...
        t.fixedReferenceTotalSeconds);
    fprintf(fid,'- Fixed-reference validation wall time: %.3f s.\n', ...
        t.fixedReferenceValidationSeconds);
    fprintf(fid,'- Numerical-inversion total wall time: %.3f s.\n\n', ...
        t.numericalReferenceTotalSeconds);
    fprintf(fid,'## Held-site reconstruction\n\n');
    for q=1:height(p)
        fprintf(fid,'- %s: latent RMSE %.6f; noisy RMSE %.6f.\n', ...
            p.method(q),p.latentRmse(q),p.noisyRmse(q));
    end
    fprintf(fid,'\n## External numerical-reference audit\n\n');
    fprintf(fid,['- Internal high-rank reference: maximum KL/site %.6g; ', ...
        'maximum relative Frobenius error %.6g.\n'], ...
        max(internal.klPerSiteVersusNumericalReference), ...
        max(internal.relativeFrobeniusError));
    fprintf(fid,['- Final adaptive rank: maximum KL/site %.6g; ', ...
        'maximum relative Frobenius error %.6g.\n'], ...
        max(adaptive.klPerSiteVersusNumericalReference), ...
        max(adaptive.relativeFrobeniusError));
    fprintf(fid,['\nThe numerical-inversion covariance was not used to ', ...
        'select the reference rank or the frequency-specific ranks.\n']);
end

function write_metadata(result,filename)
    C = result.C;
    fid = fopen(filename,'w');
    if fid<0, error('Could not create metadata file.'); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'suite_version=%s\n',C.version);
    fprintf(fid,'profile=%s\n',C.profile);
    fprintf(fid,'created_at=%s\n',result.createdAt);
    fprintf(fid,'simulation_seed=%d\n',C.simulationSeed);
    fprintf(fid,'optimizer_seed=%d\n',C.optimizerSeed);
    fprintf(fid,'boundary_factor=%.17g\n',C.c);
    fprintf(fid,'initial_M=%d\n',C.Minit);
    fprintf(fid,'candidate_M=%s\n',mat2str(C.candidateM));
    fprintf(fid,'reference_M=%s\n',mat2str(C.referenceM));
    fprintf(fid,'kl_tolerance=%.17g\n',C.klTolerance);
    fprintf(fid,'reference_tolerance=%.17g\n',C.referenceTolerance);
    fprintf(fid,'audit_rule=%s\n',C.auditRule);
    fprintf(fid,['rank_selection_data=fitting coordinates, initial fitted ', ...
        'parameters, and high-rank HSGP covariance only\n']);
    fprintf(fid,['held_out_responses_used_for=final reconstruction ', ...
        'evaluation only\n']);
    fprintf(fid,['numerical_inversion_used_for=simulation and external ', ...
        'evaluation only\n']);
end

function ensure_directory(pathValue)
    if ~isfolder(pathValue), mkdir(pathValue); end
end

