function aggregateDir = run_matched_20(workers,auditRule,klTolerance)
%RUN_MATCHED_20 Run the four-method comparison on the trusted 20 records.
%
% The simulation seeds, held-site manifest, optimizer seeds, number of
% multistarts, and optimization budget match the established Table-4 HPC
% experiment. Reported wall time includes fitting, method-specific rank
% selection when applicable, and held-site reconstruction. Common data
% generation, quadrature-table construction, and basis-bank construction are
% recorded but excluded from the four compared totals.

    if nargin < 1 || isempty(workers)
        workers = 5;
    end
    if nargin < 2 || isempty(auditRule)
        auditRule = 'exact_kl';
    end
    if nargin < 3 || isempty(klTolerance)
        klTolerance = 0.01;
    end
    workers = max(1,floor(workers));
    sourceRoot = fileparts(mfilename('fullpath'));
    experimentRoot = fileparts(sourceRoot);
    addpath(sourceRoot);
    C0 = matched_config('confirm');
    C0.auditRule = char(auditRule);
    C0.klTolerance = klTolerance;
    C0.writeReplicateFigures = false;
    C0.figureVisible = 'off';
    nReplicates = 20;

    stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    ruleTag = regexprep(lower(C0.auditRule),'[^a-z0-9]+','_');
    aggregateDir = fullfile(experimentRoot,'results', ...
        ['matched_20_',ruleTag,'_',stamp]);
    rawRoot = fullfile(aggregateDir,'raw');
    ensure_directory_local(rawRoot);

    if workers > 1
        pool = gcp('nocreate');
        if ~isempty(pool) && pool.NumWorkers ~= workers
            delete(pool);
            pool = [];
        end
        if isempty(pool)
            parpool('local',workers);
        end
    end

    metricCells = cell(nReplicates,1);
    rankCells = cell(nReplicates,1);
    outputDirs = strings(nReplicates,1);
    if workers > 1
        parfor r = 1:nReplicates
            [metricCells{r},rankCells{r},outputDirs(r)] = ...
                run_one(r,C0,rawRoot);
        end
    else
        for r = 1:nReplicates
            [metricCells{r},rankCells{r},outputDirs(r)] = ...
                run_one(r,C0,rawRoot);
        end
    end

    replicateMetrics = vertcat(metricCells{:});
    rankProfiles = vertcat(rankCells{:});
    writetable(replicateMetrics,fullfile(aggregateDir, ...
        'replicate_metrics.csv'));
    writetable(rankProfiles,fullfile(aggregateDir, ...
        'adaptive_rank_profiles.csv'));

    methodSummary = summarize_methods(replicateMetrics,C0.seed);
    pairedSummary = summarize_pairs(replicateMetrics,C0.seed);
    timingSummary = summarize_timing_contrasts(replicateMetrics,C0.seed);
    writetable(methodSummary,fullfile(aggregateDir,'four_method_summary.csv'));
    writetable(pairedSummary,fullfile(aggregateDir, ...
        'paired_rmse_differences.csv'));
    writetable(timingSummary,fullfile(aggregateDir, ...
        'paired_timing_contrasts.csv'));

    metadata = table((1:nReplicates).', ...
        C0.seed+10000+(1:nReplicates).', ...
        C0.seed+20000+(1:nReplicates).',outputDirs, ...
        'VariableNames',{'replicate','simulationSeed','optimizerSeed', ...
        'replicateOutputDirectory'});
    writetable(metadata,fullfile(aggregateDir,'run_manifest.csv'));
    write_readme(aggregateDir,C0,methodSummary,pairedSummary,timingSummary);
    save(fullfile(aggregateDir,'aggregate_result.mat'), ...
        'replicateMetrics','rankProfiles','methodSummary','pairedSummary', ...
        'timingSummary','metadata','C0','-v7.3');
    fprintf('Matched 20-replicate experiment complete: %s\n',aggregateDir);
end

function [metrics,ranks,outputDir] = run_one(r,C0,rawRoot)
    C = C0;
    C.simulationSeed = C.seed + 10000 + r;
    C.optimizerSeed = C.seed + 20000 + r;
    C.outputRoot = fullfile(rawRoot,sprintf('replicate_%02d',r));
    ensure_directory_local(C.outputRoot);
    fprintf('[%s] replicate %02d/20 started (audit=%s, epsilon=%.4g)\n', ...
        char(datetime('now','Format','HH:mm:ss')),r,C.auditRule,C.klTolerance);
    result = run_matched_replicate(C);
    fprintf('[%s] replicate %02d/20 completed\n', ...
        char(datetime('now','Format','HH:mm:ss')),r);
    outputDir = string(result.outputDir);

    sourceMethods = ["numerical_inversion_fit"; ...
        "fixed_screened_rank_fit"; "fixed_reference_rank_fit"; ...
        "algorithm_B_adaptive_fit"];
    reportMethods = ["numerical_inversion_reference"; ...
        "fixed_hsgp_M121"; "fixed_high_rank_hsgp"; ...
        "frequency_adaptive_hsgp"];
    p = result.predictionMetrics;
    s = result.fitStatus;
    latent = zeros(4,1);
    observed = zeros(4,1);
    predictionSeconds = zeros(4,1);
    exitFlag = zeros(4,1);
    for q = 1:4
        row = p(p.method==sourceMethods(q),:);
        status = s(s.method==sourceMethods(q),:);
        latent(q) = row.latentRmse;
        observed(q) = row.noisyRmse;
        predictionSeconds(q) = row.predictionSeconds;
        exitFlag(q) = status.exitFlag;
    end
    t = result.stageTimings;
    totalSeconds = [t.numericalReferenceTotalSeconds; ...
        t.screenedFixedTotalSeconds; t.fixedReferenceTotalSeconds; ...
        t.adaptiveComparableTotalSeconds];
    selectedReferenceM = result.referenceChoice.selectedReferenceM(end);
    metrics = table(repmat(r,4,1),reportMethods,latent,observed, ...
        totalSeconds,predictionSeconds,exitFlag,exitFlag>0, ...
        repmat(selectedReferenceM,4,1), ...
        repmat(C.simulationSeed,4,1),repmat(C.optimizerSeed,4,1), ...
        'VariableNames',{'replicate','method','latentRmse','observedRmse', ...
        'totalSeconds','predictionSeconds','exitFlag','converged', ...
        'selectedReferenceM','simulationSeed','optimizerSeed'});

    f = result.finalProfile;
    ranks = table(repmat(r,height(f),1),f.frequencyIndex,f.lambda, ...
        f.selectedM,f.selectedKlPerSite, ...
        f.persistentTailMaximumKlPerSite, ...
        repmat(selectedReferenceM,height(f),1), ...
        'VariableNames',{'replicate','frequencyIndex','lambda', ...
        'selectedM','klPerSite','persistentTailMaximumKlPerSite', ...
        'selectedReferenceM'});
end

function summary = summarize_methods(metrics,seed)
    order = ["numerical_inversion_reference";"fixed_hsgp_M121"; ...
        "fixed_high_rank_hsgp";"frequency_adaptive_hsgp"];
    cells = cell(numel(order),1);
    for q = 1:numel(order)
        z = metrics(metrics.method==order(q),:);
        [ll,lu] = bootstrap_mean(z.latentRmse,seed+100+q,4000);
        [ol,ou] = bootstrap_mean(z.observedRmse,seed+200+q,4000);
        [tl,tu] = bootstrap_mean(z.totalSeconds,seed+300+q,4000);
        cells{q} = table(order(q),height(z),mean(z.latentRmse),ll,lu, ...
            (lu-ll)/2,mean(z.observedRmse),ol,ou,(ou-ol)/2, ...
            mean(z.totalSeconds),tl,tu,(tu-tl)/2, ...
            mean(z.converged),median(z.selectedReferenceM), ...
            'VariableNames',{'method','nReplicates','meanLatentRmse', ...
            'latentCiLower','latentCiUpper','latentPlusMinusHalfWidth', ...
            'meanObservedRmse','observedCiLower','observedCiUpper', ...
            'observedPlusMinusHalfWidth','meanTotalSeconds', ...
            'timeCiLower','timeCiUpper','timePlusMinusHalfWidth', ...
            'convergenceRate','medianSelectedReferenceM'});
    end
    summary = vertcat(cells{:});
end

function summary = summarize_pairs(metrics,seed)
    comparisons = ["fixed_hsgp_M121","numerical_inversion_reference"; ...
        "fixed_high_rank_hsgp","numerical_inversion_reference"; ...
        "frequency_adaptive_hsgp","numerical_inversion_reference"; ...
        "frequency_adaptive_hsgp","fixed_hsgp_M121"; ...
        "frequency_adaptive_hsgp","fixed_high_rank_hsgp"];
    cells = cell(size(comparisons,1),1);
    for q = 1:size(comparisons,1)
        a = sortrows(metrics(metrics.method==comparisons(q,1),:),'replicate');
        b = sortrows(metrics(metrics.method==comparisons(q,2),:),'replicate');
        dl = a.latentRmse-b.latentRmse;
        do = a.observedRmse-b.observedRmse;
        [ll,lu] = bootstrap_mean(dl,seed+400+q,4000);
        [ol,ou] = bootstrap_mean(do,seed+500+q,4000);
        cells{q} = table(comparisons(q,1),comparisons(q,2),mean(dl),ll,lu, ...
            mean(do),ol,ou, ...
            'VariableNames',{'method','comparator','meanLatentDifference', ...
            'latentDifferenceCiLower','latentDifferenceCiUpper', ...
            'meanObservedDifference','observedDifferenceCiLower', ...
            'observedDifferenceCiUpper'});
    end
    summary = vertcat(cells{:});
end

function summary = summarize_timing_contrasts(metrics,seed)
    adaptive = sortrows(metrics(metrics.method== ...
        "frequency_adaptive_hsgp",:),'replicate');
    comparators = ["numerical_inversion_reference"; ...
        "fixed_hsgp_M121";"fixed_high_rank_hsgp"];
    cells = cell(numel(comparators),1);
    for q = 1:numel(comparators)
        b = sortrows(metrics(metrics.method==comparators(q),:),'replicate');
        saved = b.totalSeconds-adaptive.totalSeconds;
        fraction = 1-adaptive.totalSeconds./b.totalSeconds;
        ratio = b.totalSeconds./adaptive.totalSeconds;
        [sl,su] = bootstrap_mean(saved,seed+600+q,4000);
        [fl,fu] = bootstrap_mean(fraction,seed+700+q,4000);
        [rl,ru] = bootstrap_mean(ratio,seed+800+q,4000);
        cells{q} = table("frequency_adaptive_hsgp",comparators(q), ...
            mean(saved),sl,su,mean(fraction),fl,fu,mean(ratio),rl,ru, ...
            'VariableNames',{'method','comparator','meanSecondsSaved', ...
            'secondsSavedCiLower','secondsSavedCiUpper', ...
            'meanFractionSaved','fractionSavedCiLower', ...
            'fractionSavedCiUpper','meanComparatorToAdaptiveRatio', ...
            'ratioCiLower','ratioCiUpper'});
    end
    summary = vertcat(cells{:});
end

function [lo,hi] = bootstrap_mean(x,seed,B)
    x = x(:);
    rng(seed,'twister');
    n = numel(x);
    values = zeros(B,1);
    for b = 1:B
        values(b) = mean(x(randi(n,n,1)));
    end
    bounds = prctile(values,[2.5,97.5]);
    lo = bounds(1);
    hi = bounds(2);
end

function write_readme(directory,C,methods,pairs,timing)
    fid = fopen(fullfile(directory,'SUMMARY.md'),'w');
    if fid < 0, error('Could not create aggregate summary.'); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'# Matched four-method experiment\n\n');
    fprintf(fid,'Rank-audit rule: `%s`.\n\n',C.auditRule);
    fprintf(fid,['Twenty paired realizations use simulation seeds ', ...
        '`seed+10000+r`, optimizer seeds `seed+20000+r`, the trusted ', ...
        'held-site manifest, four multistarts, and the confirm-profile ', ...
        'optimization budget.\n\n']);
    fprintf(fid,['Compared wall time includes fitting, any method-specific ', ...
        'rank selection, and reconstruction. Data simulation, the shared ', ...
        'numerical-correlation table, and shared basis-bank construction ', ...
        'are excluded.\n\n']);
    fprintf(fid,'KL tolerance: %.4g. Initial screened rank: %d.\n\n', ...
        C.klTolerance,C.Minit);
    fprintf(fid,'## Four-method summary\n\n');
    for q = 1:height(methods)
        fprintf(fid,['- %s: latent RMSE %.6f [%.6f, %.6f]; ', ...
            'observed RMSE %.6f [%.6f, %.6f]; time %.3f s ', ...
            '[%.3f, %.3f]; convergence %.2f.\n'],methods.method(q), ...
            methods.meanLatentRmse(q),methods.latentCiLower(q), ...
            methods.latentCiUpper(q),methods.meanObservedRmse(q), ...
            methods.observedCiLower(q),methods.observedCiUpper(q), ...
            methods.meanTotalSeconds(q),methods.timeCiLower(q), ...
            methods.timeCiUpper(q),methods.convergenceRate(q));
    end
    fprintf(fid,'\n## Paired RMSE differences\n\n');
    for q = 1:height(pairs)
        fprintf(fid,['- %s minus %s: latent %.6f [%.6f, %.6f]; ', ...
            'observed %.6f [%.6f, %.6f].\n'],pairs.method(q), ...
            pairs.comparator(q),pairs.meanLatentDifference(q), ...
            pairs.latentDifferenceCiLower(q), ...
            pairs.latentDifferenceCiUpper(q), ...
            pairs.meanObservedDifference(q), ...
            pairs.observedDifferenceCiLower(q), ...
            pairs.observedDifferenceCiUpper(q));
    end
    fprintf(fid,'\n## Adaptive timing contrasts\n\n');
    for q = 1:height(timing)
        fprintf(fid,['- Versus %s: mean saved %.3f s [%.3f, %.3f]; ', ...
            'mean saved fraction %.4f [%.4f, %.4f].\n'], ...
            timing.comparator(q),timing.meanSecondsSaved(q), ...
            timing.secondsSavedCiLower(q),timing.secondsSavedCiUpper(q), ...
            timing.meanFractionSaved(q),timing.fractionSavedCiLower(q), ...
            timing.fractionSavedCiUpper(q));
    end
end

function ensure_directory_local(directory)
    if ~isfolder(directory)
        mkdir(directory);
    end
end
