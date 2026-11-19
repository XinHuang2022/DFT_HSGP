function runSummary=run_ghcnh_relaxed_cycle_parallel(method)
%RUN_GHCNH_RELAXED_CYCLE_PARALLEL Run relaxed folds 2--5 concurrently.
%
% Algorithm 1 (exact_kl): at most four adaptive fits and at most 1% of
% retained frequencies requesting further enlargement.
% Algorithm 2 (solve_trace): at most three adaptive fits and at most 5% of
% retained frequencies requesting further enlargement.

    method=string(method);
    if method=="exact_kl"
        maxCycles=4; maxEnlargedFraction=0.01;
    elseif method=="solve_trace"
        maxCycles=3; maxEnlargedFraction=0.05;
    else
        error('method must be exact_kl or solve_trace.');
    end

    foldIDs=2:5; requestedWorkers=numel(foldIDs);
    pool=gcp('nocreate');
    if isempty(pool) || pool.NumWorkers~=requestedWorkers
        if ~isempty(pool), delete(pool); end
        pool=parpool('local',requestedWorkers);
    end
    queue=parallel.pool.DataQueue;
    afterEach(queue,@(message) fprintf('[%s] %s\n', ...
        char(datetime('now')),message));
    tasks=cell(numel(foldIDs),1); failures=strings(numel(foldIDs),1);
    blockTimer=tic;
    parfor q=1:numel(foldIDs)
        foldID=foldIDs(q);
        try
            send(queue,sprintf(['starting fold %d method %s: maxCycles=%d ' ...
                'allowedFraction=%.3f'],foldID,method,maxCycles, ...
                maxEnlargedFraction));
            tasks{q}=run_ghcnh_relaxed_cycle_fold(method,foldID, ...
                maxCycles,maxEnlargedFraction);
            send(queue,sprintf(['finished fold %d method %s cycles=%d ' ...
                'RMSE=%.4f first-use=%.2f min'],foldID,method, ...
                tasks{q}.stopping.cyclesUsed,tasks{q}.metrics.rmse, ...
                tasks{q}.timing.equivalentFirstUseSeconds/60));
        catch ME
            failures(q)=sprintf('fold %d failed: %s',foldID,ME.message);
            send(queue,failures(q));
        end
    end
    parallelBlockSeconds=toc(blockTimer);
    delete(pool);
    if any(strlength(failures)>0)
        error('Parallel relaxed-cycle run incomplete:\n%s', ...
            strjoin(failures(strlength(failures)>0),newline));
    end

    cycleOneSeconds=zeros(numel(tasks),1);
    for q=1:numel(tasks)
        cycleOneSeconds(q)=tasks{q}.history.fitSeconds(1)+ ...
            tasks{q}.history.auditSeconds(1);
    end
    runSummary=table(method,maxCycles,maxEnlargedFraction, ...
        requestedWorkers,parallelBlockSeconds,parallelBlockSeconds/60, ...
        max(cycleOneSeconds),mean(cycleOneSeconds), ...
        'VariableNames',{'method','maxCycles','maxEnlargedFraction', ...
        'parallelWorkers','parallelBlockSeconds','parallelBlockMinutes', ...
        'maximumFoldCycleOneSeconds','meanFoldCycleOneSeconds'});

    [~,~,~,~,E]=ghcnh_expanded_setup();
    tag=sprintf('maxcycles_%d_fraction_%04d',maxCycles, ...
        round(10000*maxEnlargedFraction));
    outputDir=fullfile(E.resultsDir,'relaxed_cycle_sensitivity', ...
        ['parallel_' tag],char(method));
    if ~isfolder(outputDir), mkdir(outputDir); end
    writetable(runSummary,fullfile(outputDir,'parallel_wall_time.csv'));
    save(fullfile(outputDir,'parallel_run_summary.mat'),'runSummary', ...
        'foldIDs','maxCycles','maxEnlargedFraction');
    disp(runSummary);
end
