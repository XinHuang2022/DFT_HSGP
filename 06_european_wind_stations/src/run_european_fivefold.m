function result = run_european_fivefold(mode, numberOfWorkers)
%RUN_EUROPEAN_FIVEFOLD Run five whole-site folds with outer parallelism.
% Each worker owns one fold. Frequency-level routines remain serial inside
% a worker, which avoids nested parallel pools.

    if nargin < 1 || isempty(mode), mode = "paper"; end
    if nargin < 2 || isempty(numberOfWorkers), numberOfWorkers = 5; end
    mode = lower(string(mode));
    if numberOfWorkers < 1 || numberOfWorkers ~= floor(numberOfWorkers)
        error('numberOfWorkers must be a positive integer.');
    end
    numberOfWorkers = min(5, numberOfWorkers);

    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= numberOfWorkers
        if ~isempty(pool), delete(pool); end
        pool = parpool(numberOfWorkers);
    end
    fprintf(['European station CV: five folds on %d workers; ', ...
        'one fold is assigned to each worker when five workers are used.\n'], ...
        pool.NumWorkers);

    switch mode
        case "fixed_paper"
            run_method("fixed_reference");
            result.fixedSummary = bootstrap_ghcnh_expanded_cv("fixed_reference");
        case "algorithm1"
            ensure_shared_pilots();
            run_method("exact_kl");
            result.workflow = "Algorithm 1 completed for five folds";
        case "algorithm2"
            ensure_shared_pilots();
            run_method("solve_trace");
            result.workflow = "Algorithm 2 completed for five folds";
        case "paper"
            run_method("fixed_reference");
            run_method("exact_kl");
            run_method("solve_trace");
            result.fixedSummary = bootstrap_ghcnh_expanded_cv("fixed_reference");
            result.relaxedSummary = summarize_ghcnh_relaxed_cycle_cv();
        otherwise
            error('Unknown European five-fold mode: %s', mode);
    end

    function ensure_shared_pilots()
        [~,~,~,~,E] = ghcnh_expanded_setup();
        missing = false(5, 1);
        for foldID = 1:5
            missing(foldID) = ~isfile(fullfile(E.resultsDir, ...
                sprintf('fold_%02d', foldID), 'shared_pilot_fit.mat'));
        end
        if any(missing)
            fprintf(['Shared pilot fits are absent for folds %s. ', ...
                'Running fixed-reference folds first.\n'], ...
                mat2str(find(missing)'));
            run_method("fixed_reference");
        end
    end

    function run_method(method)
        queue = parallel.pool.DataQueue;
        afterEach(queue, @print_progress);
        fprintf('Starting %s across five folds.\n', method);
        parfor foldID = 1:5
            started = tic;
            if method == "fixed_reference"
                run_ghcnh_fixed_reference_fold(foldID);
            elseif method == "exact_kl"
                run_ghcnh_relaxed_cycle_fold("exact_kl", foldID, 4, 0.01);
            else
                run_ghcnh_relaxed_cycle_fold("solve_trace", foldID, 3, 0.05);
            end
            send(queue, struct('method', method, 'fold', foldID, ...
                'minutes', toc(started) / 60));
        end
    end
end

function print_progress(update)
    fprintf('[%s] fold %d completed in %.1f minutes.\n', ...
        update.method, update.fold, update.minutes);
end
