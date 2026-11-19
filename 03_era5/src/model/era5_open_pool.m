function [parallelUsed, poolWorkers] = era5_open_pool(cfg, nTasks)
%ERA5_OPEN_POOL Open at most one local pool; never nest parallel loops.

    parallelUsed = false; poolWorkers = 1;
    if ~cfg.parallel.useParallel, return; end
    try
        if ~license('test', 'Distrib_Computing_Toolbox')
            warning('Parallel Computing Toolbox unavailable; running serially.');
            return;
        end
        pool = gcp('nocreate');
        requested = min(cfg.parallel.numWorkers, nTasks);
        if isempty(pool) && cfg.parallel.autoStartPool
            pool = parpool(requested);
        end
        if isempty(pool), error('Parallel execution requested but no pool is open.'); end
        parallelUsed = true; poolWorkers = pool.NumWorkers;
    catch ME
        warning('ERA5:ParallelFallback', 'Parallel startup failed: %s', ME.message);
    end
end
