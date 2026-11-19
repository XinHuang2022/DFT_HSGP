function task = run_pacific_hsgp_fold_task(cfg, subset, fold, basis)
%RUN_PACIFIC_HSGP_FOLD_TASK Resume-safe wrapper for one CV fold.

    taskDir = fullfile(cfg.paths.hsgpCvDir, 'tasks');
    ensure_directory(taskDir);
    taskFile = fullfile(taskDir, sprintf('fold_%02d.mat', fold.foldID));
    if isfile(taskFile)
        loaded = load(taskFile, 'task');
        if isfield(loaded, 'task') && isfield(loaded.task, 'metrics') && ...
                isfield(loaded.task, 'fitDiagnostics') && ...
                isfield(loaded.task, 'version') && ...
                loaded.task.version == string(cfg.version) && ...
                isfield(loaded.task, 'basisFeatures') && ...
                loaded.task.basisFeatures == basis.nFeatures && ...
                loaded.task.basisBoundaryFactor == basis.boundaryFactor
            fprintf('HSGP fold %d loaded from checkpoint.\n', fold.foldID);
            task = loaded.task;
            return;
        end
    end
    fprintf('HSGP fold %d start on worker %d.\n', ...
        fold.foldID, current_worker_id());
    task = run_pacific_hsgp_fold(cfg, subset, fold.targetIndex, ...
        fold.foldID, basis);
    save(taskFile, 'task', '-v7.3');
    fprintf('HSGP fold %d complete.\n', fold.foldID);
end

function id = current_worker_id()
    id = 0;
    try
        currentTask = getCurrentTask();
        if ~isempty(currentTask), id = currentTask.ID; end
    catch
        id = 0;
    end
end
