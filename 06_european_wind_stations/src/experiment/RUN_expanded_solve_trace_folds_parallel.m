clear; clc;
foldIDs=2:5;
requestedWorkers=numel(foldIDs);
pool=gcp('nocreate');
if isempty(pool) || pool.NumWorkers~=requestedWorkers
    if ~isempty(pool), delete(pool); end
    pool=parpool('local',requestedWorkers);
end
queue=parallel.pool.DataQueue;
afterEach(queue,@(message) fprintf('[%s] %s\n', ...
    char(datetime('now')),message));
tasks=cell(numel(foldIDs),1);
failures=strings(numel(foldIDs),1);
parfor q=1:numel(foldIDs)
    foldID=foldIDs(q);
    try
        send(queue,sprintf('starting fold %d method solve_trace',foldID));
        task=run_ghcnh_expanded_method("solve_trace",foldID,6);
        tasks{q}=task;
        if ~task.rankStable
            failures(q)=sprintf('fold %d did not stabilize',foldID);
        end
        send(queue,sprintf('finished fold %d stable=%d', ...
            foldID,task.rankStable));
    catch ME
        failures(q)=sprintf('fold %d failed: %s',foldID,ME.message);
        send(queue,failures(q));
    end
end
delete(pool);
if any(strlength(failures)>0)
    error('Parallel fold run incomplete:\n%s', ...
        strjoin(failures(strlength(failures)>0),newline));
end
summary=bootstrap_ghcnh_expanded_cv(["exact_kl","solve_trace"]);
disp(summary);
