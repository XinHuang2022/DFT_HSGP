clear; clc;
methods=["exact_kl","trace_bound"];
for foldID=1:5
    for method=methods
        fprintf('[%s] starting fold %d method %s\n', ...
            char(datetime('now')),foldID,method);
        task=run_ghcnh_expanded_method(method,foldID,5);
        fprintf('[%s] finished fold %d method %s stable=%d\n', ...
            char(datetime('now')),foldID,method,task.rankStable);
        if ~task.rankStable
            error('Expanded comparison did not stabilize for fold %d, %s.', ...
                foldID,method);
        end
    end
end
summary=bootstrap_ghcnh_expanded_cv();
disp(summary);
