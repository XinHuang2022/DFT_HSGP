clear; clc;
for foldID=1:5
    fprintf('[%s] starting fold %d method exact_kl\n', ...
        char(datetime('now')),foldID);
    task=run_ghcnh_expanded_method("exact_kl",foldID,5);
    fprintf('[%s] finished fold %d method exact_kl stable=%d\n', ...
        char(datetime('now')),foldID,task.rankStable);
    if ~task.rankStable
        error('Algorithm 1 did not stabilize for fold %d.',foldID);
    end
end
summary=bootstrap_ghcnh_expanded_cv("exact_kl");
disp(summary);
