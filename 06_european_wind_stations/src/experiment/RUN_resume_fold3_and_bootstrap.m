clear; clc;
fprintf('[%s] resuming fold 3 with a six-cycle cap\n',char(datetime('now')));
task=run_ghcnh_expanded_method("exact_kl",3,6);
fprintf('[%s] fold 3 finished stable=%d\n', ...
    char(datetime('now')),task.rankStable);
if ~task.rankStable
    error('Algorithm 1 did not stabilize for fold 3 after six cycles.');
end
summary=bootstrap_ghcnh_expanded_cv("exact_kl");
disp(summary);
