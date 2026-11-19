clear; clc;
task=continue_ghcnh_algorithm1_fold(1,3);
disp(task.timing);
disp(task.history);
if task.rankStable
    disp(task.metrics);
end
