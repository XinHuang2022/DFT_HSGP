clear; clc;
continuation=continue_ghcnh_algorithm2_fold(1,3);
disp(continuation.history);
if continuation.rankStable, disp(continuation.metrics); end
