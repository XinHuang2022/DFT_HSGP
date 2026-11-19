clear; clc;
result=run_ghcnh_trace_resource_benchmark(1);
disp(rmfield(result,{'selectedJ','spotChecks','details'}));
disp(result.spotChecks);
