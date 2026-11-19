%% Fuentes simulation: choose an execution mode
% "retained" displays the preserved manuscript results without refitting.
% "preflight" checks paths, inputs, and the KL calculation.
% "paper" reruns the 20-replicate fitted-parameter rank experiment.
executionMode = "retained";
numberOfWorkers = 5;

%% Initialize paths and configuration
experimentRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(experimentRoot, 'src'));
addpath(fullfile(experimentRoot, 'src', 'model'));
cfg = fitted_adaptive_config();
cfg.nWorkers = numberOfWorkers;
cfg.useParallel = numberOfWorkers > 1;
fprintf('Fuentes simulation mode: %s\n', executionMode);

%% Validate the local experiment package
preflightReport = preflight_fuentes_simulation(cfg);
disp(preflightReport);

%% Run the selected workflow
switch lower(executionMode)
    case "retained"
        result = summarize_retained_fuentes(experimentRoot);
    case "preflight"
        result = preflightReport;
    case "paper"
        result = run_fitted_adaptive_fuentes(cfg);
    otherwise
        error('Unknown execution mode: %s', executionMode);
end

%% Inspect the principal results
if isstruct(result) && isfield(result, 'predictionSummary')
    disp(result.predictionSummary);
end
if isstruct(result) && isfield(result, 'rankSummary')
    disp(result.rankSummary);
end

%% Locate generated files
if isstruct(result) && isfield(result, 'outputDir')
    fprintf('Generated results: %s\n', result.outputDir);
else
    fprintf('No new scientific fit was run. Retained outputs are in input/.\n');
end

