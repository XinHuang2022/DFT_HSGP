%% Pacific wind application: choose an execution mode
% "retained" displays preserved results without refitting.
% "preflight" validates the local data, folds, and model interfaces.
% "fixed_paper" reruns the fixed-rank five-fold HSGP workflow.
% "adaptive_paper" reruns the frequency-adaptive five-fold workflow.
executionMode = "retained";
numberOfWorkers = 5;

%% Initialize the local reproducibility package
experimentRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(experimentRoot, 'src'));
fprintf('Pacific wind mode: %s\n', executionMode);

%% Inspect the extracted dataset
dataDirectory = fullfile(experimentRoot, 'input', 'data');
wind = readmatrix(fullfile(dataDirectory, 'ncepu.dat'), 'FileType', 'text');
coordinates = readmatrix(fullfile(dataDirectory, 'nceplatlon.dat'), ...
    'FileType', 'text');
fprintf('Pacific data files: wind %d-by-%d; coordinates %d-by-%d.\n', ...
    size(wind, 1), size(wind, 2), size(coordinates, 1), size(coordinates, 2));

%% Run the selected workflow
result = run_pacific_wind_experiment(executionMode, numberOfWorkers);

%% Inspect the principal results
if isstruct(result) && isfield(result, 'fixedMethodSummary')
    disp(result.fixedMethodSummary);
end
if isstruct(result) && isfield(result, 'adaptiveMethodSummary')
    disp(result.adaptiveMethodSummary);
end

%% Locate generated files
if isstruct(result) && isfield(result, 'outputDir')
    fprintf('Generated results: %s\n', result.outputDir);
else
    fprintf('No new fit was run. Preserved summaries are in input/retained/.\n');
end
