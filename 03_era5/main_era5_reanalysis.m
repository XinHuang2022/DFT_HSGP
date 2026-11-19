%% ERA5 North Sea reanalysis: choose an execution mode
% "retained" displays the preserved results without refitting.
% "preflight" validates the extracted data, folds, and model interfaces.
% "fixed_paper" reruns the analytic and fixed-rank five-fold workflow.
% "adaptive_paper" reruns the frequency-adaptive five-fold workflow.
executionMode = "retained";
numberOfWorkers = 5;

%% Initialize the local reproducibility package
experimentRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(experimentRoot, 'src'));
fprintf('ERA5 North Sea mode: %s\n', executionMode);

%% Inspect the extracted dataset
dataDirectory = fullfile(experimentRoot, 'input', 'data');
siteManifest = readtable(fullfile(dataDirectory, ...
    'site_manifest_17x17.csv'), 'VariableNamingRule', 'preserve');
timeManifest = readtable(fullfile(dataDirectory, ...
    'time_manifest_6hourly.csv'), 'VariableNamingRule', 'preserve');
fprintf('Extracted ERA5 data: %d sites and %d time points.\n', ...
    height(siteManifest), height(timeManifest));

%% Run the selected workflow
result = run_era5_experiment(executionMode, numberOfWorkers);

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
