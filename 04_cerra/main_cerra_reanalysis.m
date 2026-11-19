%% CERRA North Sea reanalysis: choose an execution mode
% "retained" displays the preserved 31-day and 90-day results.
% "preflight_31day" or "preflight_90day" validates one data window.
% "paper_31day" or "paper_90day" reruns the corresponding workflow.
executionMode = "retained";
numberOfWorkers = 6;

%% Initialize the local reproducibility package
experimentRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(experimentRoot, 'src'));
fprintf('CERRA North Sea mode: %s\n', executionMode);

%% Inspect the extracted dataset
dataDirectory = fullfile(experimentRoot, 'input', 'data');
siteManifest = readtable(fullfile(dataDirectory, ...
    'site_manifest_full_native.csv'), 'VariableNamingRule', 'preserve');
timeManifest = readtable(fullfile(dataDirectory, ...
    'time_manifest_6hourly.csv'), 'VariableNamingRule', 'preserve');
fprintf('Extracted CERRA data: %d marine sites and %d time points.\n', ...
    height(siteManifest), height(timeManifest));

%% Run the selected workflow
result = run_cerra_experiment(executionMode, numberOfWorkers);

%% Inspect the principal results
if isstruct(result) && isfield(result, 'summary31Day')
    disp(result.summary31Day(:, {'method','M','nTrainingSites','pooledRmse'}));
end
if isstruct(result) && isfield(result, 'summary90Day')
    disp(result.summary90Day(:, {'method','M','nTrainingSites','pooledRmse'}));
end

%% Locate generated files
if isstruct(result) && isfield(result, 'outputDir')
    fprintf('Generated results: %s\n', result.outputDir);
else
    fprintf('No new fit was run. Preserved summaries are in input/retained/.\n');
end
