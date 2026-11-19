%% European station wind experiment: choose an execution mode
% "retained" displays preserved five-fold summaries without refitting.
% "preflight" validates the extracted panel, station order, and folds.
% "fixed_paper", "algorithm1", and "algorithm2" run one procedure.
% "paper" runs fixed rank, Algorithm 1, and Algorithm 2 in that order.
executionMode = "retained";
numberOfWorkers = 5;

%% Initialize the local reproducibility package
experimentRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(experimentRoot, 'src'));
fprintf('European station wind mode: %s\n', executionMode);

%% Inspect the extracted panel
panelFile = fullfile(experimentRoot, 'input', 'data', ...
    'ghcnh_complete_panel_expanded.mat');
panel = load(panelFile, 'subset');
fprintf('Extracted station panel: %d times by %d sites.\n', ...
    size(panel.subset.U, 1), size(panel.subset.U, 2));

%% Run the selected workflow
result = run_european_wind_experiment(executionMode, numberOfWorkers);

%% Inspect the principal results
if isstruct(result) && isfield(result, 'fivefoldSummary')
    disp(result.fivefoldSummary);
end
if isstruct(result) && isfield(result, 'relaxedSummary')
    disp(result.relaxedSummary);
end

%% Locate generated files
if isstruct(result) && isfield(result, 'outputDir')
    fprintf('Generated results: %s\n', result.outputDir);
else
    fprintf('No new fit was run. Preserved summaries are in input/retained/.\n');
end
