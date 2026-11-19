%% Direct-multiplier simulation: choose an execution mode
% "retained" displays the preserved 20-replicate results.
% "preflight" checks paths, inputs, and nested rank settings.
% "smoke" runs one small end-to-end software check.
% "paper" reruns the complete 20-replicate experiment.
executionMode = "retained";
numberOfWorkers = 5;
auditRule = "exact_kl";
klTolerance = 0.01;

%% Initialize paths and configuration
experimentRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(experimentRoot, 'src'));
addpath(fullfile(experimentRoot, 'src', 'model'));
fprintf('Direct-multiplier simulation mode: %s\n', executionMode);

%% Validate the local experiment package
preflightReport = preflight_direct_multiplier(experimentRoot);
disp(preflightReport);

%% Run the selected workflow
switch lower(executionMode)
    case "retained"
        result = summarize_retained_direct_multiplier(experimentRoot);
    case "preflight"
        result = preflightReport;
    case "smoke"
        cfg = matched_config('smoke');
        cfg.auditRule = char(auditRule);
        cfg.klTolerance = klTolerance;
        result = run_matched_replicate(cfg);
    case "paper"
        result.outputDir = run_matched_20(numberOfWorkers, ...
            auditRule, klTolerance);
    otherwise
        error('Unknown execution mode: %s', executionMode);
end

%% Inspect the principal results
if isstruct(result) && isfield(result, 'methodSummary')
    disp(result.methodSummary);
end
if isstruct(result) && isfield(result, 'rankProfiles')
    fprintf('Retained rank-profile rows: %d\n', height(result.rankProfiles));
end

%% Locate generated files
if isstruct(result) && isfield(result, 'outputDir')
    fprintf('Generated results: %s\n', result.outputDir);
else
    fprintf('No new scientific fit was run. Retained outputs are in input/.\n');
end

