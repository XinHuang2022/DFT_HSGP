function summary = summarize_retained_direct_multiplier(experimentRoot)
%SUMMARIZE_RETAINED_DIRECT_MULTIPLIER Read preserved matched-run outputs.

    inputRoot = fullfile(experimentRoot, 'input');
    allMethods = readtable(fullfile(inputRoot, 'four_method_summary.csv'), ...
        'TextType', 'string');
    manuscriptMethods = ["numerical_inversion_reference"; ...
        "fixed_high_rank_hsgp"; "frequency_adaptive_hsgp"];
    summary.methodSummary = allMethods( ...
        ismember(allMethods.method, manuscriptMethods), :);
    summary.replicateMetrics = readtable(fullfile(inputRoot, ...
        'replicate_metrics.csv'), 'TextType', 'string');
    summary.rankProfiles = readtable(fullfile(inputRoot, ...
        'adaptive_rank_profiles.csv'), 'TextType', 'string');
    summary.source = "retained matched 20-replicate run";
end

