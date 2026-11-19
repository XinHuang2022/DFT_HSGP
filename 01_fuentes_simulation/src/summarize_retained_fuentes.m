function summary = summarize_retained_fuentes(experimentRoot)
%SUMMARIZE_RETAINED_FUENTES Read the preserved 20-replicate result tables.

    inputRoot = fullfile(experimentRoot, 'input');
    summary.predictionSummary = readtable(fullfile(inputRoot, ...
        'prediction_rmse_summary.csv'), 'TextType', 'string');
    summary.rankSummary = readtable(fullfile(inputRoot, ...
        'rank_profile_overall_summary.csv'), 'TextType', 'string');
    summary.rankProfiles = readtable(fullfile(inputRoot, ...
        'fitted_rank_profiles.csv'), 'TextType', 'string');
    summary.source = "retained manuscript run";
end

