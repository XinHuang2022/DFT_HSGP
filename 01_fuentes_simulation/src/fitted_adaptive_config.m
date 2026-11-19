function cfg = fitted_adaptive_config()
%FITTED_ADAPTIVE_CONFIG Prespecified fitted-parameter rank experiment.

    sourceRoot = fileparts(mfilename('fullpath'));
    experimentRoot = fileparts(sourceRoot);
    cfg.version = '0.3.0';
    cfg.experimentRoot = experimentRoot;
    cfg.sourceSuiteRoot = fullfile(sourceRoot, 'model');
    cfg.oracleResultFile = fullfile(experimentRoot, 'input', ...
        'fuentes_reference_metadata.mat');
    cfg.oracleSourceDirectory = fullfile(sourceRoot, 'model');
    cfg.fittedParameterFile = fullfile(experimentRoot, 'input', ...
        'trusted_paired_parameter_recovery.csv');
    cfg.trustedPredictionFile = fullfile(experimentRoot, 'input', ...
        'trusted_paired_prediction_rmse.csv');
    cfg.parameterFileSha256 = ...
        'B4B703B130621E8A42057DFF8ED6AFA074DA892490003B61E34B26833A8FA0BB';
    cfg.predictionFileSha256 = ...
        '06B0F8AAE4E57BCB14FA4CD6A41AE9874F74E9463D96B953CC2706DC4580637E';

    cfg.fitRegime = "free_noise";
    cfg.fitEngine = "exact";
    cfg.klPerFittingSiteTolerance = 0.01;
    cfg.candidateJ = [6, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40];
    cfg.nReplicates = 20;
    cfg.bootstrapReplicates = 10000;
    cfg.bootstrapSeed = 314159;
    cfg.useParallel = true;
    cfg.nWorkers = 5;
    cfg.figureVisible = 'off';
    cfg.outputRoot = fullfile(experimentRoot, 'results');
end
