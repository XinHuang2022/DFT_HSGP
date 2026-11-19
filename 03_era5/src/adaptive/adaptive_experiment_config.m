function C = adaptive_experiment_config(dataset)
%ADAPTIVE_EXPERIMENT_CONFIG Prespecified Algorithm-1 experiment settings.

    dataset = lower(string(dataset));
    here = fileparts(mfilename('fullpath'));
    experimentRoot = fileparts(fileparts(here));
    C.dataset = dataset;
    C.version = "0.1.0";
    C.outputRoot = fullfile(experimentRoot, 'results', 'adaptive_current');
    C.dataDir = fullfile(experimentRoot, 'input', 'data');
    C.nWorkers = 5;
    C.epsilonKlPerSite = 0.01;
    C.epsilonReferencePerSite = 0.001;
    C.initialJ = 32;
    C.candidateJ = [16 20 24 28 32 36 40 44 48 52];
    C.referenceJ = [48 56 64];
    C.maxAdaptiveFits = 3;
    C.progressEveryFrequencies = 20;
    C.bootstrapReplicates = 2000;
    C.bootstrapSeed = 94117;
    C.runFolds = 1:5;
    C.smoke = false;
    C.resume = true;
    C.reuseTrustedPilot = true;
    C.auditOnly = false;

    switch dataset
        case "era5"
            C.sourceRoot = fullfile(experimentRoot, 'src', 'model');
            C.fixedResultRoot = fullfile(experimentRoot, 'input', ...
                'retained', 'fixed_rank', 'stage3_paired_cv');
            C.fixedMethod = "hsgp_J48";
            C.initialJ = 32;
            C.trustedPilotDir = fullfile(experimentRoot, 'input', ...
                'trusted_checkpoints_not_distributed');
            C.trustedPilotName = 'hsgp_J32_fold_%02d_fit.mat';
            C.reuseTrustedPilot = false;
        case "pacific"
            error(['The Pacific implementation is packaged separately; ' ...
                'this folder contains only the ERA5 data contract.']);
        otherwise
            error('Unknown data set %s.', dataset);
    end
end
