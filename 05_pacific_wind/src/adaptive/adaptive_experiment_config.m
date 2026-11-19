function C = adaptive_experiment_config(dataset)
%ADAPTIVE_EXPERIMENT_CONFIG Prespecified Algorithm-1 experiment settings.

    dataset = lower(string(dataset));
    here = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(here));
    C.dataset = dataset;
    C.version = "0.1.0";
    outputRoot = getenv('PACIFIC_REPRO_OUTPUT_DIR');
    if isempty(outputRoot)
        outputRoot = fullfile(projectRoot, 'results', 'current');
    end
    C.outputRoot = fullfile(outputRoot, 'adaptive_rank');
    C.dataDir = fullfile(projectRoot, 'input', 'data');
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
    C.reuseTrustedPilot = false;
    C.auditOnly = false;

    switch dataset
        case "pacific"
            C.sourceRoot = fullfile(projectRoot, 'src', 'model');
            C.fixedResultRoot = fullfile(projectRoot, 'input', ...
                'retained', 'fixed_rank');
            C.fixedMethod = "hsgp_eq9";
            C.initialJ = 48;
            C.trustedPilotDir = '';
            C.trustedPilotName = 'fold_%02d.mat';
        otherwise
            error(['This package contains only the Pacific experiment. ', ...
                'Unknown data set %s.'], dataset);
    end
end
