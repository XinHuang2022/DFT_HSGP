function C = adaptive_experiment_config(dataset)
%ADAPTIVE_EXPERIMENT_CONFIG Prespecified Algorithm-1 experiment settings.

    dataset = lower(string(dataset));
    here = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(here));
    gaussianProcessRoot = fileparts(fileparts(projectRoot));
    windRoot = fullfile(gaussianProcessRoot, 'wind_field_implement', 'Comp');
    C.dataset = dataset;
    C.version = "0.1.0";
    C.outputRoot = fullfile(here, 'Results', char(dataset));
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
            C.sourceRoot = fullfile(windRoot, 'ERA5', ...
                'era5_northsea_eq9_suite_v0_1_0');
            C.fixedResultRoot = fullfile(C.sourceRoot, 'Results', ...
                'stage3_paired_cv');
            C.fixedMethod = "hsgp_J48";
            C.initialJ = 32;
            C.trustedPilotDir = fullfile(C.fixedResultRoot, 'tasks');
            C.trustedPilotName = 'hsgp_J32_fold_%02d_fit.mat';
        case "pacific"
            C.sourceRoot = fullfile(windRoot, 'Pacific', ...
                'pacific_wind_eq9_suite_v0_3_1_pdc', ...
                'pacific_wind_eq9_suite_v0_3_1_hsgp');
            C.fixedResultRoot = fullfile(C.sourceRoot, 'Results', ...
                'full17_cv5_hsgp_cv');
            C.fixedMethod = "hsgp_eq9";
            C.initialJ = 48;
            C.trustedPilotDir = fullfile(C.fixedResultRoot, 'tasks');
            C.trustedPilotName = 'fold_%02d.mat';
        otherwise
            error('Unknown data set %s.', dataset);
    end
end
