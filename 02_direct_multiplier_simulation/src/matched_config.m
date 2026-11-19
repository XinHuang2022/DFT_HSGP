function C = matched_config(profile)
%MATCHED_CONFIG Configuration for the matched four-method experiment.

    if nargin < 1 || isempty(profile)
        profile = 'smoke';
    end
    sourceRoot = fileparts(mfilename('fullpath'));
    experimentRoot = fileparts(sourceRoot);

    C.version = '0.2.0';
    C.profile = char(profile);
    C.suiteRoot = experimentRoot;
    C.baseSuiteRoot = experimentRoot;
    C.baseSourceRoot = fullfile(sourceRoot, 'model');
    C.outputRoot = fullfile(experimentRoot, 'results', 'raw');
    C.seed = 918273;
    % The replicate driver replaces these with the trusted Table-4 seeds
    % C.seed+10000+r and C.seed+20000+r, respectively.
    C.simulationSeed = C.seed + 10001;
    C.optimizerSeed = C.seed + 20001;
    C.figureVisible = 'off';
    C.writeReplicateFigures = false;
    % Per-replicate evidence is written as CSV. Saving twenty concurrent
    % HDF5 MAT files on the project drive is unnecessary for aggregation and
    % can fail under MATLAB process-based parallelism.
    C.saveReplicateMat = false;

    % theta = [vf,rho,kappa0,gamma,vn]. These definitions and bounds are
    % inherited unchanged from the established direct-multiplier suite.
    C.truth = [1.0, 0.85, 5.0, 2.0, 0.05];
    C.lower = [0.25,0.20,4.0,0.05,0.005];
    C.upper = [3.00,0.98,7.0,3.50,0.300];
    C.initial = [0.8,0.65,5.5,1.5,0.08];
    C.coverage = 0.95;

    % Fixed bounded domain and nested ranks. M_init is used only to obtain
    % the preliminary covariance-parameter estimate required by Algorithm B;
    % it is not a competing final method. The fixed high-rank comparator uses
    % the stabilized value selected from referenceM.
    C.c = 2;
    C.Minit = 121;
    C.candidateM = [25,49,81,121,169,225,324];
    C.referenceM = [400,625,900];
    C.klTolerance = 0.01;
    C.referenceTolerance = 0.001;
    % Rank-audit rule. The exact-KL choice reproduces the established
    % Algorithm-1 experiment. The trace-certificate choice implements
    % Algorithm 2 while leaving simulation, fitting, and evaluation fixed.
    C.auditRule = 'exact_kl';
    % The iteration limit is a safety cap, not a tuning parameter.  The
    % algorithm stops earlier as soon as a post-fit audit requests no rank
    % enlargement.
    C.maxAdaptiveIterations = 3;
    C.nObjectiveTimingReplicates = 3;

    % Numerical-inversion reference. It is used only to simulate and to
    % evaluate the fitted methods, never to select the Algorithm-B ranks.
    C.ref.xMax = 100;
    C.ref.panelWidth = 1;
    C.ref.order = 24;
    C.ref.dz = 0.01;
    C.ref.zMax = 50;
    C.ref.checkOrder = 32;

    C.opt.maxIter = 260;
    C.opt.maxEval = 1000;
    C.opt.tolX = 2e-4;
    C.opt.tolFun = 2e-4;
    C.opt.startScale = 0.55;
    C.starts = 4;
    % The first adaptive fit retains the full multistart design.  Subsequent
    % rank-correction fits are warm-started from the preceding estimate and
    % ordinarily need only that single start.
    C.correctionStarts = 1;

    switch lower(C.profile)
        case 'smoke'
            C.nSide = 7;
            C.T = 24;
            C.nTest = 10;
            C.Minit = 25;
            C.candidateM = [9,16,25,36,49];
            C.referenceM = [64,81,100];
            C.klTolerance = 0.05;
            C.referenceTolerance = 0.01;
            C.starts = 1;
            C.opt.maxIter = 10;
            C.opt.maxEval = 60;
            C.ref.xMax = 60;
            C.ref.order = 16;
            C.ref.checkOrder = 22;
            C.ref.dz = 0.02;
        case 'pilot'
            C.nSide = 12;
            C.T = 96;
            C.nTest = 20;
        case 'confirm'
            C.nSide = 15;
            C.T = 128;
            C.nTest = 45;
            C.opt.maxIter = 500;
            C.opt.maxEval = 2000;
        otherwise
            error('Unknown Algorithm-B profile: %s', C.profile);
    end
    C.nSites = C.nSide^2;
end

