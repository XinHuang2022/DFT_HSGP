function checks = run_pacific_unit_checks(preset)
%RUN_PACIFIC_UNIT_CHECKS Data-independent algebra and orientation checks.

    if nargin < 1 || isempty(preset)
        preset = 'smoke3';
    end
    cfg = pacific_eq9_config(preset);
    outputDir = fullfile(cfg.paths.resultsRoot, 'unit_checks');
    ensure_directory(outputDir);
    theta.signalVariance = 3.2;
    theta.alpha = 0.42;
    theta.beta = 0.28;
    theta.epsilon = 0.47;
    theta.noiseVariance = 0.35;
    theta.nu = cfg.model.nu;
    z = eq9_pack_parameters(theta, cfg);
    recovered = eq9_unpack_parameters(z, cfg);
    roundTripError = max(abs([recovered.signalVariance - ...
        theta.signalVariance, recovered.alpha - theta.alpha, ...
        recovered.beta - theta.beta, recovered.epsilon - ...
        theta.epsilon, recovered.noiseVariance - theta.noiseVariance]));

    fixedZ = pacific_pack_fixed_epsilon(theta, cfg);
    fixedRecovered = pacific_unpack_fixed_epsilon(fixedZ, 1, cfg);
    fixedRoundTripError = max(abs([fixedRecovered.signalVariance - ...
        theta.signalVariance, fixedRecovered.alpha - theta.alpha, ...
        fixedRecovered.beta - theta.beta, fixedRecovered.noiseVariance - ...
        theta.noiseVariance]));

    T = 16;
    xy = [0, 0; 1.7, 0; 0.5, 2.0];
    design = make_pacific_mean_design(T, xy);
    reconstructedDesign = ifft(design.Afrequency, [], 1) * sqrt(T);
    designError = max(abs(reconstructedDesign - design.B), [], 'all');
    lambda = eq9_fourier_frequencies(T);
    covariance = eq9_frequency_covariances(theta, xy, lambda, cfg, true);
    hermitianError = 0;
    minimumEigenvalue = Inf;
    for k = 1:T
        matrix = covariance(:, :, k);
        hermitianError = max(hermitianError, ...
            max(abs(matrix - matrix'), [], 'all'));
        minimumEigenvalue = min(minimumEigenvalue, ...
            min(real(eig((matrix + matrix') / 2))));
    end
    cross = eq9_cross_frequency_covariances(theta, xy(1, :), ...
        xy(2:3, :), lambda, cfg);
    crossShapePassed = isequal(size(cross), [1, 2, T]);

    [conditioningMeanError, conditioningVarianceError] = ...
        exact_conditioning_equivalence(theta, cfg);
    [foldCoveragePassed, foldBalancePassed, foldAdjacencyPassed] = ...
        balanced_fold_checks();

    check = ["parameter transform round trip"; ...
        "fixed-epsilon transform round trip"; ...
        "mean-design FFT round trip"; ...
        "spectral covariance Hermitian"; ...
        "spectral covariance positive definite"; ...
        "target-observed cross-covariance orientation"; ...
        "frequency/time conditional mean equivalence"; ...
        "frequency/time conditional variance equivalence"; ...
        "balanced folds cover every site once"; ...
        "balanced folds have equal target counts"; ...
        "orthogonal neighbors assigned different folds"];
    passed = [roundTripError < 1e-9; fixedRoundTripError < 1e-9; ...
        designError < 1e-10; hermitianError < 1e-10; ...
        minimumEigenvalue > 0; crossShapePassed; ...
        conditioningMeanError < 1e-9; ...
        conditioningVarianceError < 1e-9; foldCoveragePassed; ...
        foldBalancePassed; foldAdjacencyPassed];
    value = [roundTripError; fixedRoundTripError; designError; ...
        hermitianError; minimumEigenvalue; double(crossShapePassed); ...
        conditioningMeanError; conditioningVarianceError; ...
        double(foldCoveragePassed); double(foldBalancePassed); ...
        double(foldAdjacencyPassed)];
    tolerance = [1e-9; 1e-9; 1e-10; 1e-10; 0; 1; 1e-9; 1e-9; ...
        1; 1; 1];
    checks = table(check, passed, value, tolerance);
    writetable(checks, fullfile(outputDir, 'unit_checks.csv'));
    if ~all(passed)
        error('One or more Pacific unit checks failed.');
    end
    fprintf('All %d Pacific unit checks passed.\n', height(checks));
end

function [coveragePassed, balancePassed, adjacencyPassed] = ...
        balanced_fold_checks()
    [columnGrid, rowGrid] = meshgrid(1:5, 1:5);
    n = numel(rowGrid);
    sites = table((1:n)', strings(n, 1), rowGrid(:), ...
        columnGrid(:), zeros(n, 1), zeros(n, 1), ...
        'VariableNames', {'subsetIndex', 'stationID', 'gridRow', ...
        'gridColumn', 'xKm', 'yKm'});
    sites.stationID = "S" + string((1:n)');
    subset.siteTable = sites;
    [~, manifest, summary] = make_pacific_cv_folds( ...
        subset, 'balanced5');
    coveragePassed = all(manifest.foldID >= 1 & ...
        manifest.foldID <= 5) && numel(manifest.foldID) == n;
    balancePassed = all(summary.nTargetSites == 5);
    folds = reshape(manifest.foldID, 5, 5);
    adjacencyPassed = all(folds(1:end-1, :) ~= folds(2:end, :), ...
        'all') && all(folds(:, 1:end-1) ~= folds(:, 2:end), 'all');
end

function [meanError, varianceError] = exact_conditioning_equivalence( ...
        theta, cfg)
% Compare frequency conditioning with full block-circulant conditioning.

    T = 8;
    xyObserved = [0, 0; 1.2, 0.4];
    xyTarget = [0.45, 0.75];
    lambda = eq9_fourier_frequencies(T);
    Foo = eq9_frequency_covariances(theta, xyObserved, ...
        lambda, cfg, true);
    Fto = eq9_cross_frequency_covariances(theta, xyTarget, ...
        xyObserved, lambda, cfg);
    Ftt = eq9_frequency_covariances(theta, xyTarget, ...
        lambda, cfg, false);

    observed = reshape(sin((1:(2 * T))' / 3), T, 2);
    Dobserved = fft(observed, [], 1) / sqrt(T);
    Dtarget = complex(zeros(T, 1));
    frequencyVariance = zeros(T, 1);
    for k = 1:T
        Dtarget(k) = Fto(:, :, k) * (Foo(:, :, k) \ ...
            Dobserved(k, :).');
        conditional = Ftt(:, :, k) - Fto(:, :, k) * ...
            (Foo(:, :, k) \ Fto(:, :, k)');
        frequencyVariance(k) = real(conditional);
    end
    frequencyMean = real(ifft(Dtarget, [], 1) * sqrt(T));
    frequencyPointVariance = mean(frequencyVariance);

    gammaOO = ifft(Foo, [], 3);
    gammaTO = ifft(Fto, [], 3);
    gammaTT = ifft(Ftt, [], 3);
    Coo = complex(zeros(2 * T));
    Cto = complex(zeros(T, 2 * T));
    Ctt = complex(zeros(T));
    for t = 1:T
        targetRow = t;
        observedRows = (2 * (t - 1) + 1):(2 * t);
        for u = 1:T
            observedColumns = (2 * (u - 1) + 1):(2 * u);
            lag = mod(t - u, T) + 1;
            Coo(observedRows, observedColumns) = gammaOO(:, :, lag);
            Cto(targetRow, observedColumns) = gammaTO(:, :, lag);
            Ctt(t, u) = gammaTT(:, :, lag);
        end
    end
    observedVector = reshape(observed.', [], 1);
    exactMean = real(Cto * (Coo \ observedVector));
    exactConditional = Ctt - Cto * (Coo \ Cto');
    exactVariance = real(diag(exactConditional));

    meanError = max(abs(exactMean - frequencyMean));
    varianceError = max(abs(exactVariance - frequencyPointVariance));
end
