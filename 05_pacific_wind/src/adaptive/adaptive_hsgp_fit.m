function fit = adaptive_hsgp_fit(Y, xy, design, cfg, master, ...
        selectedJ, initialTheta, optimizer, seed, A)
%ADAPTIVE_HSGP_FIT Whittle REML with a frequency-specific nested basis.

    prepared = A.prepare(Y, xy, design);
    prepared.masterPhi = A.features(xy, master);
    selectedJ = selectedJ(:);
    K = floor(prepared.T / 2) + 1;
    if numel(selectedJ) ~= K
        error('selectedJ must contain one entry for each nonnegative frequency.');
    end
    objective = @(z) transformed_objective(z, prepared, cfg, master, ...
        selectedJ, A);
    z0 = A.pack(initialTheta, cfg, [], []);
    optimization = eq9_multistart_optimize(objective, z0, optimizer, seed);
    [objectiveValue, details] = objective(optimization.z);
    if isfield(details, 'errorMessage')
        error('Best adaptive fit is invalid: %s', details.errorMessage);
    end
    theta = A.unpack(optimization.z, cfg, [], []);
    [boundaryHit, minimumBoundaryDistance] = ...
        eq9_parameter_boundary_hit(theta, cfg);
    fit = struct('theta',theta,'eta',details.eta, ...
        'covarianceEta',details.covarianceEta,'objective',objectiveValue, ...
        'profileNll',details.profileNll, ...
        'remlAdjustment',details.remlAdjustment,'details',details, ...
        'exitFlag',optimization.exitFlag, ...
        'convergenceCertified',optimization.convergenceCertified, ...
        'elapsedSeconds',optimization.elapsedSeconds, ...
        'optimization',optimization,'boundaryHit',boundaryHit, ...
        'minimumBoundaryDistance',minimumBoundaryDistance, ...
        'selectedJ',selectedJ);
    fit.totalFunctionEvaluations = optimization_evaluations(optimization);
end

function [nll, details] = transformed_objective(z, prepared, cfg, ...
        master, selectedJ, A)
    details = struct();
    try
        theta = A.unpack(z, cfg, [], []);
        [nll, details] = theta_objective(theta, prepared, cfg, master, ...
            selectedJ, A);
    catch exception
        nll = realmax('double') / 1e100;
        details.errorIdentifier = exception.identifier;
        details.errorMessage = exception.message;
    end
end

function [nll, details] = theta_objective(theta, prepared, cfg, ...
        master, selectedJ, A)
    lambdaAll = eq9_fourier_frequencies(prepared.T);
    lambda = lambdaAll(prepared.usedIndices);
    weights = A.weights(theta, lambda, master, cfg);
    p = prepared.nCoefficients;
    n = prepared.nStations;
    H = zeros(p,p); g = zeros(p,1); baseNll = 0;
    maximumJitter = 0; minimumRcond = Inf;
    endpointSet = prepared.endpointIndices(:).';
    for q = 1:numel(prepared.usedIndices)
        originalIndex = prepared.usedIndices(q);
        J = selectedJ(positive_frequency_position(originalIndex, prepared.T));
        M = J^2;
        Phi = prepared.masterPhi(:,1:M);
        factor = A.factor(Phi, weights(1:M,q), theta.noiseVariance, cfg);
        d = prepared.D(originalIndex,:).';
        design = reshape(prepared.Afrequency(originalIndex,:,:),n,p);
        endpoint = ismember(originalIndex, endpointSet);
        if endpoint, d = real(d); design = real(design); end
        solveD = A.solve(factor,d); solveA = A.solve(factor,design);
        if endpoint
            H = H + real(design' * solveA);
            g = g + real(design' * solveD);
            baseNll = baseNll + 0.5 * (factor.logDeterminant + ...
                real(d' * solveD) + n * log(2*pi));
        else
            H = H + 2 * real(design' * solveA);
            g = g + 2 * real(design' * solveD);
            baseNll = baseNll + factor.logDeterminant + ...
                real(d' * solveD) + n * log(pi);
        end
        maximumJitter = max(maximumJitter, factor.jitter);
        minimumRcond = min(minimumRcond, factor.reciprocalCondition);
    end
    H = real((H+H')/2);
    [LH, meanJitter] = eq9_chol_with_jitter(H);
    eta = LH' \ (LH \ g);
    profileNll = baseNll - 0.5 * real(eta' * g);
    remlAdjustment = sum(log(diag(LH)));
    if strcmpi(cfg.likelihood.type,'whittle_reml')
        nll = profileNll + remlAdjustment;
    else
        nll = profileNll;
    end
    details.eta = eta;
    details.covarianceEta = LH' \ (LH \ eye(p));
    details.profileNll = profileNll;
    details.remlAdjustment = remlAdjustment;
    details.maximumFeatureJitter = maximumJitter;
    details.minimumFeatureReciprocalCondition = minimumRcond;
    details.meanInformationJitter = meanJitter;
end

function position = positive_frequency_position(index, T)
    K = floor(T/2)+1;
    if index <= K
        position = index;
    else
        position = T-index+2;
    end
end

function count = optimization_evaluations(o)
    count = numel(o.objectiveAtStart);
    for j=1:numel(o.outputs)
        if isfield(o.outputs{j},'funcCount'), count=count+o.outputs{j}.funcCount; end
    end
    if o.polishTriggered && isfield(o.polishOutput,'funcCount')
        count=count+o.polishOutput.funcCount;
    end
end
