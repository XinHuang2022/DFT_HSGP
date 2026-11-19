function result = eq9_multistart_optimize(objective, baseZ, optimizer, seed)
%EQ9_MULTISTART_OPTIMIZE Reproducible Nelder-Mead with selective polishing.
%
% The ordinary multistart search is followed by one longer restart from the
% best solution only when its exit flag is nonpositive.  Polishing therefore
% targets incomplete termination without multiplying the cost of every fit.

    baseZ = baseZ(:);
    nStarts = optimizer.numStarts;
    starts = zeros(numel(baseZ), nStarts);
    starts(:, 1) = baseZ;
    additionalStarts = field_or_default(optimizer, ...
        'additionalStartZ', zeros(numel(baseZ), 0));
    if isempty(additionalStarts)
        additionalStarts = zeros(numel(baseZ), 0);
    end
    if size(additionalStarts, 1) ~= numel(baseZ)
        error('optimizer.additionalStartZ has an incompatible row count.');
    end
    nAdditional = min(size(additionalStarts, 2), nStarts - 1);
    if nAdditional > 0
        starts(:, 2:(1 + nAdditional)) = ...
            additionalStarts(:, 1:nAdditional);
    end
    rng(seed, 'twister');
    for j = (2 + nAdditional):nStarts
        starts(:, j) = baseZ + optimizer.randomStartScale * ...
            randn(size(baseZ));
    end

    options = optimset('Display', optimizer.display, ...
        'TolX', optimizer.tolX, 'TolFun', optimizer.tolFun, ...
        'MaxIter', optimizer.maxIterations, ...
        'MaxFunEvals', optimizer.maxFunctionEvaluations);
    objectiveAtStart = zeros(nStarts, 1);
    objectiveAtSolution = zeros(nStarts, 1);
    exitFlags = zeros(nStarts, 1);
    solutions = zeros(numel(baseZ), nStarts);
    outputs = cell(nStarts, 1);

    timerValue = tic;
    for j = 1:nStarts
        objectiveAtStart(j) = objective(starts(:, j));
        [solutions(:, j), objectiveAtSolution(j), exitFlags(j), ...
            outputs{j}] = fminsearch(objective, starts(:, j), options);
    end
    [bestObjective, bestIndex] = min(objectiveAtSolution);
    bestZ = solutions(:, bestIndex);
    bestExitFlag = exitFlags(bestIndex);
    bestOutput = outputs{bestIndex};

    polishEnabled = field_or_default(optimizer, ...
        'selectivePolishingEnabled', false);
    polishTriggered = polishEnabled && bestExitFlag <= 0;
    polishAccepted = false;
    polishExitFlag = NaN;
    polishObjective = NaN;
    polishOutput = struct();
    prePolishObjective = bestObjective;
    prePolishExitFlag = bestExitFlag;

    if polishTriggered
        polishOptions = optimset('Display', optimizer.display, ...
            'TolX', field_or_default(optimizer, 'polishTolX', ...
                optimizer.tolX), ...
            'TolFun', field_or_default(optimizer, 'polishTolFun', ...
                optimizer.tolFun), ...
            'MaxIter', field_or_default(optimizer, ...
                'polishMaxIterations', 1000), ...
            'MaxFunEvals', field_or_default(optimizer, ...
                'polishMaxFunctionEvaluations', 5000));
        [polishZ, polishObjective, polishExitFlag, polishOutput] = ...
            fminsearch(objective, bestZ, polishOptions);
        acceptTolerance = field_or_default(optimizer, ...
            'polishAcceptObjectiveTolerance', 1e-4);
        polishAccepted = polishObjective <= bestObjective || ...
            (polishExitFlag > 0 && ...
             polishObjective <= bestObjective + acceptTolerance);
        if polishAccepted
            bestZ = polishZ;
            bestObjective = polishObjective;
            bestExitFlag = polishExitFlag;
            bestOutput = polishOutput;
        end
    end

    convergenceCertified = bestExitFlag > 0;
    if ~convergenceCertified && polishTriggered && polishExitFlag > 0
        acceptTolerance = field_or_default(optimizer, ...
            'polishAcceptObjectiveTolerance', 1e-4);
        convergenceCertified = polishObjective <= ...
            prePolishObjective + acceptTolerance;
    end

    sortedObjectives = sort(objectiveAtSolution(isfinite( ...
        objectiveAtSolution)));
    if numel(sortedObjectives) >= 2
        secondBestObjectiveGap = sortedObjectives(2) - ...
            sortedObjectives(1);
    else
        secondBestObjectiveGap = NaN;
    end
    nearBestTolerance = field_or_default(optimizer, ...
        'nearBestObjectiveTolerance', 1e-3);

    result.z = bestZ;
    result.objective = bestObjective;
    result.exitFlag = bestExitFlag;
    result.output = bestOutput;
    result.bestStartIndex = bestIndex;
    result.starts = starts;
    result.solutions = solutions;
    result.objectiveAtStart = objectiveAtStart;
    result.objectiveAtSolution = objectiveAtSolution;
    result.exitFlags = exitFlags;
    result.outputs = outputs;
    result.elapsedSeconds = toc(timerValue);
    result.improvementFromOwnStart = ...
        objectiveAtStart(bestIndex) - bestObjective;
    result.secondBestObjectiveGap = secondBestObjectiveGap;
    result.nearBestStartCount = sum(objectiveAtSolution <= ...
        min(objectiveAtSolution) + nearBestTolerance);
    result.polishTriggered = polishTriggered;
    result.polishAccepted = polishAccepted;
    result.prePolishObjective = prePolishObjective;
    result.prePolishExitFlag = prePolishExitFlag;
    result.polishObjective = polishObjective;
    result.polishExitFlag = polishExitFlag;
    result.polishOutput = polishOutput;
    result.convergenceCertified = convergenceCertified;
end

function value = field_or_default(structure, name, defaultValue)
    if isfield(structure, name) && ~isempty(structure.(name))
        value = structure.(name);
    else
        value = defaultValue;
    end
end
