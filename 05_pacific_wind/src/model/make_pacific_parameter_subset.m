function selection = make_pacific_parameter_subset(siteTable, ...
        candidateIndex, nSubset, pairFraction, seed)
%MAKE_PACIFIC_PARAMETER_SUBSET Coverage-plus-pairs training-only design.

% Four coordinate-defined boundary anchors and a farthest-point coverage
% design identify long spatial scales.  A prespecified fraction of nearest
% neighbors of those anchors supplies short-distance pairs needed to
% distinguish process variance from a nugget.  Only coordinates and the
% outer-training membership enter the selection.

    candidateIndex = unique(candidateIndex(:), 'stable');
    if nSubset > numel(candidateIndex)
        error('Requested %d parameter sites but only %d are available.', ...
            nSubset, numel(candidateIndex));
    end
    if nSubset < 8
        error('At least eight sites are required for the balanced subset.');
    end
    if pairFraction < 0 || pairFraction >= 0.5
        error('pairFraction must lie in [0, 0.5).');
    end
    xy = [siteTable.xKm(candidateIndex), siteTable.yKm(candidateIndex)];
    rng(seed, 'twister');
    selectedLocal = zeros(nSubset, 1);
    roles = strings(nSubset, 1);
    cursor = 0;

    corners = [min(xy(:, 1)), min(xy(:, 2)); ...
               min(xy(:, 1)), max(xy(:, 2)); ...
               max(xy(:, 1)), min(xy(:, 2)); ...
               max(xy(:, 1)), max(xy(:, 2))];
    for j = 1:4
        distance = sum((xy - corners(j, :)) .^ 2, 2);
        distance(selectedLocal(1:cursor)) = Inf;
        [~, chosen] = min(distance);
        cursor = cursor + 1;
        selectedLocal(cursor) = chosen;
        roles(cursor) = "boundary_anchor";
    end

    nPairs = round(pairFraction * nSubset);
    nCoverage = max(4, nSubset - nPairs);
    while cursor < nCoverage
        distance = eq9_pairwise_distances(xy, xy(selectedLocal(1:cursor), :));
        nearest = min(distance, [], 2);
        nearest(selectedLocal(1:cursor)) = -Inf;
        maximum = max(nearest);
        tied = find(abs(nearest - maximum) <= ...
            1e-12 * max(1, abs(maximum)));
        chosen = tied(randi(numel(tied)));
        cursor = cursor + 1;
        selectedLocal(cursor) = chosen;
        roles(cursor) = "coverage_anchor";
    end

    anchorOrder = randperm(nCoverage);
    anchorCursor = 0;
    while cursor < nSubset
        anchorCursor = anchorCursor + 1;
        if anchorCursor > numel(anchorOrder)
            anchorOrder = randperm(nCoverage);
            anchorCursor = 1;
        end
        anchor = selectedLocal(anchorOrder(anchorCursor));
        distance = sqrt(sum((xy - xy(anchor, :)) .^ 2, 2));
        distance(selectedLocal(1:cursor)) = Inf;
        [minimum, chosen] = min(distance);
        if ~isfinite(minimum)
            break;
        end
        cursor = cursor + 1;
        selectedLocal(cursor) = chosen;
        roles(cursor) = "local_pair";
    end
    if cursor ~= nSubset || numel(unique(selectedLocal)) ~= nSubset
        error('Could not construct the requested parameter subset.');
    end

    globalIndex = candidateIndex(selectedLocal);
    selection.index = globalIndex;
    selection.localCandidateIndex = selectedLocal;
    selection.table = siteTable(globalIndex, :);
    selection.table.selectionOrder = (1:nSubset)';
    selection.table.selectionRole = roles;
    selection.table.selectionSeed = repmat(seed, nSubset, 1);
    selection.table.isOuterTrainingSite = true(nSubset, 1);
end
