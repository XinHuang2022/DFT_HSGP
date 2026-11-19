function subset = make_pacific_subset(contract, subsetName)
%MAKE_PACIFIC_SUBSET Create documented nested Pacific grid subsets.

    subsetName = lower(char(subsetName));
    switch subsetName
        case 'local3'
            rows = 8:10;
            columns = 8:10;
            description = 'contiguous central 3-by-3 implementation subset';
        case 'local5'
            rows = 7:11;
            columns = 7:11;
            description = 'contiguous central 5-by-5 scientific pilot';
        case 'local6'
            rows = 6:11;
            columns = 6:11;
            description = 'contiguous 6-by-6 scaling pilot';
        case 'local7'
            rows = 6:12;
            columns = 6:12;
            description = 'contiguous central 7-by-7 scaling experiment';
        case 'full17'
            rows = 1:17;
            columns = 1:17;
            description = 'complete 17-by-17 Pacific grid';
        case 'sparse3'
            rows = [1, 9, 17];
            columns = [1, 9, 17];
            description = 'full-domain sparse 3-by-3 stress test';
        otherwise
            error('Unknown Pacific subset %s.', subsetName);
    end

    sites = contract.siteTable;
    selected = ismember(sites.gridRow, rows) & ...
        ismember(sites.gridColumn, columns);
    selectedSites = sortrows(sites(selected, :), ...
        {'gridRow', 'gridColumn'});
    expectedCount = numel(rows) * numel(columns);
    if height(selectedSites) ~= expectedCount
        error('Subset %s contains %d sites; expected %d.', ...
            subsetName, height(selectedSites), expectedCount);
    end
    selectedSites.subsetIndex = (1:height(selectedSites))';
    selectedSites.subsetName = repmat(string(subsetName), ...
        height(selectedSites), 1);
    selectedSites.description = repmat(string(description), ...
        height(selectedSites), 1);
    originalColumns = selectedSites.originalColumn;

    distance = eq9_pairwise_distances( ...
        [selectedSites.xKm, selectedSites.yKm], ...
        [selectedSites.xKm, selectedSites.yKm]);
    distance(1:(height(selectedSites) + 1):end) = Inf;
    nearestDistanceKm = min(distance, [], 2);
    selectedSites.nearestSubsetNeighborKm = nearestDistanceKm;

    subset.name = string(subsetName);
    subset.description = string(description);
    subset.rows = rows(:);
    subset.columns = columns(:);
    subset.siteTable = selectedSites;
    subset.U = contract.U(:, originalColumns);
    subset.xyKm = [selectedSites.xKm, selectedSites.yKm];
    subset.timeIndex = contract.timeIndex;
    subset.relativeHours = contract.relativeHours;
end
