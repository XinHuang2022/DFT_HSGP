function holdoutIndex = pacific_holdout_indices(subset, pattern)
%PACIFIC_HOLDOUT_INDICES Select documented holdouts within a subset.

    pattern = lower(char(pattern));
    sites = subset.siteTable;
    rowLevels = unique(sites.gridRow, 'sorted');
    columnLevels = unique(sites.gridColumn, 'sorted');
    middleRow = rowLevels(ceil(numel(rowLevels) / 2));
    middleColumn = columnLevels(ceil(numel(columnLevels) / 2));
    switch pattern
        case 'none'
            holdoutIndex = zeros(0, 1);
        case 'center'
            requested = [middleRow, middleColumn];
        case 'representative5'
            requested = [middleRow, middleColumn; ...
                rowLevels(1), middleColumn; rowLevels(end), middleColumn; ...
                middleRow, columnLevels(1); middleRow, columnLevels(end)];
        case 'all'
            holdoutIndex = (1:height(sites))';
            return;
        otherwise
            error('Unknown holdout pattern %s.', pattern);
    end
    holdoutIndex = zeros(size(requested, 1), 1);
    for j = 1:size(requested, 1)
        match = find(sites.gridRow == requested(j, 1) & ...
            sites.gridColumn == requested(j, 2));
        if numel(match) ~= 1
            error('Requested holdout grid point is absent or duplicated.');
        end
        holdoutIndex(j) = match;
    end
    holdoutIndex = unique(holdoutIndex, 'stable');
end
