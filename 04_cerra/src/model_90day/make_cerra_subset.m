function subset = make_cerra_subset(contract, name)
%MAKE_CERRA_SUBSET Deterministic square subsets of the fixed CERRA grid.

    if nargin < 2, name = 'full17'; end
    switch lower(char(name))
        case 'full17', keep = true(height(contract.siteTable), 1);
        case 'central9'
            keep = ismember(contract.siteTable.gridRow, 5:13) & ...
                ismember(contract.siteTable.gridColumn, 5:13);
        case 'central5'
            keep = ismember(contract.siteTable.gridRow, 7:11) & ...
                ismember(contract.siteTable.gridColumn, 7:11);
        otherwise, error('Unknown CERRA subset: %s', name);
    end
    subset.name = string(name);
    subset.U = contract.U(:, keep);
    subset.siteTable = contract.siteTable(keep, :);
    subset.xyKm = contract.xyKm(keep, :);
    subset.timeUtc = contract.timeUtc;
    subset.timeIndex = contract.timeIndex;
    subset.relativeHours = contract.relativeHours;
end

