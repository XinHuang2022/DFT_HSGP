function [folds, manifest, summary] = make_cerra_cv_folds( ...
        subset, pattern)
%MAKE_CERRA_CV_FOLDS Fixed spatial folds for CERRA reconstruction.
%
% The balanced-five design is a Latin-square coloring.  For a 5-by-5
% subset every fold contains exactly five sites and every row/column has
% one member of every fold.  The same deterministic rule extends to the
% 7-by-7 and 17-by-17 grids.

    pattern = lower(char(pattern));
    sites = subset.siteTable;
    nSites = height(sites);
    foldID = zeros(nSites, 1);

    switch pattern
        case 'none'
            % Leave foldID equal to zero.
        case 'balanced5'
            rowLevels = unique(sites.gridRow, 'sorted');
            columnLevels = unique(sites.gridColumn, 'sorted');
            [rowPresent, rowRank] = ismember(sites.gridRow, rowLevels);
            [columnPresent, columnRank] = ismember( ...
                sites.gridColumn, columnLevels);
            if ~all(rowPresent) || ~all(columnPresent)
                error('Could not rank CERRA grid rows and columns.');
            end
            foldID = 1 + mod((rowRank - 1) + 2 * ...
                (columnRank - 1), 5);
        otherwise
            targetIndex = cerra_holdout_indices(subset, pattern);
            foldID(targetIndex) = (1:numel(targetIndex))';
    end

    manifest = sites;
    manifest.foldID = foldID;
    manifest.isCvTarget = foldID > 0;
    active = unique(foldID(foldID > 0), 'sorted');
    folds = repmat(struct('foldID', 0, 'targetIndex', [], ...
        'targetStations', strings(0, 1)), numel(active), 1);
    rows = cell(numel(active), 1);
    for j = 1:numel(active)
        id = active(j);
        targetIndex = find(foldID == id);
        folds(j).foldID = id;
        folds(j).targetIndex = targetIndex;
        folds(j).targetStations = string(sites.stationID(targetIndex));
        rows{j} = table(id, numel(targetIndex), ...
            nSites - numel(targetIndex), ...
            strjoin(folds(j).targetStations, ';'), ...
            'VariableNames', {'foldID', 'nTargetSites', ...
            'nTrainingSites', 'targetStations'});
    end
    if isempty(rows)
        summary = table(zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
            strings(0, 1), 'VariableNames', {'foldID', ...
            'nTargetSites', 'nTrainingSites', 'targetStations'});
    else
        summary = vertcat(rows{:});
    end

    if strcmp(pattern, 'balanced5')
        if any(foldID == 0) || numel(unique(foldID)) ~= 5
            error('Balanced five-fold allocation did not cover all sites.');
        end
        counts = accumarray(foldID, 1, [5, 1]);
        if max(counts) - min(counts) > 1
            error('Balanced five-fold target counts differ by more than one.');
        end
    end
end


