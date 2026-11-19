function split = make_pacific_inner_split(siteTable, parameterIndex, ...
        nFolds, validationFold, seed)
%MAKE_PACIFIC_INNER_SPLIT Coordinate-only spatially dispersed validation.

% A shifted Latin coloring on the original grid supplies an inner
% validation set.  Because the split is made after removing the outer CV
% targets, no outer-target response can select the nugget.

    parameterIndex = parameterIndex(:);
    if validationFold < 1 || validationFold > nFolds
        error('validationFold must lie between 1 and nFolds.');
    end
    row = siteTable.gridRow(parameterIndex);
    column = siteTable.gridColumn(parameterIndex);
    shift = mod(seed, nFolds);
    color = 1 + mod((row - 1) + 2 * (column - 1) + shift, nFolds);
    isValidation = color == validationFold;
    if sum(isValidation) < 3 || sum(~isValidation) < 8
        % Deterministic fallback for an unusually irregular candidate set.
        [~, order] = sortrows([row, column], [1, 2]);
        isValidation = false(size(parameterIndex));
        isValidation(order(validationFold:nFolds:end)) = true;
    end
    split.fitIndex = parameterIndex(~isValidation);
    split.validationIndex = parameterIndex(isValidation);
    split.table = siteTable(parameterIndex, :);
    split.table.innerColor = color;
    split.table.innerRole = repmat("fit", numel(parameterIndex), 1);
    split.table.innerRole(isValidation) = "validation";
end
