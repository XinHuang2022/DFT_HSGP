function distances = eq9_pairwise_distances(xyA, xyB)
%EQ9_PAIRWISE_DISTANCES Euclidean distance matrix without toolboxes.

    if nargin < 2 || isempty(xyB)
        xyB = xyA;
    end
    if size(xyA, 2) ~= size(xyB, 2)
        error('Coordinate arrays must have the same number of columns.');
    end
    a2 = sum(xyA.^2, 2);
    b2 = sum(xyB.^2, 2)';
    squared = max(0, a2 + b2 - 2 * (xyA * xyB'));
    distances = sqrt(squared);
end



