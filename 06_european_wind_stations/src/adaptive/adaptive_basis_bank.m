function bank = adaptive_basis_bank(domainXy, boundaryFactor, Jmax, A, cfg)
%ADAPTIVE_BASIS_BANK Create a square-shell-ordered nested basis bank.

    if A.prefix == "era5"
        base = make_era5_hsgp_basis(domainXy, boundaryFactor, Jmax, cfg);
    else
        base = make_pacific_hsgp_basis(domainXy, boundaryFactor, Jmax, cfg);
    end
    shell = max(base.indices, [], 2);
    orderTable = table(shell, base.laplaceEigenvalue, ...
        base.indices(:,1), base.indices(:,2), (1:base.nFeatures)', ...
        'VariableNames', {'shell','eigenvalue','j1','j2','original'});
    orderTable = sortrows(orderTable, ...
        {'shell','eigenvalue','j1','j2'});
    order = orderTable.original;
    fields = {'indices','omegaByDimension','laplaceEigenvalue', ...
        'radialFrequency'};
    for j = 1:numel(fields), base.(fields{j}) = base.(fields{j})(order,:); end
    if isfield(base, 'Phi'), base.Phi = base.Phi(:, order); end
    base.modeOrdering = 'nested_square_shells';
    bank = base;
end

