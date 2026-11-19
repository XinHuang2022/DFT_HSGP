function basis = subset_adaptive_basis(master, J)
%SUBSET_ADAPTIVE_BASIS Return the first J^2 nested square-shell modes.
    M = J^2;
    if M > master.nFeatures, error('Requested J exceeds master basis.'); end
    basis = master;
    fields = {'indices','omegaByDimension','laplaceEigenvalue', ...
        'radialFrequency'};
    for j = 1:numel(fields), basis.(fields{j}) = basis.(fields{j})(1:M,:); end
    if isfield(basis, 'Phi'), basis.Phi = basis.Phi(:,1:M); end
    basis.basisPerDimension = J;
    basis.nFeatures = M;
end
