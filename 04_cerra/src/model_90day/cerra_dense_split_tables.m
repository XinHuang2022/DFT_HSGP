function [siteRoles,coarseMembership]=cerra_dense_split_tables(contract,split)
%CERRA_DENSE_SPLIT_TABLES Human-readable manifests for the fixed geometry.

    n=height(contract.siteTable);
    role=repmat("full_training_only",n,1);
    role(split.testIndex)="held_out_test";
    siteRoles=contract.siteTable;
    siteRoles.role=role;
    sizes=split.exactSizes(:);
    rows=cell(numel(sizes),1);
    for j=1:numel(sizes)
        index=split.coarseOrder(1:sizes(j));
        rows{j}=table(repmat(sizes(j),numel(index),1), ...
            (1:numel(index))',index,string(contract.siteTable.stationID(index)), ...
            contract.siteTable.latitude(index),contract.siteTable.longitude(index), ...
            'VariableNames',{'nCoarseSites','selectionOrder','fullSiteIndex', ...
            'stationID','latitude','longitude'});
    end
    coarseMembership=vertcat(rows{:});
end
