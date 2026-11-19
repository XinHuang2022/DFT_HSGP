function split=make_cerra_dense_split(contract,cfg)
%MAKE_CERRA_DENSE_SPLIT Fixed test set and nested space-filling baselines.

    nSites=height(contract.siteTable);
    if cfg.split.nTestSites+max(cfg.split.exactSizes)>nSites
        error('Requested test and exact subsets exceed the available sites.');
    end
    testOrder=farthest_point_order(contract.xyKm,cfg.split.nTestSites, ...
        cfg.seed.split);
    testIndex=sort(testOrder(:));
    fullTrainingIndex=setdiff((1:nSites)',testIndex,'stable');
    coarseLocalOrder=farthest_point_order(contract.xyKm(fullTrainingIndex,:), ...
        max(cfg.split.exactSizes),cfg.seed.split+1);
    coarseOrder=fullTrainingIndex(coarseLocalOrder);
    split.version=string(cfg.version);
    split.design=string(cfg.split.design);
    split.testIndex=testIndex;
    split.fullTrainingIndex=fullTrainingIndex;
    split.coarseOrder=coarseOrder;
    split.exactSizes=cfg.split.exactSizes(:).';
end

function order=farthest_point_order(xy,nSelect,seed)
    n=size(xy,1);
    if nSelect>n,error('Cannot select more points than are available.');end
    center=mean(xy,1);
    centerDistance=sum((xy-center).^2,2);
    [~,ranked]=sortrows([centerDistance,(1:n)'],[1 2]);
    first=ranked(1+mod(seed-1,min(n,11)));
    order=zeros(nSelect,1);order(1)=first;
    minimumDistance=sum((xy-xy(first,:)).^2,2);
    minimumDistance(first)=-Inf;
    for j=2:nSelect
        maximum=max(minimumDistance);
        candidates=find(minimumDistance==maximum);
        chosen=candidates(1+mod(seed+j-2,numel(candidates)));
        order(j)=chosen;
        distance=sum((xy-xy(chosen,:)).^2,2);
        minimumDistance=min(minimumDistance,distance);
        minimumDistance(order(1:j))=-Inf;
    end
end
