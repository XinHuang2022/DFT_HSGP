function result=cerra_spatial_block_bootstrap_difference(a,b,nBlocksPerAxis,B,seed)
%CERRA_SPATIAL_BLOCK_BOOTSTRAP_DIFFERENCE Paired block bootstrap of RMSE.

    [station,ia,ib]=intersect(string(a.stationID),string(b.stationID),'stable');
    if numel(station)<2,error('Matched site tables do not overlap.');end
    a=a(ia,:);b=b(ib,:);difference=a.rmse-b.rmse;
    x=a.longitude;y=a.latitude;
    xPad=sqrt(eps)*max(1,max(x)-min(x));
    yPad=sqrt(eps)*max(1,max(y)-min(y));
    xEdges=linspace(min(x)-xPad,max(x)+xPad,nBlocksPerAxis+1);
    yEdges=linspace(min(y)-yPad,max(y)+yPad,nBlocksPerAxis+1);
    xBin=discretize(x,xEdges);yBin=discretize(y,yEdges);
    block=(yBin-1)*nBlocksPerAxis+xBin;
    levels=unique(block);rng(seed,'twister');samples=zeros(B,1);
    for q=1:B
        selected=levels(randi(numel(levels),numel(levels),1));values=[];
        for j=1:numel(selected)
            values=[values;difference(block==selected(j))]; %#ok<AGROW>
        end
        samples(q)=mean(values);
    end
    samples=sort(samples);ci=samples([max(1,ceil(.025*B)),min(B,ceil(.975*B))]);
    result.meanDifference=mean(difference);
    result.lower95=ci(1);result.upper95=ci(2);
    result.nSites=numel(difference);result.nSpatialBlocks=numel(levels);
end
