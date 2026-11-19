function D = design(C)
%DESIGN Jittered grid and one fixed whole-site train/test split.
    rng(C.seed,'twister');
    u=linspace(0.04,0.96,C.nSide);
    [x,y]=meshgrid(u,u);
    xy=[x(:),y(:)];
    if C.nSide>4
        step=0.12/(C.nSide-1);
        xy=xy+step*(rand(size(xy))-.5);
        xy=min(.98,max(.02,xy));
    end
    % One deterministic site from each stratum, followed by a seeded fill.
    order=randperm(size(xy,1));
    test=sort(order(1:C.nTest));
    train=setdiff((1:size(xy,1))',test(:),'stable');
    D.xy=xy; D.train=train; D.test=test;
    D.site=(1:size(xy,1))';
end

