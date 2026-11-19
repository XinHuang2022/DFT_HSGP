function fit=fit_model(Y,xy,C,R,engine,B,seed)
%FIT_MODEL Reproducible serial multistart fit for one parallel task.
    base=to_z(C.initial,C); rng(seed,'twister');
    starts=zeros(numel(base),C.starts); starts(:,1)=base;
    for j=2:C.starts, starts(:,j)=base+C.opt.startScale*randn(size(base)); end
    opts=optimset('Display','off','TolX',C.opt.tolX,'TolFun',C.opt.tolFun, ...
        'MaxIter',C.opt.maxIter,'MaxFunEvals',C.opt.maxEval);
    cache=[];
    if strcmp(engine,'hsgp')
        cache.P=features(xy,B);
        cache.G=cache.P'*cache.P;
    end
    obj=@(z)nll_model(z,Y,xy,C,R,engine,B,cache);
    sol=zeros(size(starts)); val=inf(C.starts,1); flag=zeros(C.starts,1);
    out=cell(C.starts,1); t=tic;
    for j=1:C.starts
        [sol(:,j),val(j),flag(j),out{j}]=fminsearch(obj,starts(:,j),opts);
    end
    [fit.nll,b]=min(val); fit.z=sol(:,b); fit.theta=from_z(fit.z,C);
    fit.exitFlag=flag(b); fit.output=out{b}; fit.seconds=toc(t);
    fit.starts=starts; fit.solutions=sol; fit.values=val; fit.flags=flag;
    fit.engine=engine; fit.boundaryDistance=min((fit.theta-C.lower)./ ...
        (C.upper-C.lower),1-(fit.theta-C.lower)./(C.upper-C.lower));
end
