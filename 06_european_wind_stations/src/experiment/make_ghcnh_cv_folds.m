function folds = make_ghcnh_cv_folds(subset,nFolds)
%MAKE_GHCNH_CV_FOLDS Deterministic, spatially distributed whole-site folds.

    if nargin<2, nFolds=5; end
    xy=rescale(subset.xyKm,0,1);
    qx=min(255,floor(256*xy(:,1))); qy=min(255,floor(256*xy(:,2)));
    code=morton_code(uint16(qx),uint16(qy));
    [~,order]=sort(code);
    foldID=zeros(size(code));
    foldID(order)=mod((0:numel(order)-1)',nFolds)+1;
    folds=repmat(struct('foldID',0,'targetIndex',[]),nFolds,1);
    for f=1:nFolds
        folds(f).foldID=f;
        folds(f).targetIndex=find(foldID==f);
    end
end

function code=morton_code(x,y)
    code=zeros(numel(x),1,'uint32');
    for bit=0:7
        code=bitor(code,bitshift(uint32(bitget(x,bit+1)),2*bit));
        code=bitor(code,bitshift(uint32(bitget(y,bit+1)),2*bit+1));
    end
end
