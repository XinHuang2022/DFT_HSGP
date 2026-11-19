function summary = run_ghcnh_spatial_baselines(useParallel,nWorkers)
%RUN_GHCNH_SPATIAL_BASELINES Spatial baselines on the protected GHCNh folds.
%
% The suite reuses the deterministic five-fold whole-site partition from the
% adaptive-HSGP experiment. Hyperparameter selection for the Matern-3/2 RBF
% interpolator is nested within the fitting sites of each outer fold.

    if nargin<1 || isempty(useParallel), useParallel=true; end
    if nargin<2 || isempty(nWorkers), nWorkers=5; end

    suiteRoot=fileparts(mfilename('fullpath'));
    resultDir=fullfile(suiteRoot,'Results','expanded_complete_0839', ...
        'spatial_baselines');
    if ~isfolder(resultDir), mkdir(resultDir); end
    diary(fullfile(resultDir,'progress.log'));
    diaryCleanup=onCleanup(@() diary('off'));
    timerAll=tic;

    addpath(suiteRoot);
    [~,subset,folds,~,E]=ghcnh_expanded_setup();
    U=subset.U; xy=subset.xyKm;
    nSite=size(U,2); nFold=numel(folds);
    assert(nSite==839 && size(U,1)==360, ...
        'Unexpected expanded GHCNh panel dimensions.');

    foldID=zeros(nSite,1);
    for f=1:nFold, foldID(folds(f).targetIndex)=f; end
    assert(all(foldID>0),'Incomplete GHCNh fold assignment.');
    manifest=table(xy(:,1),xy(:,2),foldID, ...
        'VariableNames',{'xKm','yKm','foldID'});

    D=pairwise_distance(xy,xy); D(D==0)=Inf;
    spacingKm=median(min(D,[],2));
    rbfLengthScaleKm=spacingKm.*[0.5,1,2,4,8,16,32];
    rbfLambda=[0,1e-4,1e-3,1e-2,1,10,100];
    krigingRangeKm=4*spacingKm;
    fprintf(['GHCNh spatial baselines started %s\n', ...
        'Panel: 360 x %d; median nearest spacing %.3f km.\n'], ...
        char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
        nSite,spacingKm);

    method=["nearest_neighbour";"idw_k4_p2"; ...
        "matern32_rbf_nested";"universal_kriging_4spacing"];
    displayName=["Nearest neighbour";"IDW (k=4, p=2)"; ...
        "Matern-3/2 RBF (nested CV)"; ...
        "Universal kriging (range four spacings)"];
    prediction=NaN(size(U,1),nSite,numel(method));

    fprintf('[1/4] Nearest-neighbour and IDW interpolation.\n');
    for f=1:nFold
        target=find(foldID==f); train=find(foldID~=f);
        Dt=pairwise_distance(xy(target,:),xy(train,:));
        [sortedDistance,order]=sort(Dt,2,'ascend');
        prediction(:,target,1)=U(:,train(order(:,1)));
        localIndex=order(:,1:4);
        weights=max(sortedDistance(:,1:4),eps).^-2;
        weights=weights./sum(weights,2);
        for j=1:numel(target)
            prediction(:,target(j),2)= ...
                U(:,train(localIndex(j,:)))*weights(j,:).';
        end
        fprintf('      geometric baselines fold %d complete.\n',f);
    end

    if useParallel
        pool=gcp('nocreate');
        if isempty(pool), parpool('local',nWorkers); end
    end
    fprintf(['[2/4] Nested Matern-3/2 RBF and universal kriging ', ...
        '(range %.3f km).\n'],krigingRangeKm);
    rbfFold=cell(nFold,1); krigingFold=cell(nFold,1);
    if useParallel
        parfor f=1:nFold
            rbfFold{f}=fit_rbf_nested_fold( ...
                U,manifest,f,rbfLengthScaleKm,rbfLambda);
            krigingFold{f}=fit_kriging_fold( ...
                U,xy,foldID,f,krigingRangeKm);
            fprintf('      model-based baselines fold %d complete.\n',f);
        end
    else
        for f=1:nFold
            rbfFold{f}=fit_rbf_nested_fold( ...
                U,manifest,f,rbfLengthScaleKm,rbfLambda);
            krigingFold{f}=fit_kriging_fold( ...
                U,xy,foldID,f,krigingRangeKm);
            fprintf('      model-based baselines fold %d complete.\n',f);
        end
    end

    selectedRows=cell(nFold,1);
    for f=1:nFold
        target=rbfFold{f}.target;
        prediction(:,target,3)=rbfFold{f}.prediction;
        prediction(:,krigingFold{f}.target,4)=krigingFold{f}.prediction;
        selectedRows{f}=table(f,rbfFold{f}.lengthScaleKm, ...
            rbfFold{f}.lambda,rbfFold{f}.innerMeanSiteRmse, ...
            rbfFold{f}.innerSeconds,rbfFold{f}.outerSeconds, ...
            'VariableNames',{'foldID','lengthScaleKm','lambda', ...
            'innerMeanSiteRmse','innerSeconds','outerSeconds'});
    end
    assert(all(isfinite(prediction),'all'),'Baseline predictions incomplete.');
    writetable(vertcat(selectedRows{:}),fullfile(resultDir, ...
        'matern32_rbf_selected_hyperparameters.csv'));

    fprintf('[3/4] Metrics and paired spatial-block bootstrap.\n');
    residual=prediction-U;
    nMethod=numel(method);
    siteRmse=squeeze(sqrt(mean(residual.^2,1)));
    siteMae=squeeze(mean(abs(residual),1));
    pooledRmse=squeeze(sqrt(mean(residual.^2,[1,2])));
    pooledMae=squeeze(mean(abs(residual),[1,2]));
    [rmseCI,maeCI,nBlocks]=spatial_block_intervals( ...
        residual,xy,E.spatialBlockKm,E.bootstrapReplicates,20261115);

    summary=table(method,displayName,repmat(nSite,nMethod,1), ...
        pooledRmse,pooledMae,rmseCI(:,1),rmseCI(:,2), ...
        maeCI(:,1),maeCI(:,2),repmat(nBlocks,nMethod,1), ...
        repmat(E.bootstrapReplicates,nMethod,1), ...
        'VariableNames',{'method','displayName','nSites','rmse','mae', ...
        'rmseLower95','rmseUpper95','maeLower95','maeUpper95', ...
        'spatialBlocks','bootstrapReplicates'});
    writetable(summary,fullfile(resultDir,'baseline_summary.csv'));

    siteRows=cell(nMethod,1);
    for q=1:nMethod
        siteRows{q}=table(repmat(method(q),nSite,1), ...
            repmat(displayName(q),nSite,1),foldID, ...
            string(subset.siteTable.stationID),xy(:,1),xy(:,2), ...
            siteRmse(:,q),siteMae(:,q), ...
            'VariableNames',{'method','displayName','foldID','stationID', ...
            'xKm','yKm','rmse','mae'});
    end
    writetable(vertcat(siteRows{:}),fullfile(resultDir, ...
        'baseline_site_metrics.csv'));

    metadata=struct('created',char(datetime('now')), ...
        'panelFile',fullfile(suiteRoot,'data_processed', ...
        'ghcnh_complete_panel_expanded.mat'), ...
        'foldFunction','make_ghcnh_cv_folds', ...
        'medianNearestSpacingKm',spacingKm, ...
        'rbfLengthScaleKm',rbfLengthScaleKm, ...
        'rbfLambda',rbfLambda, ...
        'krigingRangeKm',krigingRangeKm, ...
        'bootstrapReplicates',E.bootstrapReplicates, ...
        'spatialBlockKm',E.spatialBlockKm, ...
        'nWorkers',nWorkers,'elapsedSeconds',toc(timerAll));
    save(fullfile(resultDir,'baseline_predictions.mat'), ...
        'prediction','U','foldID','method','displayName','metadata','-v7.3');
    fid=fopen(fullfile(resultDir,'run_metadata.json'),'w');
    fprintf(fid,'%s',jsonencode(metadata,PrettyPrint=true)); fclose(fid);
    disp(summary);
    fprintf('[4/4] GHCNh spatial baselines completed in %.2f seconds.\n', ...
        metadata.elapsedSeconds);
end

function result=fit_rbf_nested_fold(U,manifest,outerFold,ellGrid,lambdaGrid)
    outerTrain=find(manifest.foldID~=outerFold);
    target=find(manifest.foldID==outerFold);
    innerFolds=setdiff((1:5).',outerFold,'stable');
    [ell,lambda]=ndgrid(ellGrid,lambdaGrid);
    ell=ell(:); lambda=lambda(:); nCandidate=numel(ell);
    score=zeros(numel(innerFolds),nCandidate);
    timer=tic;
    for i=1:numel(innerFolds)
        validation=find(manifest.foldID==innerFolds(i));
        train=find(manifest.foldID~=outerFold & ...
            manifest.foldID~=innerFolds(i));
        truth=U(:,validation);
        for c=1:nCandidate
            fitted=rbf_predict(U(:,train).', ...
                [manifest.xKm(train),manifest.yKm(train)], ...
                [manifest.xKm(validation),manifest.yKm(validation)], ...
                ell(c),lambda(c));
            siteRmse=sqrt(mean((fitted-truth).^2,1));
            score(i,c)=mean(siteRmse);
        end
    end
    innerSeconds=toc(timer);
    innerMean=mean(score,1).';
    [~,selected]=min(innerMean);
    timer=tic;
    fitted=rbf_predict(U(:,outerTrain).', ...
        [manifest.xKm(outerTrain),manifest.yKm(outerTrain)], ...
        [manifest.xKm(target),manifest.yKm(target)], ...
        ell(selected),lambda(selected));
    result=struct('target',target,'prediction',fitted, ...
        'lengthScaleKm',ell(selected),'lambda',lambda(selected), ...
        'innerMeanSiteRmse',innerMean(selected), ...
        'innerSeconds',innerSeconds,'outerSeconds',toc(timer));
end

function prediction=rbf_predict(Ytrain,xyTrain,xyTarget,ell,lambda)
    center=mean(xyTrain,1); scale=std(xyTrain,0,1);
    Xtrain=[ones(size(xyTrain,1),1),(xyTrain-center)./scale];
    Xtarget=[ones(size(xyTarget,1),1),(xyTarget-center)./scale];
    beta=Xtrain\Ytrain; residual=Ytrain-Xtrain*beta;
    K=matern32(pairwise_distance(xyTrain,xyTrain),ell);
    Kcross=matern32(pairwise_distance(xyTarget,xyTrain),ell);
    [L,~]=stable_cholesky(K+lambda.*eye(size(K)));
    weight=L'\(L\residual);
    prediction=(Xtarget*beta+Kcross*weight).';
end

function result=fit_kriging_fold(U,xy,foldID,targetFold,rangeKm)
    target=find(foldID==targetFold); train=find(foldID~=targetFold);
    xyTrain=xy(train,:); xyTarget=xy(target,:);
    center=mean(xyTrain,1); scale=std(xyTrain,0,1);
    Xtrain=[ones(numel(train),1),(xyTrain-center)./scale];
    Xtarget=[ones(numel(target),1),(xyTarget-center)./scale];
    Ytrain=U(:,train).';
    R=matern32(pairwise_distance(xyTrain,xyTrain),rangeKm);
    Kcross=matern32(pairwise_distance(xyTarget,xyTrain),rangeKm);
    [L,jitter]=stable_cholesky(R);
    RinvX=L'\(L\Xtrain); information=Xtrain'*RinvX;
    [LI,~]=stable_cholesky(information);
    RinvY=L'\(L\Ytrain);
    beta=LI'\(LI\(Xtrain'*RinvY));
    residual=Ytrain-Xtrain*beta;
    alpha=L'\(L\residual);
    result=struct('target',target, ...
        'prediction',(Xtarget*beta+Kcross*alpha).', ...
        'jitter',jitter);
end

function [rmseCI,maeCI,nBlocks]=spatial_block_intervals( ...
        residual,xy,blockKm,B,seed)
    bx=floor((xy(:,1)-min(xy(:,1)))/blockKm);
    by=floor((xy(:,2)-min(xy(:,2)))/blockKm);
    [~,~,blockID]=unique([bx,by],'rows','stable');
    nBlocks=max(blockID); nMethod=size(residual,3);
    sumSq=zeros(nBlocks,nMethod); sumAbs=zeros(nBlocks,nMethod);
    count=zeros(nBlocks,1);
    for b=1:nBlocks
        index=find(blockID==b); count(b)=size(residual,1)*numel(index);
        for q=1:nMethod
            value=residual(:,index,q);
            sumSq(b,q)=sum(value.^2,'all');
            sumAbs(b,q)=sum(abs(value),'all');
        end
    end
    rng(seed,'twister'); bootRmse=zeros(B,nMethod); bootMae=zeros(B,nMethod);
    for r=1:B
        draw=randi(nBlocks,nBlocks,1); denominator=sum(count(draw));
        bootRmse(r,:)=sqrt(sum(sumSq(draw,:),1)./denominator);
        bootMae(r,:)=sum(sumAbs(draw,:),1)./denominator;
    end
    rmseCI=[prctile(bootRmse,2.5,1).',prctile(bootRmse,97.5,1).'];
    maeCI=[prctile(bootMae,2.5,1).',prctile(bootMae,97.5,1).'];
end

function value=matern32(distanceKm,rangeKm)
    scaled=sqrt(3).*distanceKm./rangeKm;
    value=(1+scaled).*exp(-scaled);
end

function distance=pairwise_distance(A,B)
    distance=sqrt(max(sum(A.^2,2)+sum(B.^2,2).'-2*A*B.',0));
end

function [L,jitter]=stable_cholesky(A)
    A=(A+A')./2; scale=max(1,mean(abs(diag(A)))); jitter=0;
    for attempt=1:12
        [L,flag]=chol(A+jitter.*eye(size(A)),'lower');
        if flag==0, return; end
        if attempt==1, jitter=1e-12.*scale; else, jitter=10.*jitter; end
    end
    error('Cholesky factorization failed after adaptive jitter.');
end
