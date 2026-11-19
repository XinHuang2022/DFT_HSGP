function summary=bootstrap_ghcnh_expanded_cv(methods)
%BOOTSTRAP_GHCNH_EXPANDED_CV Paired spatial-block intervals for OOF metrics.

    [~,subset,~,~,E]=ghcnh_expanded_setup();
    if nargin<1 || isempty(methods)
        methods=["exact_kl","trace_bound"];
    end
    methods=string(methods(:))';
    nMethods=numel(methods); predictions=cell(nMethods,1);
    for q=1:nMethods
        parts=cell(5,1);
        for f=1:5
            filePath=fullfile(E.resultsDir,sprintf('fold_%02d',f), ...
                char(methods(q)),'predictions.csv');
            if ~isfile(filePath)
                error('Missing predictions: %s',filePath);
            end
            parts{f}=readtable(filePath,'TextType','string');
        end
        predictions{q}=vertcat(parts{:});
    end

    site=string(subset.siteTable.stationID);
    x=subset.xyKm(:,1); y=subset.xyKm(:,2);
    bx=floor((x-min(x))/E.spatialBlockKm);
    by=floor((y-min(y))/E.spatialBlockKm);
    [~,~,blockID]=unique([bx,by],'rows','stable');
    nBlocks=max(blockID);
    blockMap=table(site,blockID,'VariableNames',{'stationID','blockID'});

    blockStats=cell(nMethods,1); point=zeros(nMethods,4);
    for q=1:nMethods
        P=outerjoin(predictions{q},blockMap,'Keys','stationID', ...
            'MergeKeys',true,'Type','left');
        residual=P.predictionWindMps-P.truthWindMps;
        variance=max(P.observationVariance,1e-10);
        P.sse=residual.^2; P.sae=abs(residual);
        P.nlpdContribution=0.5*(log(2*pi*variance)+P.sse./variance);
        P.covered=double(abs(residual)<=1.95996398454005*sqrt(variance));
        P.one=ones(height(P),1);
        B=groupsummary(P,'blockID','sum', ...
            {'sse','sae','nlpdContribution','covered','one'});
        B=sortrows(B,'blockID'); blockStats{q}=B;
        point(q,:)=[sqrt(sum(B.sum_sse)/sum(B.sum_one)), ...
            sum(B.sum_sae)/sum(B.sum_one), ...
            sum(B.sum_nlpdContribution)/sum(B.sum_one), ...
            sum(B.sum_covered)/sum(B.sum_one)];
    end

    rng(20261113,'twister'); R=E.bootstrapReplicates;
    boot=zeros(R,nMethods,4);
    for r=1:R
        draw=randi(nBlocks,nBlocks,1);
        for q=1:nMethods
            B=blockStats{q}; denominator=sum(B.sum_one(draw));
            boot(r,q,1)=sqrt(sum(B.sum_sse(draw))/denominator);
            boot(r,q,2)=sum(B.sum_sae(draw))/denominator;
            boot(r,q,3)=sum(B.sum_nlpdContribution(draw))/denominator;
            boot(r,q,4)=sum(B.sum_covered(draw))/denominator;
        end
    end

    metricNames=["rmse","mae","nlpd","coverage95"];
    rows=cell(nMethods*numel(metricNames),1); cursor=0;
    for q=1:nMethods
        for j=1:numel(metricNames)
            cursor=cursor+1; values=boot(:,q,j);
            rows{cursor}=table(methods(q),metricNames(j),point(q,j), ...
                prctile(values,2.5),prctile(values,97.5),R,nBlocks, ...
                E.spatialBlockKm,'VariableNames',{'method','metric', ...
                'estimate','ciLower','ciUpper','bootstrapReplicates', ...
                'spatialBlocks','blockSizeKm'});
        end
    end
    summary=vertcat(rows{:});
    writetable(summary,fullfile(E.resultsDir,'fivefold_bootstrap_summary.csv'));

    paired=table();
    if nMethods==2
        difference=zeros(R,4);
        for j=1:4, difference(:,j)=boot(:,2,j)-boot(:,1,j); end
        paired=table(metricNames',point(2,:)'-point(1,:)', ...
            prctile(difference,2.5,1)',prctile(difference,97.5,1)', ...
            'VariableNames',{'metric','secondMinusFirst','ciLower','ciUpper'});
        writetable(paired,fullfile(E.resultsDir, ...
            'fivefold_paired_differences.csv'));
    end
    save(fullfile(E.resultsDir,'fivefold_bootstrap.mat'), ...
        'summary','paired','boot','point','blockID','E','-v7.3');
end
