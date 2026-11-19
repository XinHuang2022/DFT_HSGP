function summary=summarize_ghcnh_relaxed_cycle_cv()
%SUMMARIZE_GHCNH_RELAXED_CYCLE_CV Five-fold block-bootstrap summaries.

    [~,subset,~,~,E]=ghcnh_expanded_setup();
    methods=["exact_kl","solve_trace"];
    tags=["maxcycles_4_fraction_0100","maxcycles_3_fraction_0500"];
    displayNames=["frequency_adaptive_algorithm_1_relaxed", ...
        "frequency_adaptive_algorithm_2_relaxed"];
    predictions=cell(2,1);
    foldMetrics=cell(10,1); cursor=0;
    for q=1:2
        parts=cell(5,1);
        for f=1:5
            methodDir=fullfile(E.resultsDir,'relaxed_cycle_sensitivity', ...
                sprintf('fold_%02d',f),char(tags(q)),char(methods(q)));
            predictionFile=fullfile(methodDir,'predictions.csv');
            metricFile=fullfile(methodDir,'metrics.csv');
            if ~isfile(predictionFile) || ~isfile(metricFile)
                error('Missing relaxed-cycle output under %s.',methodDir);
            end
            parts{f}=readtable(predictionFile,'TextType','string');
            metric=readtable(metricFile,'TextType','string');
            metric.method(:)=displayNames(q); cursor=cursor+1;
            foldMetrics{cursor}=metric;
        end
        predictions{q}=vertcat(parts{:});
    end
    foldMetrics=vertcat(foldMetrics{:});

    stationID=string(subset.siteTable.stationID);
    x=subset.xyKm(:,1); y=subset.xyKm(:,2);
    bx=floor((x-min(x))/E.spatialBlockKm);
    by=floor((y-min(y))/E.spatialBlockKm);
    [~,~,blockID]=unique([bx,by],'rows','stable');
    nBlocks=max(blockID);
    blockMap=table(stationID,blockID);

    nMethods=2; blockStats=cell(nMethods,1); point=zeros(nMethods,4);
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
            rows{cursor}=table(displayNames(q),metricNames(j),point(q,j), ...
                prctile(values,2.5),prctile(values,97.5),R,nBlocks, ...
                E.spatialBlockKm,'VariableNames',{'method','metric', ...
                'estimate','ciLower','ciUpper','bootstrapReplicates', ...
                'spatialBlocks','blockSizeKm'});
        end
    end
    summary=vertcat(rows{:});
    outputDir=fullfile(E.resultsDir,'relaxed_cycle_sensitivity');
    writetable(summary,fullfile(outputDir,'relaxed_fivefold_summary.csv'));
    writetable(foldMetrics,fullfile(outputDir,'relaxed_fold_metrics.csv'));
    save(fullfile(outputDir,'relaxed_fivefold_bootstrap.mat'), ...
        'summary','boot','point','blockID','foldMetrics','-v7.3');
    disp(summary);
end
