function [selectedJ,audit,summary] = solve_trace_rank_audit(theta,xy,cfg, ...
        master,C,A,logFile,label,frequencyIndices)
%SOLVE_TRACE_RANK_AUDIT Algorithm-2 rank selection by a solve-aware bound.
%
% For P=Q+R, the normalized Gaussian KL discrepancy satisfies
%   e <= c/n * {t-log(1+t)},  t=trace(Q\R).
% The calculation uses the candidate covariance Q and the positive-semidefinite
% omitted basis tail R. It does not factorize the high-rank reference P.

    T=cfg.data.expectedTimePoints; K=floor(T/2)+1;
    if nargin<9 || isempty(frequencyIndices), frequencyIndices=1:K; end
    frequencyIndices=unique(frequencyIndices(:)');
    if any(frequencyIndices<1 | frequencyIndices>K)
        error('frequencyIndices must lie between 1 and %d.',K);
    end
    lambdaAll=eq9_fourier_frequencies(T); lambda=lambdaAll(1:K);
    Phi=A.features(xy,master);
    weights=A.weights(theta,lambda,master,cfg);
    n=size(Phi,1); referenceM=C.referenceJ^2;
    candidateJ=C.candidateJ(C.candidateJ<C.referenceJ);
    noiseVariance=theta.noiseVariance;
    referenceJ=C.referenceJ;
    epsilonKlPerSite=C.epsilonKlPerSite;
    nCandidate=numel(candidateJ); nFrequency=numel(frequencyIndices);
    selectedSubset=repmat(referenceJ,nFrequency,1);
    values=nan(nFrequency,nCandidate); passes=false(nFrequency,nCandidate);
    evaluated=false(nFrequency,nCandidate);
    elapsedByFrequency=zeros(nFrequency,1);

    outerTask=getCurrentTask();
    useParallel=isempty(outerTask);
    if useParallel
        pool=gcp('nocreate');
        if isempty(pool) || pool.NumWorkers~=C.downloadWorkers
            if ~isempty(pool), delete(pool); end
            pool=parpool('local',C.downloadWorkers);
        end
        workerCount=pool.NumWorkers;
        queue=parallel.pool.DataQueue;
        afterEach(queue,@report_progress);
    else
        workerCount=1;
    end
    completed=0; auditTimer=tic;
    ghcnh_log(logFile,label,sprintf( ...
        ['solve-aware trace audit starting: %d frequencies, at most %d ' ...
         'candidates, n=%d, workers=%d'], ...
        nFrequency,nCandidate,n,workerCount));

    if useParallel
        parfor r=1:nFrequency
            k=frequencyIndices(r);
            [chosen,localValues,localPass,localEvaluated,elapsed]= ...
                audit_frequency(k,K,Phi,weights,noiseVariance,candidateJ, ...
                referenceM,referenceJ,epsilonKlPerSite,n);
            selectedSubset(r)=chosen;
            values(r,:)=localValues; passes(r,:)=localPass;
            evaluated(r,:)=localEvaluated; elapsedByFrequency(r)=elapsed;
            send(queue,[r,k,chosen,sum(localEvaluated),elapsed]);
        end
    else
        for r=1:nFrequency
            k=frequencyIndices(r);
            [chosen,localValues,localPass,localEvaluated,elapsed]= ...
                audit_frequency(k,K,Phi,weights,noiseVariance,candidateJ, ...
                referenceM,referenceJ,epsilonKlPerSite,n);
            selectedSubset(r)=chosen;
            values(r,:)=localValues; passes(r,:)=localPass;
            evaluated(r,:)=localEvaluated; elapsedByFrequency(r)=elapsed;
            report_progress([r,k,chosen,sum(localEvaluated),elapsed]);
        end
    end

    selectedJ=nan(K,1); selectedJ(frequencyIndices)=selectedSubset;
    rows=cell(nnz(evaluated),1); cursor=0;
    for r=1:nFrequency
        k=frequencyIndices(r);
        for q=find(evaluated(r,:))
            cursor=cursor+1; J=candidateJ(q);
            rows{cursor}=table(string(label),k-1,lambda(k),J,J^2, ...
                referenceJ,referenceM,values(r,q),passes(r,q), ...
                'VariableNames',{'auditLabel','frequencyIndex','lambda', ...
                'candidateJ','candidateM','referenceJ','referenceM', ...
                'solveTraceBoundPerSite','passes'});
        end
    end
    if cursor>0, audit=vertcat(rows{1:cursor}); else, audit=table(); end
    totalSeconds=toc(auditTimer);
    summary=table(nFrequency,nCandidate,n,workerCount,totalSeconds, ...
        mean(elapsedByFrequency),median(elapsedByFrequency), ...
        max(elapsedByFrequency),sum(evaluated,'all'), ...
        sum(selectedSubset==referenceJ), ...
        'VariableNames',{'nFrequencies','maximumCandidates', ...
        'nTraining','nWorkers','totalSeconds','meanFrequencySeconds', ...
        'medianFrequencySeconds','maximumFrequencySeconds', ...
        'evaluatedCandidates','referenceFallbackFrequencies'});
    ghcnh_log(logFile,label,sprintf( ...
        ['solve-aware trace audit complete in %.2f min; evaluated=%d; ' ...
         'reference fallbacks=%d'],totalSeconds/60,sum(evaluated,'all'), ...
        sum(selectedSubset==referenceJ)));

    function report_progress(message)
        completed=completed+1;
        if mod(completed,5)==0 || completed==nFrequency
            elapsed=toc(auditTimer);
            remaining=elapsed/max(completed,1)*(nFrequency-completed);
            ghcnh_log(logFile,label,sprintf( ...
                ['%d/%d complete; last k=%d J=%d candidates=%d %.2fs; ' ...
                 'estimated remaining %.1f min'],completed,nFrequency, ...
                message(2)-1,message(3),message(4),message(5),remaining/60));
        end
    end
end

function [chosen,localValues,localPass,localEvaluated,elapsed]= ...
        audit_frequency(k,K,Phi,weights,noiseVariance,candidateJ, ...
        referenceM,referenceJ,epsilonKlPerSite,n)
    timer=tic;
    factor=1; if k==1 || k==K, factor=0.5; end
    P=weighted_covariance(Phi(:,1:referenceM),weights(1:referenceM,k), ...
        noiseVariance);
    Q=noiseVariance*eye(n);
    previousM=0; chosen=referenceJ;
    nCandidate=numel(candidateJ);
    localValues=nan(1,nCandidate); localPass=false(1,nCandidate);
    localEvaluated=false(1,nCandidate);
    for q=1:nCandidate
        J=candidateJ(q); M=J^2;
        index=(previousM+1):M;
        if ~isempty(index)
            block=Phi(:,index);
            Q=Q+real(block*(weights(index,k).*block'));
            Q=real((Q+Q')/2);
        end
        R=real((P-Q+P'-Q')/2);
        [LQ,jitterUsed]=eq9_chol_with_jitter(Q);
        if jitterUsed>0
            error(['The solve-aware certificate requires an unjittered ' ...
                'candidate covariance; Cholesky used jitter %.3g.'],jitterUsed);
        end
        traceTerm=max(0,real(trace(LQ'\(LQ\R))));
        value=max(0,factor*(traceTerm-log1p(traceTerm))/n);
        localValues(q)=value; localPass(q)=value<=epsilonKlPerSite;
        localEvaluated(q)=true;
        if localPass(q)
            chosen=J;
            break
        end
        previousM=M;
    end
    elapsed=toc(timer);
end

function K=weighted_covariance(Phi,weight,noise)
    K=real(Phi*(weight(:).*Phi'))+noise*eye(size(Phi,1));
    K=real((K+K')/2);
end
