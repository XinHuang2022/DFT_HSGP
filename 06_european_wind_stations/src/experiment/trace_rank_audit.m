function [selectedJ,audit,spotChecks] = trace_rank_audit(theta,xy,cfg, ...
        master,C,A,logFile,label)
%TRACE_RANK_AUDIT Algorithm-2 rank selection using the trace certificate.

    T=cfg.data.expectedTimePoints; K=floor(T/2)+1;
    lambdaAll=eq9_fourier_frequencies(T); lambda=lambdaAll(1:K);
    Phi=A.features(xy,master);
    weights=A.weights(theta,lambda,master,cfg);
    n=size(Phi,1); referenceM=C.referenceJ^2;
    normSquared=sum(abs(Phi(:,1:referenceM)).^2,1)';
    candidateJ=C.candidateJ(C.candidateJ<C.referenceJ);
    selectedJ=repmat(C.referenceJ,K,1);
    rows=cell(K*numel(candidateJ),1); cursor=0;
    for k=1:K
        contribution=weights(1:referenceM,k).*normSquared;
        reverseTail=flipud(cumsum(flipud(contribution)));
        factor=1;
        if k==1 || k==K, factor=0.5; end
        found=false;
        for q=1:numel(candidateJ)
            J=candidateJ(q); M=J^2;
            if M<referenceM, tail=reverseTail(M+1); else, tail=0; end
            bound=factor*tail/(n*theta.noiseVariance);
            cursor=cursor+1;
            rows{cursor}=table(string(label),k-1,lambda(k),J,M, ...
                C.referenceJ,referenceM,bound,bound<=C.epsilonKlPerSite, ...
                'VariableNames',{'auditLabel','frequencyIndex','lambda', ...
                'candidateJ','candidateM','referenceJ','referenceM', ...
                'traceBoundPerSite','passes'});
            if ~found && bound<=C.epsilonKlPerSite
                selectedJ(k)=J; found=true;
            end
        end
        if mod(k,C.progressEveryFrequencies)==0 || k==K
            ghcnh_log(logFile,label,sprintf( ...
                'trace screen %d/%d; selected J=%d; reference fallbacks=%d', ...
                k,K,selectedJ(k),sum(selectedJ(1:k)==C.referenceJ)));
        end
    end
    audit=vertcat(rows{1:cursor});
    spotChecks=exact_spot_checks(Phi,weights,selectedJ,theta.noiseVariance, ...
        C.referenceJ,lambda,n);
end

function checks=exact_spot_checks(Phi,weights,selectedJ,noise,referenceJ,lambda,n)
    K=numel(lambda);
    frequencyIndex=unique([1,round(K/4),round(K/2),round(3*K/4),K]);
    rows=cell(numel(frequencyIndex),1);
    referenceM=referenceJ^2;
    normSquared=sum(abs(Phi(:,1:referenceM)).^2,1)';
    for r=1:numel(frequencyIndex)
        k=frequencyIndex(r); J=selectedJ(k); M=J^2;
        P=real(Phi(:,1:referenceM)*(weights(1:referenceM,k).* ...
            Phi(:,1:referenceM)'))+noise*eye(n);
        Q=real(Phi(:,1:M)*(weights(1:M,k).*Phi(:,1:M)'))+noise*eye(n);
        P=(P+P')/2; Q=(Q+Q')/2;
        LP=chol(P,'lower'); LQ=chol(Q,'lower');
        factor=1; if k==1 || k==K, factor=0.5; end
        exact=max(0,factor*(2*sum(log(diag(LQ)))- ...
            2*sum(log(diag(LP)))+trace(LQ'\(LQ\P))-n)/n);
        tail=sum(weights(M+1:referenceM,k).*normSquared(M+1:referenceM));
        bound=factor*tail/(n*noise);
        rows{r}=table(k-1,lambda(k),J,M,exact,bound,bound-exact, ...
            'VariableNames',{'frequencyIndex','lambda','selectedJ', ...
            'selectedM','exactKlPerSite','traceBoundPerSite','boundGap'});
    end
    checks=vertcat(rows{:});
end
