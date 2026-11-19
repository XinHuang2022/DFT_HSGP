function [selectedJ, audit, reference] = adaptive_rank_audit(theta, ...
        xy, cfg, master, C, A, logFile, label)
%ADAPTIVE_RANK_AUDIT Stabilize the internal reference and select J_k.

    T = cfg.data.expectedTimePoints;
    K = floor(T/2)+1;
    lambdaAll = eq9_fourier_frequencies(T);
    lambda = lambdaAll(1:K);
    Phi = A.features(xy, master);
    weights = A.weights(theta, lambda, master, cfg);
    n = size(xy,1);

    referenceRows = cell(numel(C.referenceJ)-1,1);
    selectedReferenceJ = C.referenceJ(end);
    referenceStable = false;
    for pair = 1:(numel(C.referenceJ)-1)
        lowJ = C.referenceJ(pair); highJ = C.referenceJ(pair+1);
        maximum = 0; meanValue = 0;
        for k = 1:K
            P = covariance_at(Phi, weights(:,k), lowJ, ...
                theta.noiseVariance, cfg, A);
            Q = covariance_at(Phi, weights(:,k), highJ, ...
                theta.noiseVariance, cfg, A);
            factor = endpoint_factor(k,K);
            value = max(gaussian_kl(P,Q,factor), ...
                gaussian_kl(Q,P,factor)) / n;
            maximum = max(maximum,value); meanValue=meanValue+value/K;
        end
        stable = maximum <= C.epsilonReferencePerSite;
        referenceRows{pair}=table(string(label),lowJ,highJ,meanValue, ...
            maximum,stable,'VariableNames',{'auditLabel','lowerJ', ...
            'higherJ','meanBidirectionalKlPerSite', ...
            'maximumBidirectionalKlPerSite','stable'});
        log_progress(logFile,label,sprintf(['reference J=%d vs %d: ', ...
            'max KL/site %.4g (%s)'],lowJ,highJ,maximum, ...
            string(stable)));
        if stable
            selectedReferenceJ=lowJ; referenceStable=true; break;
        end
    end
    reference = vertcat(referenceRows{1:pair});
    if ~referenceStable
        log_progress(logFile,label,sprintf(['reference ladder did not ', ...
            'stabilize; using highest J=%d'],selectedReferenceJ));
    end

    candidates = C.candidateJ(C.candidateJ < selectedReferenceJ);
    if isempty(candidates), error('No candidate rank is below reference rank.'); end
    selectedJ = repmat(selectedReferenceJ,K,1);
    rows = cell(K*numel(candidates),1); cursor=0;
    unresolved = false(K,1);
    for k = 1:K
        P = covariance_at(Phi,weights(:,k),selectedReferenceJ, ...
            theta.noiseVariance,cfg,A);
        found=false;
        for q = 1:numel(candidates)
            J = candidates(q);
            Q = covariance_at(Phi,weights(:,k),J, ...
                theta.noiseVariance,cfg,A);
            value = gaussian_kl(P,Q,endpoint_factor(k,K))/n;
            cursor=cursor+1;
            rows{cursor}=table(string(label),k-1,lambda(k),J,J^2, ...
                selectedReferenceJ,selectedReferenceJ^2,value, ...
                value<=C.epsilonKlPerSite,'VariableNames', ...
                {'auditLabel','frequencyIndex','lambda','candidateJ', ...
                'candidateM','referenceJ','referenceM','klPerSite','passes'});
            if ~found && value<=C.epsilonKlPerSite
                selectedJ(k)=J; found=true;
            end
        end
        unresolved(k)=~found;
        if mod(k,C.progressEveryFrequencies)==0 || k==K
            log_progress(logFile,label,sprintf(['rank audit %d/%d; ', ...
                'current J=%d; reference fallbacks=%d'],k,K, ...
                selectedJ(k),sum(unresolved(1:k))));
        end
    end
    audit=vertcat(rows{1:cursor});
    reference.selectedReferenceJ=repmat(selectedReferenceJ,height(reference),1);
    reference.referenceStable=repmat(referenceStable,height(reference),1);
    reference.nReferenceFallbackFrequencies= ...
        repmat(sum(unresolved),height(reference),1);
end

function K = covariance_at(Phi,weight,J,noise,cfg,A)
    M=J^2;
    K=A.covariance(Phi(:,1:M),Phi(:,1:M),weight(1:M),noise,true);
    K=real((K+K')/2);
end

function value = gaussian_kl(P,Q,factor)
    [LP,~]=eq9_chol_with_jitter(P); [LQ,~]=eq9_chol_with_jitter(Q);
    logP=2*sum(log(diag(LP))); logQ=2*sum(log(diag(LQ)));
    traceTerm=real(trace(LQ'\(LQ\P)));
    value=max(0,factor*(logQ-logP+traceTerm-size(P,1)));
end

function value=endpoint_factor(k,K)
    if k==1 || k==K, value=0.5; else, value=1; end
end

function log_progress(file,label,message)
    stamp=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    line=sprintf('%s\t%s\t%s\n',stamp,char(label),message);
    fid=fopen(file,'a'); if fid>=0, fprintf(fid,'%s',line); fclose(fid); end
    fprintf('[%s] %s: %s\n',stamp,char(label),message);
end
