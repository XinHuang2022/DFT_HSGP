function fit=evaluate_adaptive_hsgp_at_theta(Y,xy,design,cfg,master, ...
        selectedJ,theta,A)
%EVALUATE_ADAPTIVE_HSGP_AT_THETA Recover profiled mean quantities.
% This evaluates the fitted covariance parameters at a supplied rank profile
% without re-optimizing them.  It is used only to reconstruct capped runs
% under the same rank profile used for their final parameter fit.

    prepared=A.prepare(Y,xy,design);
    prepared.masterPhi=A.features(xy,master);
    selectedJ=selectedJ(:);
    lambdaAll=eq9_fourier_frequencies(prepared.T);
    lambda=lambdaAll(prepared.usedIndices);
    weights=A.weights(theta,lambda,master,cfg);
    p=prepared.nCoefficients; n=prepared.nStations;
    H=zeros(p,p); g=zeros(p,1);
    endpointSet=prepared.endpointIndices(:).';
    for q=1:numel(prepared.usedIndices)
        originalIndex=prepared.usedIndices(q);
        position=positive_frequency_position(originalIndex,prepared.T);
        M=selectedJ(position)^2;
        Phi=prepared.masterPhi(:,1:M);
        factor=A.factor(Phi,weights(1:M,q),theta.noiseVariance,cfg);
        d=prepared.D(originalIndex,:).';
        X=reshape(prepared.Afrequency(originalIndex,:,:),n,p);
        endpoint=ismember(originalIndex,endpointSet);
        if endpoint, d=real(d); X=real(X); end
        solveD=A.solve(factor,d); solveX=A.solve(factor,X);
        if endpoint
            H=H+real(X'*solveX); g=g+real(X'*solveD);
        else
            H=H+2*real(X'*solveX); g=g+2*real(X'*solveD);
        end
    end
    H=real((H+H')/2);
    [LH,~]=eq9_chol_with_jitter(H);
    eta=LH'\(LH\g);
    covarianceEta=LH'\(LH\eye(p));
    fit=struct('theta',theta,'eta',eta,'covarianceEta',covarianceEta, ...
        'selectedJ',selectedJ);
end

function position=positive_frequency_position(index,T)
    K=floor(T/2)+1;
    if index<=K, position=index; else, position=T-index+2; end
end
