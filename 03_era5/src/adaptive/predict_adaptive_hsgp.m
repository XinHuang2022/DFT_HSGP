function prediction = predict_adaptive_hsgp(Yobserved, xyObserved, ...
        xyTarget, observedDesign, targetDesign, fit, cfg, master, A)
%PREDICT_ADAPTIVE_HSGP Whole-record reconstruction with selected J_k.

    [T,nObserved]=size(Yobserved); nTarget=size(xyTarget,1);
    fittedMean=A.evaluateMean(observedDesign,fit.eta);
    D=fft(Yobserved-fittedMean,[],1)/sqrt(T);
    lambda=eq9_fourier_frequencies(T);
    PhiO=A.features(xyObserved,master); PhiT=A.features(xyTarget,master);
    weights=A.weights(fit.theta,lambda,master,cfg);
    p=observedDesign.nCoefficients;
    Dprediction=complex(zeros(T,nTarget));
    sensitivityFrequency=complex(zeros(T,nTarget,p));
    conditionalVarianceFrequency=zeros(T,nTarget);
    maximumJitter=0; minimumRcond=Inf;
    for k=1:T
        position=positive_position(k,T); J=fit.selectedJ(position); M=J^2;
        factor=A.factor(PhiO(:,1:M),weights(1:M,k), ...
            fit.theta.noiseVariance,cfg);
        solveD=A.solve(factor,D(k,:).');
        Kto=A.covariance(PhiT(:,1:M),PhiO(:,1:M), ...
            weights(1:M,k),0,false);
        Dprediction(k,:)=(Kto*solveD).';
        Ao=reshape(observedDesign.Afrequency(k,:,:),nObserved,p);
        At=reshape(targetDesign.Afrequency(k,:,:),nTarget,p);
        solveA=A.solve(factor,Ao);
        sensitivity=At-Kto*solveA;
        sensitivityFrequency(k,:,:)=reshape(sensitivity,1,nTarget,p);
        solveCross=A.solve(factor,Kto');
        Ktt=A.covariance(PhiT(:,1:M),PhiT(:,1:M), ...
            weights(1:M,k),0,true);
        conditional=real((Ktt-Kto*solveCross));
        conditional=real((conditional+conditional')/2);
        conditionalVarianceFrequency(k,:)=max(0,diag(conditional)).';
        maximumJitter=max(maximumJitter,factor.jitter);
        minimumRcond=min(minimumRcond,factor.reciprocalCondition);
    end
    residualComplex=ifft(Dprediction,[],1)*sqrt(T);
    residualMean=real(residualComplex);
    modeledMean=A.evaluateMean(targetDesign,fit.eta);
    conditionalVariance=repmat(mean(conditionalVarianceFrequency,1),T,1);
    sensitivityTime=real(ifft(sensitivityFrequency,[],1)*sqrt(T));
    covarianceEta=real((fit.covarianceEta+fit.covarianceEta')/2);
    fixedEffectVariance=zeros(T,nTarget);
    for j=1:nTarget
        S=reshape(sensitivityTime(:,j,:),T,p);
        fixedEffectVariance(:,j)=max(0,real(sum((S*covarianceEta).*S,2)));
    end
    latentVariance=max(cfg.hsgp.minimumPredictiveVariance, ...
        conditionalVariance+fixedEffectVariance);
    observationVariance=max(cfg.hsgp.minimumPredictiveVariance, ...
        latentVariance+fit.theta.noiseVariance);
    prediction.mean=modeledMean+residualMean;
    prediction.meanOnly=modeledMean;
    prediction.residualMean=residualMean;
    prediction.latentVariance=latentVariance;
    prediction.observationVariance=observationVariance;
    prediction.fixedEffectVariance=fixedEffectVariance;
    prediction.conditionalProcessVariance=conditionalVariance;
    prediction.maximumFeatureJitter=maximumJitter;
    prediction.minimumFeatureReciprocalCondition=minimumRcond;
    prediction.maximumImaginaryPrediction=max(abs(imag(residualComplex)),[],'all');
end

function position=positive_position(k,T)
    K=floor(T/2)+1;
    if k<=K, position=k; else, position=T-k+2; end
end
