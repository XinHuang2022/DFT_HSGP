function reconstruct_relaxed_algorithm1_capped_folds()
%RECONSTRUCT_RELAXED_ALGORITHM1_CAPPED_FOLDS Correct capped post-fit ranks.
% Folds 3--5 reached the four-fit cap.  Their original prediction step used
% the post-audit requested ranks rather than the ranks used in fit 4.  This
% function re-evaluates the profiled mean at the saved covariance parameters
% and reconstructs under the actual fit-4 ranks, without re-optimization.

    C=ghcnh_trace_config(); [cfg,subset,folds,A,E]=ghcnh_expanded_setup();
    for foldID=3:5
        targetIndex=folds(foldID).targetIndex(:);
        observedIndex=setdiff((1:size(subset.U,2))',targetIndex,'stable');
        Yo=subset.U(:,observedIndex); Yt=subset.U(:,targetIndex);
        [observed,target,domainXy]=fold_geometry(subset,observedIndex, ...
            targetIndex,cfg,A);
        master=adaptive_basis_bank(domainXy,cfg.hsgp.boundaryFactor, ...
            C.referenceJ,A,cfg);
        methodDir=fullfile(E.resultsDir,'relaxed_cycle_sensitivity', ...
            sprintf('fold_%02d',foldID),'maxcycles_4_fraction_0100','exact_kl');
        parameterTable=readtable(fullfile(methodDir,'parameters.csv'));
        thetaFields={'signalVariance','alpha','beta','logSpatialScaleContrast', ...
            'epsilon','noiseVariance','nu'};
        theta=struct();
        for j=1:numel(thetaFields)
            theta.(thetaFields{j})=parameterTable.(thetaFields{j})(1);
        end
        fittedProfile=readtable(fullfile(methodDir,'cycle_3_rank_profile.csv'));
        selectedJ=fittedProfile.selectedJ;
        fit=evaluate_adaptive_hsgp_at_theta(Yo,observed.xyCovariance, ...
            observed.design,cfg,master,selectedJ,theta,A);
        prediction=predict_adaptive_hsgp(Yo,observed.xyCovariance, ...
            target.xyCovariance,observed.design,target.design,fit,cfg,master,A);
        [metrics,siteMetrics,predictions]=score_ghcnh_predictions( ...
            Yt,prediction,subset,targetIndex,foldID,cfg, ...
            "relaxed_exact_kl_fitted_rank");
        outputDir=fullfile(methodDir,'fitted_rank_reconstruction');
        if ~isfolder(outputDir), mkdir(outputDir); end
        writetable(metrics,fullfile(outputDir,'metrics.csv'));
        writetable(siteMetrics,fullfile(outputDir,'site_metrics.csv'));
        writetable(predictions,fullfile(outputDir,'predictions.csv'));
        writetable(fittedProfile,fullfile(outputDir,'fitted_rank_profile.csv'));
        fprintf('Corrected Algorithm 1 fold %d RMSE %.6f.\n', ...
            foldID,metrics.rmse);
    end
end

function [observed,target,domainXy]=fold_geometry(subset,observedIndex, ...
        targetIndex,cfg,A)
    allCenter=mean(subset.xyKm,1);
    domainXy=(subset.xyKm-allCenter)/cfg.spatial.coordinateScaleKm;
    observed.xyCovariance=domainXy(observedIndex,:);
    target.xyCovariance=domainXy(targetIndex,:);
    center=mean(subset.xyKm(observedIndex,:),1);
    scale=std(subset.xyKm(observedIndex,:),0,1); scale(scale==0)=1;
    observed.design=A.makeDesign(size(subset.U,1), ...
        (subset.xyKm(observedIndex,:)-center)./scale);
    target.design=A.makeDesign(size(subset.U,1), ...
        (subset.xyKm(targetIndex,:)-center)./scale);
end
