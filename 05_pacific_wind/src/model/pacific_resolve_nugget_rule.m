function rule = pacific_resolve_nugget_rule(cfg)
%PACIFIC_RESOLVE_NUGGET_RULE Read the training-only nugget decision.

% If no completed decision is supplied, the statistically conservative
% default is the original free nugget with its small numerical lower bound.

    rule.mode = "free_current";
    rule.residualVarianceFraction = NaN;
    rule.source = "default_current_model";
    if ~isfield(cfg, 'nugget') || ~cfg.nugget.useDecisionFile
        return;
    end
    path = cfg.nugget.decisionFile;
    if ~isfile(path)
        warning(['Nugget decision file was not found; retaining the ', ...
            'current free nugget: %s'], path);
        return;
    end
    decision = readtable(path, 'TextType', 'string');
    if height(decision) ~= 1
        error('Nugget decision file must contain exactly one row.');
    end
    if decision.recommendedMode == ...
            "larger_scale_adaptive_floor_supported"
        rho = decision.recommendedResidualVarianceFraction;
        if ~(isfinite(rho) && rho > 0 && rho < 1)
            error('Recommended residual-variance fraction is invalid.');
        end
        rule.mode = "scale_adaptive_floor";
        rule.residualVarianceFraction = rho;
        rule.source = "nested_validation_decision";
    else
        rule.source = "nested_validation_retained_current_floor";
    end
end
