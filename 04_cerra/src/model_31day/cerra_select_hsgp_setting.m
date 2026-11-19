function [boundaryFactor, J, passed] = cerra_select_hsgp_setting(audit, cfg)
%CERRA_SELECT_HSGP_SETTING Select the smallest response-free passing basis.

    pass = cerra_logical_column(audit.passesAudit);
    eligible = audit(pass, :);
    passed = ~isempty(eligible);
    if passed
        eligible = sortrows(eligible, {'M','boundaryFactor','J'}, ...
            {'ascend','ascend','ascend'});
        boundaryFactor = eligible.boundaryFactor(1);
        J = eligible.J(1);
        return;
    end
    if ~cfg.hsgp.allowFailedAuditContinuation
        error(['No (boundary factor, J) pair passed Stage 1. This is a ' ...
            'numerical result, not a missing value. Inspect the audit or set ' ...
            'cfg.hsgp.allowFailedAuditContinuation=true for an explicitly ' ...
            'labelled exploratory continuation.']);
    end
    audit.score = audit.p90LatentCovarianceError + ...
        abs(audit.nllDifferencePerObservation);
    audit = sortrows(audit, {'score','M','boundaryFactor'}, ...
        {'ascend','ascend','ascend'});
    boundaryFactor = audit.boundaryFactor(1);
    J = audit.J(1);
end

