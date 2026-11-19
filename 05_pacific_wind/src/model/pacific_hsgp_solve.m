function solution = pacific_hsgp_solve(factor, rightHandSide)
%PACIFIC_HSGP_SOLVE Apply a factor returned by pacific_hsgp_factor.

    if factor.mode == "feature_qr"
        solution = pacific_rff_solve(factor.core, rightHandSide);
    elseif factor.mode == "site_cholesky"
        solution = factor.L' \ (factor.L \ rightHandSide);
    else
        error('Unknown Pacific HSGP factorization mode.');
    end
end
